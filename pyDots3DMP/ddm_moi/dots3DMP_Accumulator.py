# %% ----------------------------------------------------------------

import numpy as np
import pandas as pd
from typing import Union, Optional
from sklearn.base import BaseEstimator
from sklearn.utils.validation import check_X_y, check_array, check_is_fitted
from scipy.signal import convolve
from scipy.stats import norm, truncnorm, skewnorm

from pybads import BADS
from scipy.optimize import minimize

from behavior.preprocessing import (
    dots3DMP_create_trial_list, dots3DMP_create_conditions,
    data_cleanup, format_onetargconf,
    )

import itertools
from functools import partial
from copy import deepcopy, copy

import logging

from Accumulator import Accumulator

# TODO list

# - make commits!
# - get optimization running
# - get visualizations of model predictions
# - implement simulate method
# - bstract common parts of predict methods
# - cleanup RT part
# - add logging 
# - add mlflow

# %% ----------------------------------------------------------------

logger = logging.getLogger()
logging.basicConfig(level=logging.INFO)

def main():
    accum_kw = {'grid': np.arange(-3, 0, 0.05),
                'tvec': np.arange(0, 2, 0.05),
                'kmult': [0.3, 0.3],
                'bound': [0.5, 0.5, 0.5],
                'non_dec_time': [0.3],
                'wager_thr': [0.8, 0.8, 0.8],
                'wager_alpha': [0]}
    accum = dots3dmpAccumulator(**accum_kw)

    datafilepath = "/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/dataStructs/lucio_20220512-20230606.csv"
    data = data_cleanup(datafilepath)  # this is a hack function to quickly load and clean the data, could be improved/generalized with options
    data = format_onetargconf(data, remove_one_targ=True)

    X = data[['modality', 'coherence', 'delta', 'heading']]
    y = data[['choice', 'PDW', 'RT']].astype({'choice': 'int', 'PDW': 'int'})

    accum.fit(X, y, fixed_params=['kmult', 'wager_alpha'])
    #accum.predict_proba(X)

# %% ----------------------------------------------------------------


class AccumulatorModel2D(BaseEstimator):

    PARAM_NAMES = ('kmult', 'bound', 'non_dec_time', 'wager_thr', 'wager_alpha') # class attribute

    def __init__(self, grid, tvec, kmult, bound, non_dec_time=0.3, 
                 wager_thr=1.0, wager_alpha=0.0, wager_axis=None):
        self.grid = grid
        self.tvec = tvec
        self.kmult = kmult
        self.bound = bound
        self.non_dec_time = non_dec_time
        self.wager_thr = wager_thr
        self.wager_alpha = wager_alpha
        self.wager_axis = None # None means both axes i.e. time and evidence (log odds)
        
    def fit(self, X, y, fixed_params=None):
        pass
        
    def predict(self, X):
        pass

    def simulate(self, X):
        pass
    
    # TODO define other generic AccumulatorModel methods


class dots3dmpAccumulator(AccumulatorModel2D):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)  # Inherit parent class parameters

    # https://stackoverflow.com/questions/51430484/how-to-subclass-a-vectorizer-in-scikit-learn-without-repeating-all-parameters-in
    def get_params(self, deep=True):
        # Hack to make get_params return base class params...
        params = super().get_params(deep)
        cp = copy(self)
        cp.__class__ = AccumulatorModel2D
        params.update(AccumulatorModel2D.get_params(cp, deep))
        return params

    def _get_fit_params(self):

        # convert kwargs to dict,
        # subset out fittable parameters and fixed parameters
        # make fittable parameters, return as array, with lens
        self.model_params = {k: getattr(self, k) for k in self.PARAM_NAMES}
        self.fit_param_names = self.model_params.keys()
        if self.fixed_params:
            self.fixed_params = {k: self.model_params[k] for k in self.fixed_params}
            self.fit_param_names = [k for k in self.fit_param_names if k not in self.fixed_params.keys()]

        self.fitted_params_ = {k: self.model_params[k] for k in self.fit_param_names}
        params_list = list({k: v for k,v in self.model_params.items() if k in self.fit_param_names}.values())

        return params_list    
            
    def fit(self, X: pd.DataFrame, y: pd.DataFrame, fixed_params=None):
        
        self.n_features_in_ = len(X.columns)
        self.fixed_params = fixed_params

        params_list  = self._get_fit_params()
        params_array = np.array(list(itertools.chain(*params_list)))

        # pass fixed model_params and data as fixed inputs to objective function
        optim_fcn_part = lambda params: self._objective_fcn(params, X, y)

        lb = params_array * 0.3
        ub = params_array * 2
        plb = params_array * 0.5
        pub = params_array * 1.5

        bads_bounds = (lb, ub, plb, pub)
        bads = BADS(optim_fcn_part, params_array, *bads_bounds)
        result = bads.optimize()

        self.fitted_params_ = dict(zip(fit_param_names, result.x))

        return self

    def _objective_fcn(self, params_array, X, y):

        params_list = self._get_fit_params()
        end_inds = list(itertools.accumulate(map(len, params_list)))

        # from params_array and fixed_params, reconstruct full params_dict
        start_ind = 0
        for p, param in enumerate(self.fit_param_names):
            if p > 0:
                start_ind = end_inds[p-1]
            self.fitted_params_[param] = params_array[start_ind:end_inds[p]].tolist()

        self.params_ = self.fitted_params_.copy()
        if self.fixed_params:
            self.params_ = self.params_ | self.fixed_params
        # self.params = {k: self.params[k] for k in self.PARAM_NAMES}

        y_pred = self.predict_proba(X, y['RT'].to_numpy())

        y_obs = np.array(y[['choice', 'PDW', 'RT']])
        y_pred = np.array(y_pred[['choice', 'PDW', 'RT']])

        log_lik_choice = self.log_lik_bin(y_obs[:,0], y_pred[:,0])
        log_lik_pdw    = self.log_lik_bin(y_obs[:,1], y_pred[:,1])
        log_lik_rt     = self.log_lik_cont(y_obs[:,2], y_pred[:,2], self.tvec)

        self.log_lik_ = {'choice': log_lik_choice, 'pdw': log_lik_pdw, 'rt': log_lik_rt}

        self.neg_llh_ = -sum(self.log_lik_.values())
        # neg_llh = -sum([self.log_lik[v]*w for v, w in zip(outputs, llh_scaling) if v in self.log_lik])

        logger.info('Total loss:\t%.2f', self.neg_llh_)
        logger.info('Individual log likelihoods: %s', {key: round(val, 2) for key,val in self.log_lik_.items()})

        return self.neg_llh_

    def predict(self, X):
        # run new data through (fitted) model parameters, generate sampled predictions
        # This is going to have redundancies with predict_proba, so make it DRY
        pass

    @staticmethod
    def _handle_kmult(kmult, cohs, k_scale=1):
        """break down kmult into (kves, kvis)"""
        kmult *= k_scale
        kves = kmult[0]
        if len(kmult) == 1:
            kvis = kmult * cohs
            kves = np.mean(kvis)  # ves lies between vis cohs
        if len(kmult) == 2:
            kvis = kmult[1] * cohs
        else:
            kvis = kmult[1:]  # 3 independent kmults

        return kves, kvis
            
    def predict_proba(self, X, observed_RTs):
        """make a probabilistic prediction, given model params"""

        mods = np.unique(X['modality'])
        cohs = np.unique(X['coherence'])
        deltas = np.unique(X['delta'])
        hdgs = np.unique(X['heading'])

        cond_groups = X[['modality','coherence','delta']].drop_duplicates(
            ).sort_values(by=['modality','coherence','delta']).reset_index(drop=True)

        kves, kvis = self._handle_kmult(self.params_['kmult'], cohs.T)     

        # TODO Fix this, need to use separate cohs/
        # maybe assign k and t_vals to rows on cond_groups, and loop over that
        # so we make sure that we are only using the relevant ones for mods

        # also this should maybe also be stored in self

        b_ves, b_vis = np.ones_like(self.tvec), np.ones_like(self.tvec)
        # divide by len for time course of sensitivity
        b_ves /= len(b_ves)
        b_vis /= len(b_vis)

        b_vals = [b_ves, b_vis, np.vstack((b_ves, b_vis)).T]
        k_vals = [kves, np.mean(kvis), [kves, np.mean(kvis)]]

        wager_maps, wager_above_threshold = self.get_log_odds(X, b_vals, k_vals)

        predictions = np.full((X.shape[0], 3), fill_value=np.nan)
        predictions = pd.DataFrame(predictions, columns=['choice', 'PDW', 'RT'])

        for m, mod in enumerate(mods):

            bound = self.params_['bound']
            if len(bound) == len(mods):
                bound = bound[m]

            non_dec_time = self.params_['non_dec_time']
            if len(non_dec_time) == len(mods):
                non_dec_time = non_dec_time[m]

            for c, coh in enumerate(cohs):
                k_vals = [kves, kvis[c], [kves, kvis[c]]]
                for d, delta in enumerate(deltas):

                    # non-valid conditions
                    if (delta != 0 and mod < 3) or (c > 0 and mod == 1):
                        continue       
                        
                    accumulator = Accumulator(grid_vec=self.grid, tvec=self.tvec)
                    accumulator.bound = bound
                    orig_tvec = accumulator.tvec # for later
                    
                    drifts, accumulator.tvec = calc_3dmp_drift_rates(
                        b_vals[m], k_vals[m], self.tvec[-1], hdgs, delta=delta, return_abs=False
                        )
                    drifts_list = np.split(drifts, len(hdgs), axis=1)
                    accumulator.set_drifts(drifts_list, hdgs)
                    accumulator.dist(return_pdf=True)

                    for h, hdg in enumerate(hdgs):
                        trial_index = (X['modality'] == mod) & (X['coherence'] == coh) & \
                                    (X['heading'] == hdg) & (X['delta'] == delta)
                        trial_index = trial_index.to_numpy()
                        if trial_index.sum() == 0:
                            continue

                        p_right = accumulator.p_corr_[h].item()
                        p_right = np.clip(p_right, 1e-100, 1-1e-100)
                        p_choice = np.array([p_right, 1 - p_right])
                        
                        # ====== CHOICE ======
                        predictions.loc[trial_index, 'choice'] = p_right
                        
                        # ====== WAGER ======
                        # select pdf for losing race, given correct or incorrect
                        pxt_up = np.squeeze(accumulator.up_lose_pdf_[h, :, :])
                        pxt_lo = np.squeeze(accumulator.lo_lose_pdf_[h, :, :])
                        total_p = np.sum(pxt_up + pxt_lo)
                        pxt_up /= total_p
                        pxt_lo /= total_p

                        p_choice_and_wager = np.array(
                            [[np.sum(pxt_up[wager_above_threshold[m]]),     # pRight+High
                                np.sum(pxt_up[~wager_above_threshold[m]])],   # pRight+Low
                                [np.sum(pxt_lo[wager_above_threshold[m]]),     # pLeft+High
                                np.sum(pxt_lo[~wager_above_threshold[m]])]],  # pLeft+Low
                        )

                        # calculate p_wager using Bayes rule, then factor in base rate of low bets ("alpha")
                        p_choice_given_wager, p_wager = _margconds_from_intersection(p_choice_and_wager, p_choice)
                        p_wager += np.array([-self.params_['wager_alpha'][0], self.params_['wager_alpha'][0]]) * p_wager[0]
                        p_wager = np.clip(p_wager, 1e-100, 1-1e-100)

                        predictions.loc[trial_index, 'PDW'] = p_wager[0]

                        # to avoid log(0) issues when doing log-likelihoods, replace zeros and ones
                        predictions.loc[:, ['choice', 'PDW']] = predictions.loc[:, ['choice', 'PDW']].replace(to_replace=0, value=1e-10)
                        predictions.loc[:, ['choice', 'PDW']] = predictions.loc[:, ['choice', 'PDW']].replace(to_replace=1, value=1 - 1e-10)
                        
                        # ====== RT ======

                        # first convolve model RT distribution with non-decision time
                        ndt_dist = norm.pdf(orig_tvec, loc=non_dec_time, scale=0.001) #scale=self.params_['sigma_ndt'])
                        rt_dist = np.squeeze(accumulator.rt_dist_[h, :])

                        rt_dist = np.clip(rt_dist, np.finfo(np.float64).eps, a_max=None)
                        rt_dist = convolve(rt_dist, ndt_dist / ndt_dist.sum())
                        rt_dist = rt_dist[:len(accumulator.tvec)]
                        rt_dist /= rt_dist.sum()  # renormalize

                        # NOTE 12-2023 I think we need to sample from the orig_tvec here, not the changed tvec!

                        # hmm, we need access to original data here even though its a pure model prediction
                        actual_rts = observed_RTs[trial_index] + 0.3 # TODO fix this hard-coded value
                        dist_inds = [np.argmin(np.abs(orig_tvec - rt)) for rt in actual_rts]
                        rt_dist /= rt_dist.max()  # rescale to 0-1, makes likelihood total closer in magnitude to choice and PDW                        
                        predictions.loc[trial_index, 'RT'] = rt_dist[dist_inds]
                        
        return predictions
        

    def get_log_odds(self, X, b_vals, k_vals):

        mods = np.unique(X['modality'])
        hdgs = np.unique(X['heading'])

        cond_groups = X['modality'].drop_duplicates().sort_values().reset_index(drop=True)

        wager_maps = []
        thetas = self.params_['wager_thr']
        
        for m, mod in enumerate(mods):

            bound = self.params_['bound']
            if len(bound) == len(mods):
                bound = bound[m]

            accumulator = Accumulator(grid_vec=self.grid, tvec=self.tvec)
            accumulator.bound = bound
            abs_drifts, accumulator.tvec = calc_3dmp_drift_rates(
                b_vals[m], k_vals[m], self.tvec[-1], hdgs, delta=0, return_abs=True
                )
            drifts_list = np.split(abs_drifts, abs_drifts.shape[1], axis=1)

            # run the method of images - diffusion to bound to extract pdfs, cdfs, and LPO
            accumulator.set_drifts(drifts_list, hdgs[hdgs >= 0])
            accumulator.dist(return_pdf=True)

            if self.wager_axis is None:
                accumulator.log_posterior_odds()
            else:
                raise NotImplementedError('alternatives to log odds not yet implemented')

            wager_maps.append(accumulator.log_odds_)

        wager_above_threshold = [p >= theta for p, theta in zip(wager_maps, thetas)]

        return wager_maps, wager_above_threshold



            

            

            

        

        


    def simulate(self, X):
        check_is_fitted(self, 'fitted_params_')
        # dv simulation for each trial in X, given model params
        
    def loss(self):
        # self.loss = -sum([self.log_lik[v]*w for v, w in zip(outputs, llh_scaling) if v in model_llh])
        pass

    @staticmethod
    def log_lik_bin(y, y_hat):
        return np.sum(np.log(y_hat[y==1])) + np.sum(np.log(1 - y_hat[y==0]))

    @staticmethod
    def log_lik_cont(y, y_hat, t):
        y_inds = np.searchsorted(t, y)
        y_hat /= y_hat.max()
        return np.sum(np.log(y_hat[y_inds]))
        


def calc_3dmp_drift_rates(b_t: np.ndarray, b_k: Union[float, tuple[float, float]],
    tmax: float, hdgs: np.ndarray, delta: Optional[float] = 0.0, 
    return_abs: Optional[bool] = False,
    cue_weights: Optional[tuple[float, float]] = None) -> tuple[np.ndarray, np.ndarray]:
    """
    Calculate drift rates
    if cue_weights is None, optimal cue weights calculated from k's
    """

    sin_uhdgs = np.sin(np.deg2rad(hdgs))
    if return_abs:
        sin_uhdgs = sin_uhdgs[sin_uhdgs >= 0]

    cumul_bt = np.cumsum((b_t**2)/(b_t**2).sum(axis=0), axis=0)

    if isinstance(b_k, (int, float)):
        # only one sensitivity and time-course - ves or vis (logic is the same)
        b_t = b_t.reshape(-1, 1)
        tvec = cumul_bt * tmax
        drifts = np.cumsum(b_t**2 * b_k * sin_uhdgs, axis=0)
        # drifts = b_k * sin_uhdgs   # w/o stim scaling, reduces to this

    elif len(b_k) == 2:
        # two sensitivities/time-courses - combined condition

        b_k = np.array(b_k, dtype=float)
        if cue_weights is None:
            k2 = b_k**2
            cue_weights = np.sqrt(k2 / k2.sum())
        else:
            cue_weights = cue_weights.toarray()

        w_ves, w_vis = cue_weights
        if return_abs:
            drift_ves = np.cumsum(b_t[:, [0]]**2 * b_k[0] * sin_uhdgs, axis=0)
            drift_vis = np.cumsum(b_t[:, [1]]**2 * b_k[1] * sin_uhdgs, axis=0)
        else:
            # +ve delta means ves to the left, vis to the right
            # Drugo eq suggests cumsum each modality separately first, then do the weighted sum     
            drift_ves = np.cumsum(b_t[:, [0]]**2 * b_k[0] * np.sin(np.deg2rad(hdgs - delta / 2)), axis=0)
            drift_vis = np.cumsum(b_t[:, [1]]**2 * b_k[1] * np.sin(np.deg2rad(hdgs + delta / 2)), axis=0)

        # Eq 14
        tvec = w_ves**2 * cumul_bt[:,0] + w_vis**2 * cumul_bt[:,1]
        tvec *= tmax
        drifts = w_ves * drift_ves + w_vis * drift_vis

    drifts /= tvec[:, None]

    return drifts, tvec






    ## % ----------------------------------------------------------------
## PRIVATE FUNCTIONS

def _margconds_from_intersection(prob_ab, prob_a):
    """
    :param prob_ab: joint probability of A and B
    :param prob_a: marginal probability of A
    :return: 
        a_given_b: conditional probability of A given B
        prob_b: marginal probability of B
    """
    prob_a = prob_a.reshape(-1, 1)  # make it 2-D, for element-wise and matrix mults below

    # conditional probability
    b_given_a = (prob_ab / np.sum(prob_ab)) / prob_a

    prob_b = prob_a.T @ b_given_a

    if np.any(prob_b == 0):
        a_given_b = b_given_a * prob_a
    else:
        a_given_b = b_given_a * prob_a / prob_b

    return a_given_b, prob_b.flatten()


def _intersection_from_margconds(a_given_b, prob_a, prob_b):
    """
    Recover intersection of a and b using conditionals and marginals, according to Bayes theorem
    Essentially the inverse of margconds_from_intersection, and the two can be used together e.g.
    to update the intersections after adding a base_rate to prob_b

    :param a_given_b:
    :param prob_a:
    :param prob_b:
    :return:
        prob_ab: joint probability of A and B
        b_given_a
    """

    prob_ab = a_given_b * prob_b
    b_given_a = prob_ab / prob_a

    return prob_ab, b_given_a

# %%

if __name__ == '__main__':
    logger = logging.getLogger(__name__)
    main()