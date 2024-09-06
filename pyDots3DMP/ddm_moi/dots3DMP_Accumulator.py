# %% ----------------------------------------------------------------

import numpy as np
import pandas as pd
from typing import Union, Optional
from sklearn.base import BaseEstimator
from scipy.signal import convolve
from scipy.stats import norm, skewnorm

from pybads import BADS
from scipy.optimize import minimize

from behavior.preprocessing import (
    dots3DMP_create_trial_list, dots3DMP_create_conditions,
    data_cleanup, format_onetargconf,
    )

import itertools
from copy import copy, deepcopy
import logging

from Accumulator import Accumulator

# TODO list

# - get optimization running
# - get visualizations of model predictions
# - implement simulate method
# - abstract common parts of predict and simulate methods together
# - if no wager, drop wager params and handle
# - cleanup RT part
# - add logging 
# - add mlflow



# %% ----------------------------------------------------------------

logger = logging.getLogger()
logging.basicConfig(level=logging.INFO)

def main():
    accum_kw = {'grid': np.arange(-3, 0, 0.05),
                'tvec': np.arange(0, 2, 0.01),
                'kmult': [0.4, 0.7],
                'bound': [0.7, 0.7, 0.7],
                'non_dec_time': [0.3],
                'wager_thr': [.5, .5, .5],
                'wager_alpha': [0]}
    accum = dots3dmpAccumulator(stim_scaling=True, **accum_kw)

    datafilepath = "/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/dataStructs/lucio_20220512-20230606.csv"
    data = data_cleanup(datafilepath)  # this is a hack function to quickly load and clean the data, could be improved/generalized with options
    data = format_onetargconf(data, remove_one_targ=True)
    data = data.reset_index(drop=True)

    X = data[['modality', 'coherence', 'delta', 'heading']]
    y = data[['choice', 'PDW', 'RT']].astype({'choice': 'int', 'PDW': 'int'})
    y['RT'] += 0.3 # to reverse the downshift from motion platform latency

    # accum.fit(X, y, fixed_params=['kmult', 'wager_alpha'])
    y_pred, y_pred_samp = accum.predict(X, y, n_samples=1)

    print(y_pred.head())
    print(y_pred_samp.head())

# %% ----------------------------------------------------------------

class AccumulatorModel2D(BaseEstimator):

    PARAM_NAMES = ('kmult', 'bound', 'non_dec_time', 'wager_thr', 'wager_alpha') # class attribute

    def __init__(self, grid: np.ndarray, tvec: np.ndarray, kmult: list, bound: list, non_dec_time: list=[0.3], 
                 return_wager: bool=True, wager_thr: list=[1], wager_alpha: list=[0.0], wager_axis=None):
        self.grid = grid
        self.tvec = tvec
        self.kmult = kmult
        self.bound = bound
        self.non_dec_time = non_dec_time
        self.return_wager = return_wager
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

    def __init__(self, stim_scaling=False, **kwargs):
        super().__init__(**kwargs)  # Inherit parent class parameters
        self.stim_scaling = stim_scaling

    def get_params(self, deep=True):
        """
        This is a workaround to make the scikit-learn get_params call actually return all parent class parameters
        # https://stackoverflow.com/questions/51430484/how-to-subclass-a-vectorizer-in-scikit-learn-without-repeating-all-parameters-in
        """
        params = super().get_params(deep)
        cp = copy(self)
        cp.__class__ = AccumulatorModel2D
        params.update(AccumulatorModel2D.get_params(cp, deep))
        return params

    def _get_fit_params(self):

        # convert kwargs to dict,
        # subset out fittable parameters and fixed parameters
        # make fittable parameters, return as list
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

        # TODO allow user to specify these instead or with init params
        lb = params_array * 0.3
        ub = params_array * 2
        plb = params_array * 0.5
        pub = params_array * 1.5

        bads_bounds = (lb, ub, plb, pub)
        bads = BADS(optim_fcn_part, params_array, *bads_bounds)
        result = bads.optimize()

        self.fitted_params_ = dict(zip(fit_param_names, result.x))

        return self

    def set_params(self):
        
        params_list  = self._get_fit_params()
        params_array = np.array(list(itertools.chain(*params_list)))
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

        
    def _objective_fcn(self, params_array, X, y):

        self.set_params()
        print(self.params_)

        y_pred = self.predict(X, y['RT'].to_numpy())

        log_lik_choice = log_lik_bin(y['choice'].to_numpy(), y_pred['choice'].to_numpy())
        log_lik_pdw    = log_lik_bin(y['PDW'].to_numpy(), y_pred['PDW'].to_numpy())
        log_lik_rt     = log_lik_cont(y['RT'].to_numpy(), y_pred['RT'].to_numpy(), self.tvec)

        self.log_lik_ = {'choice': log_lik_choice, 'pdw': log_lik_pdw, 'rt': log_lik_rt}
        self.log_lik = {k:v/len(y) for k, v in self.log_lik_.items()}
        
        # self.neg_llh_ = -sum(self.log_lik_.values())
        # self.neg_llh_ = -sum([self.log_lik[v]*w for v, w in zip(outputs, llh_scaling) if v in self.log_lik])

        self.neg_llh_ = self.log_lik_choice 

        logger.info('Total loss:\t%.2f', self.neg_llh_)
        logger.info('Individual log likelihoods: %s', {key: round(val, 2) for key,val in self.log_lik_.items()})

        return self.neg_llh_

    # def predict(self, X, y, nsamples: int=10, seed = None):
    #     # run new data through (fitted) model parameters, generate sampled predictions
    #     # This is going to have redundancies with predict_proba, so make it DRY
            
    #     rng = np.random.RandomState(seed)
    #     predictions = self.predict_proba(X, y)
    #     predictions['choice'] = predictions['choice'].apply(lambda x: rng.binomial(nsamples, x)/nsamples)
    #     predictions['PDW'] = predictions['PDW'].apply(lambda x: rng.binomial(nsamples, x)/nsamples)

    #     return predictions

    def predict(self, X, y, n_samples: int=0, seed=None):
        """make a probabilistic prediction, given model params, or sample from that prediction"""

        rng = np.random.RandomState(seed)

        # this seems messy...and probably not the right approach
        if not hasattr(self, 'neg_llh_'):
            self.fixed_params = None
            self.set_params()
            
        mods = np.unique(X['modality'])
        cohs = np.unique(X['coherence'])
        deltas = np.unique(X['delta'])
        hdgs, hdg_inds = np.unique(X['heading'], return_inverse=True)

        cond_groups = X[['modality','coherence','delta']].drop_duplicates(
            ).sort_values(by=['modality','coherence','delta'])

        kves, kvis = self._handle_kmult(self.params_['kmult'], cohs.T) 
        bound = self._handle_param_mod(self.params_['bound'], mods)  
        non_dec_time = self._handle_param_mod(self.params_['non_dec_time'], mods)  
        thetas = self._handle_param_mod(self.params_['wager_thr'], mods)  

        if not self.stim_scaling:
            b_ves, b_vis = np.ones_like(self.tvec), np.ones_like(self.tvec)
            # divide by len for time course of sensitivity
            b_ves /= len(b_ves)
            b_vis /= len(b_vis)
        elif isinstance(self.stim_scaling, tuple):
            b_ves, b_vis = self.stim_scaling
        else:
            b_ves, b_vis = self.get_stim_urgs(self.tvec)

        b_vals = [b_ves, b_vis, np.vstack((b_ves, b_vis)).T]

        predictions = pd.DataFrame(np.full((X.shape[0], 3), fill_value=np.nan), 
                                   columns=['choice', 'PDW', 'RT'])
        pred_sample = deepcopy(predictions)

        self.wager_maps = []

        k_vals_fixed = [kves, np.mean(kvis), [kves, np.mean(kvis)]]

        for m, mod in enumerate(mods):

            # set accumulators with absolute drifts, for log odds mappings
            
            accumulator = Accumulator(grid_vec=self.grid, tvec=self.tvec)
            accumulator.bound = bound[m]
            abs_drifts, accumulator.tvec = self.calc_3dmp_drift_rates(
                b_vals[m], k_vals_fixed[m], self.tvec[-1], hdgs[hdgs>=0], delta=0,
                )

            # run the method of images - diffusion to bound to extract pdfs, cdfs, and LPO
            accumulator.set_drifts(abs_drifts, hdgs[hdgs>=0])
            accumulator.compute_distrs(return_pdf=True)

            if self.wager_axis is None:
                log_odds_map = accumulator.log_posterior_odds()
                self.wager_maps.append(log_odds_map)
            else:
                raise NotImplementedError('alternatives to log odds not yet implemented')

            wager_above_threshold = [p >= theta for p, theta in zip(self.wager_maps, thetas)]

        # now loop over coherences and deltas with one accumulator each for actual predictions
        
            for c, coh in enumerate(cohs):
                k_vals = [kves, kvis[c], [kves, kvis[c]]]
                
                for d, delta in enumerate(deltas):

                    # non-valid conditions
                    if (delta != 0 and mod < 3) or (c > 0 and mod == 1):
                        continue       
                        
                    accumulator = Accumulator(grid_vec=self.grid, tvec=self.tvec)
                    accumulator.bound = bound[m] # TODO deal with this in initialization?

                    drifts, accumulator.tvec = self.calc_3dmp_drift_rates(
                        b_vals[m], k_vals[m], self.tvec[-1], hdgs, delta=delta,
                        )
                    accumulator.set_drifts(drifts, hdgs) # TODO include this in set_drift_rates, basically wrap it
                    accumulator.compute_distrs(return_pdf=True)

                    trial_index = (X['modality'] == mod) & (X['coherence'] == coh) & (X['delta'] == delta)
                    trial_index = trial_index.to_numpy()

                    # ====== CHOICE ======
                    p_right = np.clip(accumulator.p_corr_.T, 1e-100, 1-1e-100)
                    predictions.loc[trial_index, 'choice'] = p_right[hdg_inds[trial_index]]

                    # ====== WAGER ======

                    # TODO figure out how to vectorize over heading for PDW and RT
                    
                    for h, hdg in enumerate(hdgs):
                        
                        trial_index = (X['modality'] == mod) & (X['coherence'] == coh) & \
                                    (X['heading'] == hdg) & (X['delta'] == delta)
                        trial_index = trial_index.to_numpy()
                        if trial_index.sum() == 0:
                            continue

                        p_choice = np.array([p_right[h], 1 - p_right[h]])
                        
                        # # ====== CHOICE ======                        
                        if n_samples:
                            pred_sample.loc[trial_index, 'choice'] = rng.binomial(n_samples, p_right[h], trial_index.sum()) / n_samples

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

                        # print(f'heading={hdg}, p_wager={p_wager}')
                        predictions.loc[trial_index, 'PDW'] = p_wager[0]

                        if n_samples:
                            pred_sample.loc[trial_index, 'PDW'] = rng.binomial(n_samples, p_wager[0], trial_index.sum()) / n_samples


                        # ====== RT ======

                        # first convolve model RT distribution with non-decision time
                        ndt_dist = norm.pdf(self.tvec, loc=non_dec_time[m], scale=0.001) #scale=self.params_['sigma_ndt'])
                        rt_dist = np.squeeze(accumulator.rt_dist_[h, :])

                        rt_dist = np.clip(rt_dist, np.finfo(np.float64).eps, a_max=None)
                        rt_dist = convolve(rt_dist, ndt_dist / ndt_dist.sum())
                        rt_dist = rt_dist[:len(accumulator.tvec)]
                        rt_dist /= rt_dist.sum()  # renormalize to get posterior

                        # need original data here for RT to get the likelihood in the predictions
                        actual_rts = y['RT'][trial_index] + 0.3
                        dist_inds = [np.argmin(np.abs(self.tvec - rt)) for rt in actual_rts]
                        predictions.loc[trial_index, 'RT'] = rt_dist[dist_inds]

                        if n_samples:
                            pred_sample.loc[trial_index, 'RT'] = rng.choice(self.tvec, trial_index.sum(), replace=True, p=rt_dist) - 0.3

                        
        if n_samples:
            return predictions, pred_sample
        else:
            return predictions
        

    def loss(self):
        # self.loss = -sum([self.log_lik[v]*w for v, w in zip(outputs, llh_scaling) if v in model_llh])
        pass

    def simulate(self, X):
        # dv simulation for each trial in X, given model params
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

    @staticmethod
    # for repeating param for each mod if per-mod is not specificed
    def _handle_param_mod(param, mods):
        if len(param) == 1:
            param *= len(mods)
        return param

    @staticmethod
    def get_stim_urgs(tvec: np.ndarray, pos=None, skew_params=None):

        if pos is None:
            ampl = 0.16

            # pos = norm.cdf(tvec, 0.9, 0.3) * ampl
            if skew_params is None:
                pos = skewnorm.cdf(tvec, 2, 0.8, 0.4) * ampl  # emulate tf
            else:
                pos = skewnorm.cdf(tvec, **skew_params) * ampl

        vel = np.gradient(pos)
        acc = np.gradient(vel)

        vel /= vel.max()
        acc /= acc.max()
        
        return acc, vel
        
    @staticmethod
    def calc_3dmp_drift_rates(b_t: np.ndarray, b_k: Union[float, tuple[float, float]],
        tmax: float, hdgs: np.ndarray, delta: Optional[float] = 0.0, 
        #return_abs: Optional[bool] = False,
        cue_weights: Optional[tuple[float, float]] = None) -> tuple[np.ndarray, np.ndarray]:
        """
        Calculate drift rates
        if cue_weights is None, optimal cue weights calculated from k's
        """

        sin_uhdgs = np.sin(np.deg2rad(hdgs))
        # if return_abs:
        #     sin_uhdgs = sin_uhdgs[sin_uhdgs >= 0]

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
            # if return_abs:
            #     drift_ves = np.cumsum(b_t[:, [0]]**2 * b_k[0] * sin_uhdgs, axis=0)
            #     drift_vis = np.cumsum(b_t[:, [1]]**2 * b_k[1] * sin_uhdgs, axis=0)
            # else:
            
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



# %% ----------------------------------------------------------------
## HELPER FUNCTIONS

def log_lik_bin(y, y_hat):
    return np.sum(y * np.log(y_hat) + (1 - y) * np.log(y_hat))
    #return np.sum(np.log(y_hat[y==1])) + np.sum(np.log(1 - y_hat[y==0]))

def log_lik_cont(y, y_hat, t):
    y_inds = np.searchsorted(t, y)
    return np.sum(np.log(y_hat[y_inds]))

def _margconds_from_intersection(prob_ab, prob_a):
    """
    :param prob_ab: joint probability of A and B
    :param prob_a: marginal probability of A
    :return: 
        a_given_b: conditional probability of A given B
        prob_b: marginal probability of B
    """
    prob_a = np.expand_dims(prob_a, axis=0)

    # conditional probability
    b_given_a = (prob_ab / np.sum(prob_ab, axis=(0, 1))) / prob_a
    prob_b = prob_a.t  @ b_given_a

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