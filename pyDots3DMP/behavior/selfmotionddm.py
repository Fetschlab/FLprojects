# %% ----------------------------------------------------------------

from collections import namedtuple
from copy import deepcopy
from datetime import datetime
import logging
import itertools
import time
from typing import Union, Optional, Any

import numpy as np
import pandas as pd
from pybads import BADS
from scipy.signal import convolve
from scipy.stats import norm, skewnorm

from preprocessing import data_cleanup, format_onetargconf
from accumulator import Accumulator

# TODO list

# - get optimization running
# - get visualizations of model predictions
# - implement simulate method
# - abstract common parts of predict and simulate methods together
# - if no wager, drop wager params and handle
# - cleanup RT part
# - add logging 


# low-level to do 
# handle fixed params/params when calling predict before any fitting?
# maybe the way is to call fit with all fixed params...

# ISSUE - are params being set correctly on each fitting iteration - should be concatenating 
# fit and fixed for predict function # use print print print

# %% ----------------------------------------------------------------

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Create handlers
c_handler = logging.StreamHandler()
f_handler = logging.FileHandler(f'ddm_fitting_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log')
c_handler.setLevel(logging.DEBUG) 
f_handler.setLevel(logging.DEBUG) 

# Create formatters and add them to handlers, and add handlers to logger
c_format = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
f_format = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
c_handler.setFormatter(c_format)
f_handler.setFormatter(f_format)
logger.addHandler(c_handler)
logger.addHandler(f_handler)


OptimBounds = namedtuple('OptimBounds', ['lb', 'ub', 'plb', 'pub'])

def main():
    
    grid_vec = np.arange(-3, 0, 0.05)
    time_vec = np.arange(0, 2, 0.05)
       
    init_params = {
        'kmult': [0.3, 0.3], 
        'bound': [1., 1., 1.],
        'non_dec_time': [0.3],
        'wager_thr': [1, 1, 1],
        'wager_alpha': [0.05]
    }
    accum = SelfMotionDDM(grid_vec=grid_vec, tvec=time_vec, **init_params, 
                                stim_scaling=False, return_wager=False)

    datafilepath = "/Users/stevenjerjian/Desktop/Academia/FetschLab/PLDAPS_data/dataStructs/lucio_20220512-20230606.csv"
    data = data_cleanup(datafilepath)  # this is a hack function to quickly load and clean the data, could be improved/generalized with options
    data = format_onetargconf(data, remove_one_targ=True)
    data = data.reset_index(drop=True)

    X = data[['modality', 'coherence', 'delta', 'heading']]
    y = data[['choice', 'PDW', 'RT']].astype({'choice': 'int', 'PDW': 'int'})
    y['RT'] += 0.3 # to reverse the downshift from motion platform latency

    accum.fit(X, y)
    y_pred, y_pred_samp = accum.predict(X, y, n_samples=10, seed=1)
    # y_pred_samp['RT'] -= 0.3

    # print(y.head())
    # print(y_pred.head())
    # print(y_pred_samp.head())

# %% ----------------------------------------------------------------
    
class SelfMotionDDM:
    PARAM_NAMES = ('kmult', 'bound', 'non_dec_time', 'wager_thr', 'wager_alpha')

    def __init__(
        self,
        grid_vec: np.ndarray,
        tvec: np.ndarray,
        kmult: list = [0.3, 0.3],
        bound: list | int | float = 1.,
        non_dec_time: list | int | float = 0.3,
        wager_thr: list | int | float = 1.,
        wager_alpha: list | int | float = 0.05,
        return_wager: bool = True,
        wager_axis: Optional[int] = None,
        stim_scaling: bool | tuple[np.ndarray, np.ndarray] = True,
        ):
        """3DMP accumulator model with confidence/wager readout
        :param grid_vec: vector of DV grid points
        :param tvec: vector of time points
        :param kmult: list of k multipliers [k_ves, k_vis]
        :param bound: list of bounds per modality [ves, vis, comb]      
        :param non_dec_time: list of non-decision times per modality
        :param wager_thr: list of wager thresholds per modality
        :param wager_alpha: list of wager alpha parameters per modality
        :param return_wager: whether to compute wager predictions
        :param wager_axis: axis for wager calculation (None = log odds)
        :param stim_scaling: whether to use stimulus-driven urgency signals
        """
        self.grid_vec = grid_vec
        self.tvec = tvec
        self.kmult = kmult
        self.bound = bound
        self.non_dec_time = non_dec_time
        self.wager_thr = wager_thr
        self.wager_alpha = wager_alpha
        self.return_wager = return_wager
        self.wager_axis = wager_axis
        self.stim_scaling = stim_scaling

        # initialize internal containers used by fit/predict
        self.init_params = {k: getattr(self, k) for k in self.PARAM_NAMES}
        self.params_ = self.init_params.copy() # will hold all params after fitting
        self.fixed_params = {}

       
    def fit(
        self,
        X: pd.DataFrame,
        y: pd.DataFrame,
        fixed_params: Optional[list[str]]=None,
        # optim_bounds: Optional[Union[OptimBounds, dict]]=None,
        ) -> 'SelfMotionDDM':
        """fit model to data in X and y, with optional fixed parameters"""

        logger.info('Starting model fitting')
        
        self.n_features_in_ = len(X.columns)

        # get list of fittable parameters, if any
        self.params_ = self.init_params.copy() 
        params_list = self._get_fit_params(fixed_params)

        if params_list:
            # concatenate into single array for passing to optimization function
            # but store original list lengths for reconstructing dict later
            params_array = np.array(list(itertools.chain(*params_list)))
            self.param_end_inds = list(itertools.accumulate(map(len, params_list)))

            # pass data as fixed inputs to objective function
            optim_fcn_part = lambda params: self._objective_fcn(params, X, y)

            # TODO allow user to specify these instead or with init params
            lb = params_array * 0.3
            ub = params_array * 2
            plb = params_array * 0.5
            pub = params_array * 1.5
            bads_bounds = (lb, ub, plb, pub)
            
            bads = BADS(
                optim_fcn_part, 
                params_array, 
                *bads_bounds, 
                options = {'random_seed': 1234}
                )
            result = bads.optimize()

            # at the end, store fitted params back into dict, with fixed ones
            self._build_params_dict(result.x, self.param_end_inds)

        return self

    def _objective_fcn(
        self,
        params_array: np.ndarray,
        X: pd.DataFrame,
        y: pd.DataFrame,
        fixed_params: Optional[dict[str, Any]]=None,
        ) -> float:
        """
        objective function for optimization - negative log likelihood of data given model params
        return float for minimization
        """
        
        # use params array and fixed params to reconstruct full params
        self._build_params_dict(params_array, self.param_end_inds, fixed_params)
        logger.info('Current params: %s', {k: [round(vv, 2) for vv in v] for k, v in self.params_.items()})

        t0_pred = time.perf_counter()
        y_pred, _ = self.predict(X, y)
        t1_pred = time.perf_counter() - t0_pred
        logger.info(f'single objective function prediction run took {t1_pred:.2f} seconds')
        
        # calculate log likelihoods for each output
        log_lik_choice = log_lik_bin(y['choice'].to_numpy(), y_pred['choice'].to_numpy()) / len(y)
        log_lik_pdw    = log_lik_bin(y['PDW'].to_numpy(), y_pred['PDW'].to_numpy()) / len(y)
        log_lik_rt     = log_lik_cont(y['RT'].to_numpy(), y_pred['RT'].to_numpy(), self.tvec) / len(y)

        self.log_lik_ = {
            'choice': log_lik_choice,
            'pdw': log_lik_pdw,
            'rt': log_lik_rt
        }
        if self.return_wager:
            logger.info('Log likelihoods - choice: %.2f, PDW: %.2f, RT: %.2f', 
                        log_lik_choice, log_lik_pdw, log_lik_rt)
            self.neg_llh_ = -sum([log_lik_choice, log_lik_pdw, log_lik_rt])
        else:
            logger.info('Log likelihoods - choice: %.2f, RT: %.2f', 
                        log_lik_choice, log_lik_rt)
            self.neg_llh_ = -sum([log_lik_choice, log_lik_rt])
        logger.info('Total loss:\t%.2f', self.neg_llh_)

        return self.neg_llh_
    

    def predict(self, X, y, n_samples: int=0, seed=None):
        """
        generate model predictions for data in X
        :param X: DataFrame with columns modality, coherence, delta, heading
        :param y: DataFrame with columns choice, PDW, RT (for RT likelihood calculation)
        :param n_samples: number of samples to draw for probabilistic predictions (0 = none)
        :param seed: random seed for sampling
        :return: 
            predictions - DataFrame with columns choice, PDW, RT (predicted likelihoods)
            pred_sample - DataFrame with sampled predictions (n_samples > 0 -> how many draws from binomial for choice/PDW)
        """
        
        rng = np.random.RandomState(seed)
        
        mods = np.unique(X['modality'])
        cohs = np.unique(X['coherence'])
        deltas = np.unique(X['delta'])
        hdgs, hdg_inds = np.unique(X['heading'], return_inverse=True)

        cond_groups = X[['modality','coherence','delta']].drop_duplicates(
            ).sort_values(by=['modality','coherence','delta'])

        # TODO clean these up better
        kves, kvis = self._handle_kmult(self.params_['kmult'], cohs.T, k_scale=1) 
        bound = self._handle_param_mod(self.params_['bound'], mods)  
        non_dec_time = self._handle_param_mod(self.params_['non_dec_time'], mods)  
        thetas = self._handle_param_mod(self.params_['wager_thr'], mods)  
        alphas = self._handle_param_mod(self.params_['wager_alpha'], mods)  

        if not self.stim_scaling:
            kves, kvis = self._handle_kmult(self.params_['kmult'], cohs.T, k_scale=100) 

            b_ves, b_vis = np.ones_like(self.tvec), np.ones_like(self.tvec)

            b_ves /= len(b_ves)
            b_vis /= len(b_vis)
        elif isinstance(self.stim_scaling, tuple):
            b_ves, b_vis = self.stim_scaling
        else:
            b_ves, b_vis = self.get_stim_urgs(self.tvec)

        b_vals = [b_ves, b_vis, np.vstack((b_ves, b_vis)).T]

        predictions = pd.DataFrame(np.full((X.shape[0], 3), fill_value=np.nan), 
                                   columns=['choice', 'PDW', 'RT'])
        pred_sample = deepcopy(predictions) if n_samples else None

        if self.return_wager:
            logger.info('Calculating wager maps')
            self.wager_maps = []
            k_vals_fixed = [kves, np.mean(kvis), [kves, np.mean(kvis)]]
            for m, mod in enumerate(mods):

                # set accumulators with absolute drifts, for log odds mappings
                
                accumulator = Accumulator(grid_vec=self.grid_vec, tvec=self.tvec, bound=bound[m])

                abs_drifts, accumulator.tvec = self.calc_3dmp_drift_rates(
                    b_vals[m], k_vals_fixed[m], self.tvec[-1], hdgs[hdgs>=0], delta=0,
                    )

                # run the method of images - diffusion to bound to extract pdfs, cdfs, and LPO
                # positive headings only (assume symmetry around 0 for log odds mapping)
                accumulator.apply_drifts(abs_drifts, hdgs[hdgs>=0])
                accumulator.compute_distrs(return_pdf=True) # get the pdfs necessarily

                if self.wager_axis is None:
                    log_odds_map = accumulator.log_posterior_odds()
                    self.wager_maps.append(log_odds_map)
                else:
                    raise NotImplementedError('alternatives to log odds not yet implemented')

                wager_is_high = [p >= theta for p, theta in zip(self.wager_maps, thetas)]

        # now loop over coherences and deltas with one accumulator each for actual predictions
        for c, coh in enumerate(cohs):
            k_vals = [kves, kvis[c], [kves, kvis[c]]]

            for m, mod in enumerate(mods):
                for d, delta in enumerate(deltas):

                    trial_index = (X['modality'] == mod) & (X['coherence'] == coh) & (X['delta'] == delta)
                    trial_index = trial_index.to_numpy()

                    # non-valid conditions
                    if (delta != 0 and mod < 3) or (c > 0 and mod == 1) or trial_index.sum()==0:
                        continue       

                    # set up accumulator for this condition
                    accumulator = Accumulator(grid_vec=self.grid_vec, tvec=self.tvec, bound=bound[m])
                    drifts, accumulator.tvec = self.calc_3dmp_drift_rates(
                        b_vals[m], k_vals[m], self.tvec[-1], hdgs, delta=delta,
                        )
                    # this time use signed headings
                    accumulator.apply_drifts(drifts, hdgs) 

                    # run the method of images - diffusion to bound to extract pdfs, cdfs, and LPO
                    accumulator.compute_distrs(return_pdf=self.return_wager)

                    # get predictions for all trials in this condition
                    
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
                        if self.return_wager:
                            # select pdf for losing race, given correct or incorrect
                            pxt_up = np.squeeze(accumulator.up_lose_pdf_[h, :, :])
                            pxt_lo = np.squeeze(accumulator.lo_lose_pdf_[h, :, :])
                            total_p = np.sum(pxt_up + pxt_lo)
                            pxt_up /= total_p
                            pxt_lo /= total_p

                            p_choice_and_wager = np.array(
                                [
                                    [
                                        np.sum(pxt_up[wager_is_high[m]]),   # pRight+High
                                        np.sum(pxt_up[~wager_is_high[m]])   # pRight+Low
                                    ],   
                                    [
                                        np.sum(pxt_lo[wager_is_high[m]]),   # pLeft+High
                                        np.sum(pxt_lo[~wager_is_high[m]])   # pLeft+Low
                                    ]
                                ]
                            )

                            # calculate p_wager using Bayes rule, then factor in base rate of low bets ("alpha")
                            p_choice_given_wager, p_wager = _margconds_from_intersection(p_choice_and_wager, p_choice)
                            p_wager += np.array([-self.params_['wager_alpha'][m], self.params_['wager_alpha'][m]]) * p_wager[0]
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
                        actual_rts = y.loc[trial_index, 'RT'].values
                        dist_inds = [np.argmin(np.abs(self.tvec - rt)) for rt in actual_rts]
                        predictions.loc[trial_index, 'RT'] = rt_dist[dist_inds]

                        if n_samples:
                            pred_sample.loc[trial_index, 'RT'] = rng.choice(
                                self.tvec, trial_index.sum(), replace=True, p=rt_dist
                                )
        
        return predictions, pred_sample
    

    def _build_params_dict(
        self,
        params_array: np.ndarray,
        end_inds: list,
        fixed_params: Optional[dict[str, Any]]=None,
        ) -> None:
        """reconstruct full params dictionary from params array and fixed params dictionary"""

        fixed_params = fixed_params if fixed_params is not None else self.fixed_params

        start_ind = 0
        for p, param in enumerate(self.fit_param_names):
            if p > 0:
                start_ind = end_inds[p-1]
            self.params_[param] = params_array[start_ind:end_inds[p]].tolist()
        self.params_ = self.params_ | fixed_params

        # for pn, p in self.params_.items():
        #     print(f'{pn}:', sep='\t')
        #     for elem in p:
        #         print(f'{elem:.2f}', sep=',')


    def _get_fit_params(self, fixed_params) -> list:
        """return list of fittable parameter values (from initial parameters dict)"""

        if fixed_params:
            self.fixed_params = {k: self.init_params[k] for k in fixed_params}
            
        self.fit_param_names = [k for k in self.init_params.keys() if k not in self.fixed_params.keys()]
        logger.info(f'Fixed parameters: {self.fixed_params}')
        logger.info(f'Fitting parameters: {self.fit_param_names}')

        return [self.init_params[k] for k in self.fit_param_names]


    def simulate(self, X, n_samples: Optional[int] = 100, seed=None):
        """
        simulate data from model given X
        :param X: DataFrame with columns modality, coherence, delta, heading
        :return: simulated DataFrame with columns choice, PDW, RT
        """
        _, preds = self.predict(X, None, n_samples=1)
        
        raise NotImplementedError('simulate method not yet implemented')
        
    
    @staticmethod
    def _handle_kmult(kmult, cohs, k_scale=1):
        """break down kmult into (kves, kvis)"""
        kmult = [k*k_scale for k in kmult]
        kves = kmult[0]
        if len(kmult) == 1:
            kvis = kmult * cohs
            kves = np.mean(kvis)  # ves lies between vis cohs
        if len(kmult) == 2:
            kvis = kmult[1] * cohs # vis scaled by coherence
        else:
            kvis = kmult[1:]  # 3 independent kmults
        return kves, kvis

    @staticmethod
    def _handle_param_mod(param, mods):
        """broadcast param to match number of mods"""
        n_mods = len(mods)
        
        if np.isscalar(param):
            return [param] * n_mods

        param_list = list(param)
        if len(param_list) == n_mods:
            return param_list
        if len(param_list) == 1:
            return param_list * n_mods

        raise ValueError("Length of param list does not match number of modalities")


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
    def calc_3dmp_drift_rates(
        b_t: np.ndarray,
        b_k: Union[float, tuple[float, float]],
        tmax: float,
        hdgs: np.ndarray,
        delta: Optional[float] = 0.0, 
        cue_weights: Optional[tuple[float, float]] = None
        ) -> tuple[np.ndarray, np.ndarray]:
        """
        Calculate drift rates
        if cue_weights is None, optimal cue weights calculated from k's
        """

        sin_hdgs = np.sin(np.deg2rad(hdgs))

        cumul_bt = np.cumsum((b_t**2)/(b_t**2).sum(axis=0), axis=0)

        if isinstance(b_k, (int, float)):
            # only one sensitivity and time-course - ves or vis (logic is the same)
            b_t = b_t.reshape(-1, 1)
            tvec = cumul_bt * tmax
            drifts = np.cumsum(b_t**2 * b_k * sin_hdgs, axis=0)
            # drifts2 = b_k * sin_uhdgs   # w/o stim scaling, reduces to this

        elif len(b_k) == 2:
            # two sensitivities/time-courses - combined condition

            b_k = np.array(b_k, dtype=float)
            if cue_weights is None:
                k2 = b_k**2
                cue_weights = np.sqrt(k2 / k2.sum())
                
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
## GENERAL HELPR FUNCTIONS

def log_lik_bin(y, y_hat):
    # return np.sum(y * np.log(y_hat) + (1 - y) * np.log(y_hat))
    return np.sum(np.log(y_hat[y==1])) + np.sum(np.log(1 - y_hat[y==0]))

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

    prob_a = prob_a.reshape(-1, 1)  # make it 2-D, for element-wise and matrix mults below
    b_given_a = (prob_ab / np.sum(prob_ab)) / prob_a
    prob_b = prob_a.T @ b_given_a
    prob_b[prob_b==0] = np.finfo(np.float64).tiny
    a_given_b = b_given_a * prob_a / prob_b

    # assert np.allclose(prob_b, 1.0) and np.allclose(np.sum(a_given_b, axis=0), [1.0, 1.0])
    # assert np.allclose(prob_a, 1.0) and np.allclose(prob_ab, 1.0)
    
    return a_given_b, prob_b.flatten()
    
    # TODO attempt to do for multiple headings at once...work in progress
    # proba_full = np.expand_dims(prob_a, axis=1)
    # b_given_a = (prob_ab / np.sum(prob_ab, axis=(0, 1))) / proba_full
    # prob_b = np.zeros_like(prob_a)
    # for d in range(prob_b.shape[1]):
    #     prob_b[:, d] = prob_a[:, d].T @ b_given_a[:, :, d]
    
    # a_given_b = b_given_a * proba_full / np.expand_dims(prob_b, axis=0)
    # return a_given_b, prob_b


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
    main()