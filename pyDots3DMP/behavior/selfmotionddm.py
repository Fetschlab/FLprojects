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
from scipy.optimize import minimize

from .utils import log_lik_bin, log_lik_cont, margconds_from_intersection
from .Accumulator import Accumulator

logger = logging.getLogger(__name__)

# %% ----------------------------------------------------------------
    
class SelfMotionDDM:
    PARAM_NAMES = ('kmult', 'bound', 'non_dec_time', 'wager_thr', 'wager_alpha')

    def __init__(
        self,
        grid_vec: np.ndarray,
        tvec: np.ndarray,
        kmult: list = [0.3, 0.3],
        bound: list = [1., 1., 1.],
        non_dec_time: list = [0.3],
        wager_thr: list = [1.],
        wager_alpha: list = [0.05],
        return_wager: bool = True,
        wager_axis: Optional[int] = None,
        stim_scaling: Union[tuple[np.ndarray, np.ndarray], bool] = True,
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
        # TODO add init option to set whether to use vectorized accumulator cdf/pdf calculations

        # initialize internal containers used by fit/predict
        self.init_params = {k: getattr(self, k) for k in self.PARAM_NAMES}
        self.params_ = self.init_params.copy() # will hold all params after fitting
        self.fixed_params = {}

        self.accumulators_ = {}  # cache of accumulators per condition if desired

       
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

            min_method = 'bads'
            # min_method = 'Nelder-Mead'
            
            if min_method == 'bads':

                lb = params_array * 0.25
                ub = params_array * 3.0
                plb = params_array * 0.5
                pub = params_array * 2.0
                bads_bounds = (lb, ub, plb, pub)

                options = {
                        "random_seed": 42,
                        # "uncertainty_handling": True,
                        "max_fun_evals": 300,
                        # "noise_final_samples": 100,
                        "display": "full"
                    }
                
                bads = BADS(
                    optim_fcn_part, 
                    params_array, 
                    *bads_bounds, 
                    options=options
                    )
                result = bads.optimize()

            else:
                result = minimize(
                    self._objective_fcn,
                    params_array,
                    args=(X, y),
                    method=min_method,
                    options={'maxiter': 5000, 'disp': True}
                    )
                
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

        params_array = np.asarray(params_array)
            
        # combine params array passed to objective function with fixed params
        # to reconstruct full params dict expected by custom predict method
        self._build_params_dict(params_array, self.param_end_inds, fixed_params)
        logger.info(params_array)
        # logger.info('Current params: %s', {k: [round(vv, 2) for vv in v] for k, v in self.params_.items()})

        t0_pred = time.perf_counter()
        y_pred, _ = self.predict(X, y)
        t1_pred = time.perf_counter() - t0_pred
        logger.info(f'single objective function prediction run took {t1_pred:.2f} seconds')
        
        # calculate log likelihoods for each output
        log_lik_choice = log_lik_bin(y['choice'].to_numpy(), y_pred['choice'].to_numpy()) / len(y)
        log_lik_pdw    = log_lik_bin(y['PDW'].to_numpy(), y_pred['PDW'].to_numpy()) / len(y)
        log_lik_rt     = log_lik_cont(y_pred['RT'].to_numpy()) / len(y)
        
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


    def predict(
        self,
        X,
        y=None,
        n_samples: int=1,
        cache_accumulators=False,
        seed=None
        ):
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
        
        mods = np.unique(X['modality']).astype(float)
        cohs = np.unique(X['coherence']).astype(float)
        deltas = np.unique(X['delta']).astype(float)
        hdgs, hdg_inds = np.unique(X['heading'], return_inverse=True)
        hdgs = hdgs.astype(float)

        k_scale = 1e3
        # get stimulus-driven urgency signals if specified, for scaling drifts
        if not self.stim_scaling:
            b_ves, b_vis = np.ones_like(self.tvec), np.ones_like(self.tvec)
            b_ves /= len(b_ves)
            b_vis /= len(b_vis)
        elif isinstance(self.stim_scaling, tuple):
            b_ves, b_vis = self.stim_scaling
        else:
            b_ves, b_vis = self.get_stim_urgs(self.tvec)
            k_scale = 10

        b_vals = [b_ves, b_vis, np.vstack((b_ves, b_vis)).T]

        # handle parameters per modality
        kves, kvis = self._handle_kmult(self.params_['kmult'], cohs.T, k_scale=k_scale) 
        bound = self._handle_param_mod(self.params_['bound'], mods)  
        non_dec_time = self._handle_param_mod(self.params_['non_dec_time'], mods)  
        thetas = self._handle_param_mod(self.params_['wager_thr'], mods)  
        alphas = self._handle_param_mod(self.params_['wager_alpha'], mods)  


        # initialize predictions dataframes
        predictions = pd.DataFrame(
            np.full((X.shape[0], 3), fill_value=np.nan), columns=['choice', 'PDW', 'RT']
            )
        pred_sample = deepcopy(predictions) if n_samples else None

        self.wager_maps = []
        if self.return_wager:

            # ves, vis, comb overall sensitivities
            k_vals_fixed = [kves, kvis.mean().item(), [kves, kvis.mean().item()]]
            
            for m, mod in enumerate(mods):

                # set accumulators with absolute drifts, for log odds mappings
                accumulator = Accumulator(grid_vec=self.grid_vec, tvec=self.tvec, bound=bound[m])

                abs_drifts, accumulator.tvec = self.calc_3dmp_drift_rates(
                    b_vals[m], k_vals_fixed[m], self.tvec[-1], hdgs[hdgs>=0], delta=0,
                    )

                # run the method of images - diffusion to bound to extract pdfs, cdfs, and LPO
                accumulator.apply_drifts(abs_drifts, hdgs[hdgs>=0])
                accumulator.compute_distrs(return_pdf=True) # get the pdfs for wager calculation

                if cache_accumulators:
                    self.accumulators_[('wager', mod)] = accumulator

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

                    if cache_accumulators:
                        self.accumulators_[(mod, coh, delta)] = accumulator
                    # get predictions for all trials in this condition
                    
                    # ====== CHOICE ======
                    p_right = np.clip(accumulator.p_corr_.T, 1e-10, 1-1e-10)
                    predictions.loc[trial_index, 'choice'] = p_right[hdg_inds[trial_index]]

                    # ====== WAGER ======
                
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
                            p_choice_given_wager, p_wager = margconds_from_intersection(
                                p_choice_and_wager, p_choice
                                )
                            p_wager += np.array([-alphas[m], alphas[m]]) * p_wager[0]
                            p_wager = np.clip(p_wager, 1e-100, 1-1e-100)

                            predictions.loc[trial_index, 'PDW'] = p_wager[0]

                            if n_samples:
                                pred_sample.loc[trial_index, 'PDW'] = rng.binomial(n_samples, p_wager[0], trial_index.sum()) / n_samples


                        # ====== RT ======

                        # first convolve model RT distribution with non-decision time
                        ndt_dist = norm.pdf(self.tvec, loc=non_dec_time[m], scale=0.2) #scale=self.params_['sigma_ndt'])
                        rt_dist = np.squeeze(accumulator.rt_dist_[h, :])
                        rt_dist = convolve(rt_dist, ndt_dist / ndt_dist.sum())
                        rt_dist = np.clip(rt_dist, 1e-10, a_max=None)

                        # trim to original length of tvec and renormalize to get posterior
                        rt_dist = rt_dist[:len(accumulator.tvec)]
                        rt_dist /= rt_dist.sum()

                        if y is not None:
                            # need original data here to get RT likelihoods
                            actual_rts = y.loc[trial_index, 'RT'].values
                            # dist_inds = [np.argmin(np.abs(self.tvec - rt)) for rt in actual_rts]
                            dist_inds = np.searchsorted(self.tvec, actual_rts)
                            dist_inds[dist_inds >= len(self.tvec)] = len(self.tvec) - 1  # cap at max index
                            predictions.loc[trial_index, 'RT'] = rt_dist[dist_inds]

                        if n_samples:
                            sampled_RTs = np.random.choice(
                                self.tvec, (trial_index.sum(), n_samples), replace=True, p=rt_dist
                                )
                            pred_sample.loc[trial_index, 'RT'] = sampled_RTs.mean(axis=1)

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


    def simulate(
        self,
        X,
        n_samples: int = 1,
        sample_dvs: bool = False,
        seed=None):
        """
        simulate data from model given X
        :param X: DataFrame with columns modality, coherence, delta, heading
        :param n_samples: number of samples to draw per trial (default 1)
        :param return_dvs: whether to sample DV trajectories
        :param seed: random seed for sampling
        :return: simulated DataFrame with columns choice, PDW, RT

        
        """
        
        _, preds = self.predict(
            X,
            y=None,
            n_samples=n_samples, 
            cache_accumulators=sample_dvs, 
            seed=seed
            )

        if sample_dvs:
            # build n_samples DV trajectories for unique conditions (from cached accumulator objects)

            # NOTE: IMCOMPLETE
            
            Xu = X.drop_duplicates().values
            preds = pd.DataFrame(np.repeat(Xu, n_samples, axis=0), columns=X.columns)
            
            for i, row in X.iterrows():
                mod = row['modality']
                coh = row['coherence']
                delta = row['delta']

                accum = self.accumulators_[(mod, coh, delta)]

                for hdg, drift in enumerate(zip(accum.drift_labels, accum.drift_rates)):

                    trial_inds = preds.index[
                        (preds["modality"]==mod) & (preds["coherence"]==coh) & \
                            (preds["delta"]==delta) & (preds["heading"]==hdg)
                    ]

                    for itr in range(n_samples):
                        dv = accum.dv(drift, sigma=np.array([1., 1.]))

                        is_hit_bnd = (dv >= accum.bound).any(axis=0)
                        t_bnd_cross = np.argmax((dv >= accum.bound) == 1, axis=0)

                        if is_hit_bnd.all():
                            # both accumulators hit the bound - choice and RT determined by whichever hit first
                            choice = np.argmin(t_bnd_cross)
                            rt_ind = t_bnd_cross[choice]
                            final_v = dv[rt_ind, choice ^ 1]  # losing accumulator at RT

                        elif ~is_hit_bnd.any():
                            # neither accumulator hits the bound
                            rt_ind = -1
                            choice = np.argmax(np.argmax(dv, axis=0))  # winner is whoever has max dv value
                            final_v = accum.bound[choice] - (dv[rt_ind, choice] - dv[rt_ind, choice ^ 1])
                            # wager odds map accounts for the distance between winner and loser, so this shifts up both accumulators
                            # as if the 'winner' did hit the bound, so we can do a consistent look-up on the wager odds map

                        else:
                            # only one hits the bound
                            choice = np.argmax(is_hit_bnd)
                            rt_ind = t_bnd_cross[choice]
                            final_v = dv[rt_ind, choice ^ 1]

                        if self.return_wager:

                            # look-up log odds threshold
                            grid_ind = np.argmin(np.abs(accum.grid_vec - final_v))

                            wager_accum = self.accumulators_[('wager', mod)]
                            
                            # log_odds = wager_odds_maps[m][rt_ind, grid_ind]
                            wager = int(wager_odds_above_threshold[m][rt_ind, grid_ind])
                            wager *= (np.random.random() > params['alpha'])  # incorporate base-rate of low bets
                            preds.loc[trial_inds[itr], 'PDW'] = wager

                        # flip choice result so that left choices = 0, right choices = 1 in the output
                        preds.loc[these_trials[tr], 'choice'] = choice ^ 1

                        # RT = decision time + non-decision time - motion onset latency
                        preds.loc[these_trials[tr], 'RT'] = \
                            orig_tvec[rt_ind] + non_dec_time[tr] - 0.3

        return preds
        
    
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