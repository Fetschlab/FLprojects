#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 14:19:37 2023

@author: stevenjerjian
"""

# generate psth over trials(fixed alignements, or time-warping)
# average psth over conditions
# fano factor calculations
# multiple alignments, use map/apply? test in main

# %% imports

import numpy as np
import pandas as pd
import xarray as xr

import matplotlib.pyplot as plt
from matplotlib import cm
import seaborn as sns

import pdb

# move corr functions to separate module

from scipy.ndimage import convolve1d #, gaussian_filter1d
from scipy.signal import gaussian


# import matplotlib.pyplot as plt
# import seaborn as sns
# custom_params = {"axes.spines.right": False, "axes.spines.top": False}
# sns.set_theme(context="notebook", style="ticks", rc=custom_params)

# %% per-trial spike counts/firing rates

def trial_psth(spiketimes, align_ev, trange,
               binsize=0.05, sm_params={},
               all_trials=False, normalize=True):

    nTr = align_ev.shape[0]

    if nTr == 0:
        return
    else:
        if align_ev.ndim == 2:
            align_ev = np.sort(align_ev, axis=1)
            ev_order = np.argsort(align_ev, axis=1)
            # assert check here that ev_order is consistent on every trial?
            which_ev = ev_order[0, 0]  # align to this event
        else:  # only one event provided, duplicate it for below
            align_ev = np.expand_dims(align_ev, axis=1)
            align_ev = np.tile(align_ev, (1, 2))
            which_ev = 0

        nantrs = np.any(np.isnan(align_ev), axis=1)

        if np.sum(nantrs) > 0:
            print(f'Dropping {np.sum(nantrs)} trials with missing event (NaN)')

        align_ev = align_ev[~nantrs, :]

        nTr = align_ev.shape[0]  # recalculate after bad trs removed

        tr_starts = align_ev[:, 0] + trange[0]
        tr_ends = align_ev[:, 1] + trange[1]
        durs = tr_ends - tr_starts

        # compute 'corrected' tStart and tEnd based on align_ev input
        # 'unified'
        if which_ev == 1:
            tstarts_new = tr_starts - align_ev[:, 1]
            tstart_new = np.min(tstarts_new)
            tend_new = trange[1]
            tends_new = tend_new.repeat(tstarts_new.shape[0])
            # TODO maybe just use np.repeat here, rather than itertools?

        elif which_ev == 0:
            tends_new = tr_ends - align_ev[:, 0]
            tend_new = np.max(tends_new)
            tstart_new = trange[0]
            tstarts_new = tstart_new.repeat(tends_new.shape[0])

        if binsize > 0:

            if trange[0] < 0 and trange[1] > 0:
                # ensure that time '0' is in between two bins exactly
                x0 = np.arange(0, tstart_new-binsize, -binsize)
                x1 = np.arange(0, tend_new+binsize, binsize)
                x = np.hstack((x0[::-1, ], x1[1:, ]))
            else:
                x = np.arange(tstart_new, tend_new+binsize, binsize)

            fr_out = np.full([nTr, x.shape[0]-1], np.nan)
        else:
            fr_out = np.full(nTr, np.nan)
            x = durs

        spktimes_aligned = []
        if spiketimes.any():
            itr_start, itr_end = 1, nTr

            if not all_trials:
                itr_start = np.argmin(np.abs(tr_starts - spiketimes[0]))
                itr_end = np.argmin(np.abs(tr_ends - spiketimes[-1]))

            for itr in range(itr_start, itr_end+1):
                spk_inds = np.logical_and(spiketimes >= tr_starts[itr],
                                          spiketimes <= tr_ends[itr])

                if binsize == 0:
                    fr_out[itr] = np.sum(spk_inds)

                else:
                    inds_t = spiketimes[spk_inds] - align_ev[itr, which_ev]
                    fr_out[itr, :], _ = np.histogram(inds_t, x)

                    spktimes_aligned.append(inds_t)

                    # set nans outside the range of align/trange for each trial
                    if which_ev == 0:
                        end_pos = np.argmin(abs(x - tends_new[itr]))
                        fr_out[itr, end_pos:] = np.nan
                    elif which_ev == 1:
                        start_pos = np.argmin(abs(x - tstarts_new[itr]))
                        fr_out[itr, :start_pos] = np.nan

            if binsize > 0:
                x = x[:-1] + np.diff(x)/2  # shift x values to bin centers

                if sm_params:
                    fr_out = smooth_counts(fr_out, params=sm_params)

                if normalize:
                    fr_out /= binsize
            elif binsize == 0:
                if normalize:
                    fr_out /= x

        return fr_out, x, spktimes_aligned


def smooth_counts(raw_fr, params={'kind': 'boxcar', 'binsize': 0.05,
                                  'width': 0.4, 'sigma': 0.25}):

    N = int(np.ceil(params['width'] / params['binsize']))  # width, in bins

    if params['kind'] == 'boxcar':

        win = np.ones(N) / N

    elif params['kind'] == 'gaussian':

        alpha = (N - 1) / (2 * (params['sigma'] / params['binsize']))
        win = gaussian(N, std=alpha)
        # win /= np.sum(win)  # win is already normalized to win in scipy

        # smoothed_fr = gaussian_filter1d(raw_fr, sigma=alpha)

    elif params['kind'] == 'CHG':  # causal half-gaussian

        alpha = (N - 1) / (2 * (params['sigma'] / params['binsize']))
        win = gaussian(N, std=alpha)
        win[:(N//2)-1] = 0
        win /= np.sum(win)  # renormalize here

    smoothed_fr = convolve1d(raw_fr, win, axis=0, mode='nearest')
    return smoothed_fr


# %% Extract firing rates for population

# now defined as a method in Population dataclass

def get_aligned_rates(popn, align=['stimOn'], trange=np.array([[-2, 3]]),
                      binsize=0.05, sm_params={'kind': 'boxcar',
                                               'binsize': 0.05,
                                               'width': 0.4, 'sigma': 0.25},
                      condlabels=['modality', 'coherence', 'heading'],
                      return_Dataset=False):

    good_trs = popn.events['goodtrial'].to_numpy(dtype='bool')
    condlist = popn.events[condlabels].loc[good_trs, :]

    rates = []
    tvecs = []
    align_lst = []

    for al, t_r in zip(align, trange):

        if popn.events.loc[good_trs, al].isna().all(axis=0).any():
            raise ValueError(al)

        align_ev = popn.events.loc[good_trs, al].to_numpy(dtype='float64')

        # get spike counts and relevant t_vec for each unit - thanks chatGPT!
        # trial_psth in list comprehesnsion is going to generate a list of
        # tuples, the zip(*iter) syntax allows us to unpack the tuples into
        # separate variables
        spike_counts, t_vec, _ = \
            zip(*[(trial_psth(unit.spiketimes, align_ev,
                              t_r, binsize, sm_params=sm_params))
                  for unit in popn.units])

        rates.append(np.asarray(spike_counts))
        tvecs.append(np.asarray(t_vec[0]))
        # previously tried dict with key as alignment event

        align_lst.append([al]*len(t_vec[0]))

    align_arr = np.concatenate(align_lst)
    unitlabels = np.array([u.clus_group for u in popn.units])
    # unit_ids  = np.array([u.clus_id for u in popn.units])

    if return_Dataset:
        # now construct a dataset
        arr = xr.DataArray(np.concatenate(rates, axis=2),
                           coords=[unitlabels,
                                   condlist.index,
                                   np.concatenate(tvecs)],
                           dims=['unit', 'trial', 'time'])

        cond_coords = {condlabels[c]: ('trial', condlist.iloc[:, c])
                       for c in range(len(condlist.columns))}

        ds = xr.Dataset({'firing_rate': arr},
                        coords={'align_event': ('time', align_arr),
                                **cond_coords})

        ds.attrs['rec_date'] = popn.rec_date
        ds.attrs['area'] = popn.area

        return ds

    else:
        # return separate vars
        return rates, tvecs, condlist, align_lst


def concat_aligned_rates(f_rates, tvecs=None):
    """

    Parameters
    ----------
    f_rates : list
        list of firing rates, where each element is a different interval
        containing units x trials or units x trials x binned time
    tvecs : list or None
        corresponding time vector (only for time-resolved f_rates)
        default: None

    Returns
    -------
    rates_cat : TYPE
        concatenated rates
    len_intervals : TYPE
        array containing the length of each interval, for future split

    """

    if tvecs is not None:
        # concatenate all different alignments...
        # but store the lens for later splits
        len_intervals = [np.asarray(list(map(lambda x: x.size, t)),
                         dtype='int').cumsum() for t in tvecs]
        rates_cat = list(map(lambda x: np.concatenate(x, axis=2), f_rates))
    else:
        # each 'interval' is length 1 if binsize was set to 0
        len_intervals = [np.ones(len(r), dtype=int).cumsum() for r in f_rates]
        rates_cat = list(map(np.dstack, f_rates))

    return rates_cat, len_intervals


def lowfr_units(f_rates, minfr=0):

    mean_fr = np.squeeze(np.mean(f_rates, axis=1))
    lowfr_units = np.logical_or(np.isnan(mean_fr), mean_fr <= minfr)

    return lowfr_units

    

# %% trial condition helpers

def condition_index(condlist, cond_groups=None):

    assert isinstance(condlist, pd.DataFrame)

    if cond_groups is None:
        cond_groups, ic = np.unique(condlist.to_numpy('float64'), axis=0,
                                    return_inverse=True)
        cond_groups = pd.DataFrame(cond_groups, columns=condlist.columns)

    else:
        # cond_groups user-specified
        assert isinstance(cond_groups, pd.DataFrame)
        cond_groups = cond_groups.loc[:, condlist.columns]

        # fill with nan?
        ic = np.full(condlist.shape[0], fill_value=-1, dtype=int)
        for i, cond in enumerate(cond_groups.values):
            ic[(condlist == cond).all(axis=1)] = i

    nC = len(cond_groups.index)
    return ic, nC, cond_groups


def condition_averages_ds(ds, *args):
    # xarray groupby functions seem a bit clunky...
    pass


def condition_averages(f_rates, condlist, cond_groups=None):

    ic, nC, cond_groups = condition_index(condlist, cond_groups)

    cond_fr = np.full((f_rates.shape[0], nC, f_rates.shape[2]), np.nan)
    cond_sem = np.full((f_rates.shape[0], nC, f_rates.shape[2]), np.nan)

    for c in range(nC):
        if np.sum(ic == c):
            cond_fr[:, c, :] = np.mean(f_rates[:, ic == c, :], axis=1)
            cond_sem[:, c, :] = np.std(f_rates[:, ic == c, :],
                                       axis=1) / np.sqrt(np.sum(ic == c))

    return cond_fr, cond_sem, cond_groups


# %% plot functions


# plot raster
# dataframe containing conditions on each trial, and list of spikes on each trial
# list of spikes can be multiple columns for multiple alignments
# then should loop over these and put them on subplots within same subplot


# align_ev = events.loc[good_trs, align_ev].to_numpy(dtype='float64')

# fr, x, spks = trial_psth(spiketimes, align_ev, trange,
#                          binsize, sm_params)

# df = condlist[['modality','coherence','heading']]
# df['spikes'] = spks

# cmap = mpl.cm.get_cmap('PiYG')
# cvals = cmap(list(np.around((hdgs/90+1)/2,2)))


def spks_eventplot(df, row, col, hue):

    if row and col:
        split_conds = [row, col]
        g = sns.FacetGrid(df, row=row, col=col, sharex=True)
    elif row:
        split_conds = row
        g = sns.FacetGrid(df, row=row, sharex=True)
    else:
        split_conds = col
        g = sns.FacetGrid(df, col=col, sharey=True)

    for c, cond_df in df.groupby(split_conds):

        ax = g.axes_dict[c]
        ic, nC, cond_groups = condition_index(cond_df[split_conds])

        spks_c = cond_df['spikes'].to_numpy(dtype='object')
        ic2 = np.sort(ic)
        order = np.argsort(ic)

        # spks_c = list(spks_c[order])

        # colors = [palette[i] for i in ic2]
        # ax.eventplot(spks_c, lineoffsets=list(range(len(spks_c))),
        #              colors=colors)

        ax.eventplot(spks_c, lineoffsets=list(order))
        ax.invert_yaxis()
    plt.show()

    return g


def plot_psth():
    ...

# plot psth (time-res)

# plot condition averages (non-time res, heading on x-axis)

# plot tuning curves (with or without error bars/dots to show rsc)
# plot noise correlation (rsc) scatter plot


#def plot_psth(df, palette, type=1):

    # TODO
    # code for rasterized spike plots
    # rescale x-axes according to length of time vector
    # show coherences

    # if type == 'avg':
    #     sns.lineplot(data=df,
    #                  x='heading', y='spike counts',
    #                  estimator='mean', errorbar='se', err_style='bars',
    #                  hue=df[condlabels[:-1]].apply(tuple, axis=1))
    # elif type == 'time':
    #     g = sns.relplot(data=df,
    #                     x='time', y='firing_rate',
    #                     hue='heading', row='modality', col='align',
    #                     palette=palette,
    #                     facet_kws=dict(sharex=False),
    #                     kind='line', aspect=.75)
