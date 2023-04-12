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

# %%

import pdb

import numpy as np
import pandas as pd
import xarray as xr

import pickle
from pathlib import PurePath
from itertools import repeat
from scipy.ndimage import convolve1d, gaussian_filter1d


import seaborn as sns
custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(context="notebook", style="ticks", rc=custom_params)

# %% trial_psth function


def trial_psth(spiketimes, align_ev, trange=np.array([-1, 1]),
               binsize=0.05, smooth_params=('boxcar', 5),
               all_trials=False, normalize=True):

    nTr = align_ev.shape[0]

    if nTr == 0:
        return
    else:
        if align_ev.ndim == 2:
            align_ev = np.sort(align_ev, axis=1)
            ev_order = np.argsort(align_ev, axis=1)
            which_ev = ev_order[0]  # align to this event
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

        elif which_ev == 0:
            tends_new = tr_ends - align_ev[:, 0]
            tend_new = np.max(tends_new)
            tstart_new = trange[0]
            tstarts_new = tstart_new.repeat(tends_new.shape[0])

        if binsize > 0:

            if trange[0] < 0 and trange[1] > 0:
                # ensure that time '0' is in between two bins exactly
                x0 = np.arange(0, tstart_new-1e-3, -binsize)
                x1 = np.arange(0, tend_new, binsize)
                x = np.hstack((x0[::-1, ], x1[1:, ]))
            else:
                x = np.arange(tstart_new, tend_new+1e-3, binsize)

            fr_out = np.full([nTr, x.shape[0]-1], np.nan)
        else:
            fr_out = np.full(nTr, np.nan)

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

                    x = durs
                    if normalize:
                        fr_out /= x

                else:
                    inds_t = spiketimes[spk_inds] - align_ev[itr, which_ev]
                    fr_out[itr, :], _ = np.histogram(inds_t, x)

                    if which_ev == 0:
                        end_pos = np.argmin(abs(x - tends_new[itr]))
                        fr_out[itr, end_pos:] = np.nan
                    elif which_ev == 1:
                        start_pos = np.argmin(abs(x - tstarts_new[itr]))
                        fr_out[itr, :start_pos] = np.nan

                    # shift x values to bin centers
                    x = x[:-1] + np.diff(x)/2

                    fr_out = smooth_counts(fr_out, smooth_params=smooth_params)

                    if normalize:
                        fr_out /= binsize

        return fr_out, x, which_ev


# %%

def smooth_counts(raw_fr, smooth_params=('boxcar', 5)):

    if smooth_params[0] == 'boxcar':

        kw = smooth_params[1]  # width, in bins
        kernel = np.ones(kw) / kw
        smoothed_fr = convolve1d(raw_fr, kernel,
                                     axis=0, mode='reflect')

    elif smooth_params[0] == 'gaussian':

        sigma = smooth_params[1]
        smoothed_fr = gaussian_filter1d(raw_fr, sigma=sigma)

    # TODO causal filter

    return smoothed_fr

# %% Extract firing rates for population


def get_aligned_rates(popn, align=['stimOn'], trange=np.array([-2, 2]),
                      binsize=0.05, smooth_params=('boxcar', 5),
                      condlabels=['modality', 'coherence', 'heading'],
                      return_as_Dataset=False):

    # TODO somewhere need to get an update to tvec to be all relative to one event?
    # so that we can plot on the same time axis if we want

    good_trs = popn.events[align].notna().to_numpy(dtype='bool').all(axis=1)
    condlist = popn.events[condlabels].loc[good_trs, :]

    avg_rates = []
    durs = []
    rates = {}
    tvecs = {}
    align_lst = []

    for al, t_r in zip(align, trange):

        # good trials for this alignment event
        align_ev = popn.events.loc[good_trs, al].to_numpy(dtype='float64')

        # thanks chatGPT!
        # get spike counts and t_vec for each unit
        spike_counts, t_vec, total_counts, dur, which_ev = \
            zip(*[(trial_psth(unit.spiketimes, align_ev, t_r, binsize,
                              smooth_params))
                  for unit in popn.units])
            
        avg_rates.append(np.asarray(total_counts))
        durs.append(np.asarray(dur[0]))
        
        if return_as_Dataset:
            rates.append(np.asarray(spike_counts))
            tvecs.append(np.asarray(t_vec[0]))
        else:
            # do it as a dictionary, so alignment 'event' is saved as key
            rates[al] = np.asarray(spike_counts)
            tvecs[al] = np.asarray(t_vec[0])

        align_lst.append(np.asarray(list(repeat(al, len(t_vec[0])))))

    align_arr = np.concatenate(align_lst)
    u_labels = np.array([u.clus_group for u in popn.units])

    avg_rates = np.dstack(avg_rates)
    durs = np.dstack(durs)

    if return_as_Dataset:
        # now construct a dataset
        arr = xr.DataArray(np.concatenate(rates, axis=2),
                           coords=[u_labels,
                                   condlist.index,
                                   np.concatenate(tvecs)],
                           dims=['unit', 'trial', 'time'])

        cond_coords = {condlabels[c]: ('trial', condlist.iloc[:, c])
                       for c in range(len(condlist.columns))}

        ds = xr.Dataset({'firing_rate': arr},
                        coords={'align_event': ('time', align_arr),
                                **cond_coords})

        # TODO store some more ds.attrs
        ds.attrs['rec_date'] = popn.rec_date
        ds.attrs['area'] = popn.area

        return ds

    else:
        # return separate vars
        return rates, avg_rates, u_labels, condlist, tvecs, durs


# %% Average across conditions

def condition_average_ds(ds, *args):
    # xarray groupby functions seem a bit clunky...
    pass


def condition_averages(f_rates, condlist, cond_groups=None):
    
    assert isinstance(condlist, pd.DataFrame)

    if cond_groups is None:
        cond_groups, ic = np.unique(condlist.to_numpy('float64'), axis=0,
                                    return_inverse=True)
        cond_groups = pd.DataFrame(cond_groups, columns=condlist.columns)

    else:
        # cond_groups user-specified
        assert isinstance(cond_groups, pd.DataFrame)
        cond_groups = cond_groups.loc[:, condlist.columns]

        ic = np.full(condlist.shape[0], fill_value=-1, dtype=int)
        for i, c in enumerate(cond_groups.values):
            ic[(condlist == c).all(axis=1)] = i

    
    # nC = cond_groups.shape[0]
    nC = len(cond_groups.index)

    # average across trials of each unique condition
    cond_fr = np.dstack([f_rates[:, ic == c, :].mean(axis=1)
                        for c in range(nC)])
    cond_fr = np.swapaxes(cond_fr, 1, 2)

    # cond_sem = np.array([spike_counts[ic == c, :].std(axis=1) /
    #                      np.sqrt(np.sum(ic == c))
    #                      for c in range(nC)])
    # cond_sem = np.swapaxes(cond_sem, 1, 2)

    return cond_fr, cond_groups



# %% plot functions

def plot_rates(cond_fr):
    
    pass


def plot_psth(df, palette, type=1):

    # TODO
    # code for rasterized spike plots
    # rescale x-axes according to length of time vector
    # show coherences

    if type == 'avg':
        sns.lineplot(data=df,
                     x='heading', y='spike counts',
                     estimator='mean', errorbar='se', err_style='bars',
                     hue=df[condlabels[:-1]].apply(tuple, axis=1))
    elif type == 'time':
        g = sns.relplot(data=df,
                        x='time', y='firing_rate',
                        hue='heading', row='modality', col='align',
                        palette=palette,
                        facet_kws=dict(sharex=False),
                        kind='line', aspect=.75)
