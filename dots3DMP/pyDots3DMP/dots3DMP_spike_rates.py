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

import numpy as np
import pandas as pd
import xarray as xr
from flox.xarray import xarray_reduce
from scipy.stats import norm
from scipy.ndimage import convolve1d, gaussian_filter1d
from pathlib import PurePath

from itertools import repeat, starmap
import pdb

import pickle
from dots3DMP_build_dataset import Population, Neuron

import seaborn as sns
custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(context="notebook", style="ticks", rc=custom_params)

# %%

def trial_psth(spiketimes, align_ev, trange, binsize=0.05,
               smooth_params=('boxcar', 5),
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

        nTr = align_ev.shape[0]

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

        if trange[0] < 0 and trange[1] > 0:
            # ensure that time '0' is in between two bins exactly
            x0 = np.arange(0, tstart_new-1e-3, -binsize)
            x1 = np.arange(0, tend_new, binsize)
            x = np.hstack((x0[::-1, ], x1[1:, ]))
        else:
            x = np.arange(tstart_new, tend_new+1e-3, binsize)

        Yt = np.full([nTr, x.shape[0]-1], np.nan)
        Y = np.full(nTr, np.nan)

        if spiketimes.any():
            itr_start, itr_end = 1, nTr

            if not all_trials:
                itr_start = np.argmin(np.abs(tr_starts - spiketimes[0]))
                itr_end = np.argmin(np.abs(tr_ends - spiketimes[-1]))

            for itr in range(itr_start, itr_end+1):
                spk_inds = np.logical_and(spiketimes >= tr_starts[itr],
                                          spiketimes <= tr_ends[itr])

                Y[itr] = np.sum(spk_inds)

                inds_t = spiketimes[spk_inds] - align_ev[itr, which_ev]
                Yt[itr, :], _ = np.histogram(inds_t, x)

                if which_ev == 0:
                    end_pos = np.argmin(abs(x - tends_new[itr]))
                    Yt[itr, end_pos:] = np.nan
                elif which_ev == 1:
                    start_pos = np.argmin(abs(x - tstarts_new[itr]))
                    Yt[itr, :start_pos] = np.nan

            Yt = smooth_counts(Yt, smooth_params=smooth_params)

            if normalize:
                Yt /= binsize
                Y /= durs

        # shift x values to bin centers
        x = x[:-1] + np.diff(x)/2

        return Yt, x, Y, durs, which_ev



# %%

def smooth_counts(spike_counts, smooth_params=('boxcar', 5)):

    if smooth_params[0] == 'boxcar':

        kw = smooth_params[1]  # width, in bins
        kernel = np.ones(kw) / kw
        smoothed_counts = convolve1d(spike_counts, kernel,
                                     axis=0, mode='reflect')

    elif smooth_params[0] == 'gaussian':

        sigma = smooth_params[1]
        smoothed_counts = gaussian_filter1d(spike_counts, sigma=sigma)

    # TODO causal filter

    return smoothed_counts

# %%

def get_interval_rates(popn, align=['stimOn','stimOff'],
                       trange=np.array([0, 0]),
                       condlabels=['modality', 'coherence', 'heading']):
    pass


def get_aligned_rates(popn, align=['stimOn'], trange=np.array([-2, 2]),
                      binsize=0.05, smooth_params=('boxcar', 5),
                      condlabels=['modality', 'coherence', 'heading'],
                      return_as_Dataset=False):
    # TODO somewhere need to get an update to tvec to be all relative to one event?
    # so that we can plot on the same time axis

    good_trs = popn.events[align].notna().to_numpy(dtype='bool').all(axis=1)
    condlist = popn.events[condlabels].loc[good_trs, :]

    rates = {}
    tvecs = {}
    align_lst = []

    for al, t_r in zip(align, trange):

        # good trials for this alignment event
        align_ev = popn.events.loc[good_trs, al].to_numpy(dtype='float64')

        # thanks chatGPT!
        # get spike counts and t_vec for each unit
        # slice using [:2] since these are the first two outputs of trial_psth
        spike_counts, t_vec = zip(*[(trial_psth(unit.spiketimes, align_ev,
                                                t_r, binsize,
                                                smooth_params)[:2])
                                    for unit in popn.units])

        # previously just appended as a list
        # rates.append(np.asarray(spike_counts))
        # tvecs.append(np.asarray(t_vec[0]))

        # do it as a dictionary, so alignment 'event' is saved as key
        rates[al] = np.asarray(spike_counts)
        tvecs[al] = np.asarray(t_vec[0])

        align_lst.append(np.asarray(list(repeat(al, len(t_vec[0])))))

    align_arr = np.concatenate(align_lst)
    u_labels = np.array([u.clus_group for u in popn.units])

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

        # TODO store some more ds.attrs?
        ds.attrs['rec_date'] = popn.rec_date
        ds.attrs['area'] = popn.area

        return ds

    else:
        # return separate vars
        return rates, u_labels, condlist, tvecs


def condition_average_ds(ds, *args):
    # xarray groupby functions seem a bit clunky...
    pass


def condition_averages(f_rates, condlist, t_vec, cond_groups=None):

    colnames = list(condlist.columns) + ['time', 'firing_rate']
    if f_rates.ndim == 2:
        # TODO
        # single count in interval, f_rates should be units x trial
        # return dataframe?

        data = pd.DataFrame(np.hstack([condlist,
                                       np.expand_dims(t_vec, axis=1),
                                       f_rates]),
                            columns=colnames)

        # grp_df = df.groupby(by=condlabels, axis=0).agg(
        #     meanFR=pd.NamedAgg(column='spike counts', aggfunc="mean"),
        #     semFR=pd.NamedAgg(column='spike counts',
        #                       aggfunc=lambda x: np.std(x) / np.sqrt(len(x))))
        # grp_df = grp_df.reset_index()  # remove multi-level index

    else:
        # TODO test assertions
        # f_rate should be units x trial x time
        # condlist/cond_groups should be a pandas Dataframe
        # output cond_fr will be units x time x condition

        if cond_groups is None:
            cond_groups, ic = np.unique(condlist.to_numpy('float32'), axis=0,
                                        return_inverse=True)
            cond_groups = pd.DataFrame(cond_groups, columns=condlist.columns)

        # nC = cond_groups.shape[0]
        nC = len(cond_groups.index)

        # average across trials of each unique condition
        cond_fr = np.dstack([f_rates[:, ic == c, :].mean(axis=1)
                            for c in range(nC)])

        return cond_fr, cond_groups

        # cond_sem = np.array([spike_counts[ic == c, :].std(axis=0) /
        #                      np.sqrt(np.sum(ic == c))
        #                      for c in range(nC)]).T

    #     cond_groups = np.tile(cond_groups, (len(t_vec), 1))
    #     timestamps = np.repeat(t_vec, nC)

    #     df = pd.DataFrame(np.hstack([cond_groups,
    #                                 np.expand_dims(timestamps, axis=1),
    #                                 cond_mu.T.reshape(-1, 1)]),
    #                       columns=colnames)

    # df.insert(0, 'align', align)

    #return data

# %%

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

if __name__ == '__main__':

    data_folder = '/Users/stevenjerjian/Desktop/FetschLab/Analysis/data/lucio_neuro_datasets'
    filename = PurePath(data_folder, 'test_data.pkl')
    with open(filename, 'rb') as file:
        this_df = pickle.load(file)

    # use a very broad trange here, we can sub-slice it later as needed!
    align = ['stimOn', 'saccOnset']
    trange = np.array([[-2, 0.8], [-0.5, 2]])

    condlabels = ['modality', 'coherence', 'heading']
    par = 'Task'

    binsize = 0.05
    smbinsize = 0.2
    sm_params = ('boxcar', int(smbinsize/binsize))

    data = this_df[this_df[par].notna()][par]

    # res = data.apply(get_aligned_rates,
    #                  args=(align, trange, binsize,
    #                        sm_params, condlabels))

    # firing rate over time across units and trials, per session
    rates, ids, conds, tvecs \
        = zip(*data.apply(get_aligned_rates, args=(align, trange, binsize,
                                                   sm_params, condlabels)))

    # avg firing rate over time across units and each "condition",
    # per session

    cond_frs = []
    cond_groups = [] # may also be different across sessions!

    # outer loop over session
    for fr, cond, tvec in zip(rates, conds, tvecs):
        frs = {}
        # inner loop over different alignments
        for k, f in fr.items():
            frs[k], cg = condition_averages(f, cond, tvec[k])
        cond_frs.append(frs)
        cond_groups.append(cg)

    # TODO plotting helpers
    # TODO plotting and 
    # for pseudopop, should check/force cond_groups to be the same across all units/sessions

    
    #palette = sns.color_palette('PRGn', len(np.unique(condlist['heading'])))
