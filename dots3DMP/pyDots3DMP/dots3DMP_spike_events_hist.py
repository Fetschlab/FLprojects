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

import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import norm
from scipy.ndimage import convolve1d
from pathlib import PurePath

import pickle
from dots3DMP_build_population import Population, Neuron

import seaborn as sns
custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(context="notebook", style="ticks", rc=custom_params)


def trial_psth(spiketimes, align_ev, trange, binsize=0.05, all_trials=False):
    """

    Parameters
    ----------
    spiketimes : 1-D numpy array
        spike times
    align_ev : numpy array
        nTr x 1 or nTr x 2 array with event times (in same units as spiketimes)
        spike times for psth are always aligned to first column
    trange : numpy array
        start and end times relative to align_ev (in same units as spiketimes)
    binsize : int, optional
        size of bins for spike counts (same units). The default is 0.05.
    all_trials : boolean, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    Yt : numpy array
        DESCRIPTION.
    x : 1-D numpy array
        DESCRIPTION.
    Y : numpy array
        DESCRIPTION.
    durs : numpy array
        DESCRIPTION.
    which_ev : int
        DESCRIPTION.

    """

    nTr = align_ev.shape[0]

    if nTr == 0:
        return
    else:
        if align_ev.ndim == 2:
            align_ev = np.sort(align_ev, axis=1)
            ev_order = np.argsort(align_ev, axis=1)
            which_ev = ev_order[0]  # align to this event
        else:
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

        itr_start = 1
        itr_end = nTr

        if not all_trials and spiketimes.any():
            itr_start = np.argmin(np.abs(tr_starts - spiketimes[0]))
            itr_end = np.argmin(np.abs(tr_ends - spiketimes[-1]))

        # compute 'corrected' tStart and tEnd
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
            # ensure that there 0 time is in between two bins exactly
            x0 = np.arange(0, tstart_new-1e-3, -binsize)
            x1 = np.arange(0, tend_new, binsize)
            x = np.hstack((x0[::-1, ], x1[1:, ]))
        else:
            x = np.arange(tstart_new, tend_new+1e-3, binsize)

        Yt = np.full([nTr, x.shape[0]-1], np.nan)
        Y = np.full([nTr, 1], np.nan)

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

        # shift x values to bin centers
        x = x[:-1] + np.diff(x)/2

        return Yt, x, Y, durs, which_ev


def condition_average(spike_counts, conds, condlabels, t_vec, align='',
                      normalize=False):

    colnames = condlabels + ['time', 'firing_rate']
    if spike_counts.ndim == 1:
        if normalize:
            spike_counts /= t_vec

        data = pd.DataFrame(np.hstack([conds,
                                     np.expand_dims(t_vec, axis=1),
                                     spike_counts]),
                          columns=colnames)

        # grp_df = df.groupby(by=condlabels, axis=0).agg(
        #     meanFR=pd.NamedAgg(column='spike counts', aggfunc="mean"),
        #     semFR=pd.NamedAgg(column='spike counts',
        #                       aggfunc=lambda x: np.std(x) / np.sqrt(len(x))))
        # grp_df = grp_df.reset_index()  # remove multi-level index

    else:
        cond_groups, ic = np.unique(conds.to_numpy('float32'), axis=0,
                                    return_inverse=True)
        nC = cond_groups.shape[0]

        if normalize:
            spike_counts /= np.diff(t_vec[:2])

        cond_mu = np.array([spike_counts[ic == c, :].mean(axis=0)
                            for c in range(nC)]).T

        # cond_sem = np.array([spike_counts[ic == c, :].std(axis=0) /
        #                      np.sqrt(np.sum(ic == c))
        #                      for c in range(nC)]).T

        cond_groups = np.tile(cond_groups, (len(t_vec), 1))
        timestamps = np.repeat(t_vec, nC)

        df = pd.DataFrame(np.hstack([cond_groups,
                                    np.expand_dims(timestamps, axis=1),
                                    cond_mu.T.reshape(-1, 1)]),
                          columns=colnames)

    df.insert(0, 'align', align)

    return data


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
        


def smooth_psth(spike_counts, method='boxcar', **kwargs):

    if method == 'boxcar':

        kw = kwargs['width']
        kernel = np.ones(kw) / kw

    elif method == 'gaussian':

        mu = kwargs['mu']
        sigma = kwargs['sigma']
        kw = kwargs['width']

        kernel_size = int(kw * sigma)  # choose kernel size based on sigma
        kernel = norm.pdf(np.arange(-kernel_size, kernel_size+1),
                          loc=mu, scale=sigma)
        kernel /= kernel.sum()  # normalize kernel

    smoothed_counts = convolve1d(spike_counts, kernel, axis=0, mode='reflect')
    return smoothed_counts, kernel


if __name__ == '__main__':

    data_folder = '/Users/stevenjerjian/Desktop/FetschLab/Analysis/data/lucio_neuro_datasets'
    filename = PurePath(data_folder, 'test_data.pkl')
    with open(filename, 'rb') as file:
        this_df = pickle.load(file)

    popn = this_df.iloc[29]['data']

    spiketimes = popn.units[0].spiketimes

    align = ['stimOn', 'saccOnset']
    trange = np.array([[-2, 0.8], [-0.5, 2]])

    condlabels = ['modality', 'coherence', 'heading']

    good_trs = popn.events[align].notna().to_numpy(dtype='bool').all(axis=1)

    condlist = popn.events[condlabels].loc[good_trs, :]

    palette = sns.color_palette('PRGn', len(np.unique(condlist['heading'])))

    all_df = pd.DataFrame()

    for a, t_r in list(zip(align, trange)):
        align_ev = popn.events.loc[good_trs, a].to_numpy(dtype='float64')

        spike_counts, t_vec, _, _, _ = trial_psth(spiketimes, align_ev, t_r,
                                                  binsize=0.1)

        df = condition_average(spike_counts, condlist, condlabels, t_vec,
                               align=a)

        all_df = pd.concat([all_df, df], axis=0, ignore_index=True)
            
        plot_psth(all_df, palette=palette, type='time')