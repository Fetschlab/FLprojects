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
import matplotlib.pyplot as plt


def trial_psth(spiketimes, align_ev, trange, binsize=0.05,
               all_trials=False, normalize=False):
    """
    Inputs:
        spiketimes - 1-D array of spike times
        align_ev   - array of times on each trial (nTr*2 or nTr*2)
        trange     - length 2 array, start and end times wrt align_ev
        binsize    - size of bin size count (default=0.05)
        all_trials - compute for all trials in align_ev, or just trials
        spanned by range of spike times
        normalize  - firing rate (divide by binsize), or just count
                (default=False)

    spiketimes, align_ev, trange, and binsize should all be in the same units
    (default seconds)
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

        if normalize:
            Yt = Yt / binsize
            Y = Y / durs

        return Yt, x, Y, durs, which_ev
