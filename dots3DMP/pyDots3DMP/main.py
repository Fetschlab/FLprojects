#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 08:49:37 2023

@author: stevenjerjian
"""

# %% imports and data load

import numpy as np
import pandas as pd
from pathlib import PurePath
import pickle as pkl

# custom imports
from NeuralDataClasses import Population, ksUnit, Unit
from dots3DMP_behavior import dots3DMP_create_trial_list
import dots3DMP_FRutils as FRutils
import tuning_utils as tuning

data_folder = '/Users/stevenjerjian/Desktop/FetschLab/Analysis/data'
filename = PurePath(data_folder, 'lucio_neuro_datasets',
                    'lucio_20220512-20230411_neuralData.pkl')

with open(filename, 'rb') as file:
    this_df = pkl.load(file)

# %% define conditions and data of interest

# define conditions of interest
condlabels = ['modality', 'coherenceInd', 'heading', 'delta']
mods = np.array([1, 2, 3])
cohs = np.array([1, 2])
# hdgs = np.array([-60, -45, -25, -22.5, -12, 0, 12, 22.5, 25, 45, 60])
hdgs = np.array([-12, -6, -3, -1.5, 0, 1.5, 3, 6, 12])
deltas = np.array([0])

trial_table, ntrials = dots3DMP_create_trial_list(hdgs, mods, cohs, deltas, 1,
                                                  shuff=False)

# parameters for extracting firing rates
par = 'Task'
align = [['stimOn', 'stimOn'], ['stimOn', 'stimOff']]
trange = np.array([[-0.5, 0], [0, 0]])

# time-resolved
binsize = 0.05
smbinsize = 0.2
sm_params = ('boxcar', int(smbinsize/binsize))

# aggregated spike rates in intervals
binsize = 0
sm_params = ('boxcar', 0)

data = this_df[this_df[par].notna()][par]

# %% get all unit firing rates across trials/conds


# firing rate over time across units and trials, per session
rates, unitlabels, tvecs, align_lst, conds = zip(*data.apply(
    FRutils.get_aligned_rates, args=(align, trange, binsize, sm_params,
                                     condlabels)))

# avg firing rate over time across units, each cond, per session
cond_frs, ucond_all = [], []
for f_in, cond in zip(rates, conds):
    # for k, f in f_in.items() # if using dict for alignment events
    f_out, _, cg = FRutils.condition_averages(f_in, cond,
                                              cond_groups=trial_table)
    cond_frs.append(f_out)
    ucond_all.append(cg)
au_frs = np.vstack(cond_frs)


# %% tuning curves

# this needs work

xhdgs = np.linspace(hdgs.min(), hdgs.max(), 200)
nUnits = au_frs.shape[0]

# for u in range(nUnits):
#     for m in mods:
#         for c in cohs:

#             temp = au_frs[u, (u_conds['modality'] == m) &
#                              (u_conds['coherenceInd'] == c), :]
#             inds = ~(np.any(np.isnan(temp), axis=1))
#             temp = temp[inds, :]
#             thesehdgs = hdgs[inds]

#             y_pred, popt, p0, perr = \
#                 zip(*np.apply_along_axis(tuning.fit_predict_vonMises,
#                                          axis=0, arr=temp,
#                                          x=thesehdgs, x_pred=xhdgs).T)
#             y_pred = np.vstack(y_pred)

# %% 'noise' correlations

# within recording areas

# TODO use map to reduce loops here
cond_corrs, cond_pvals = [], []
for f_in, cond in zip(rates, conds):
    corrs_t, pvals_t = [], []
    for f_t in f_in:
        corrs, pvals = FRutils.rsc_within(f_t, cond,
                                          cond_groups=trial_table)
        corrs_t.append(corrs)
        pvals_t.append(pvals)
    cond_corrs.append(corrs_t)
    cond_pvals.append(pvals_t)

# across areas




