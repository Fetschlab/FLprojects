#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  6 18:39:16 2023

@author: stevenjerjian
"""

# plot psths for all units

import numpy as np
from pathlib import PurePath
import pandas as pd

# custom imports
from NeuralDataClasses import Population, ksUnit, Unit
from dots3DMP_behavior import dots3DMP_create_trial_list
import dots3DMP_FRutils as FRutils

data_folder = '/Users/stevenjerjian/Desktop/FetschLab/Analysis/data'
filename = PurePath(data_folder, 'lucio_neuro_datasets',
                    'lucio_20220512-20230411_neuralData.pkl')

with open(filename, 'rb') as file:
    # this_df = pkl.load(file)
    this_df = pd.read_pickle(file)

par = ['Tuning', 'Task']
data = this_df[this_df[par].notna().all(axis=1)][par]

# %%

condlabels = ['modality', 'coherenceInd', 'heading', 'delta']
mods = np.array([1, 2, 3])
cohs = np.array([1, 2])
deltas = np.array([0])

# hdgs = np.array([-60, -45, -25, -22.5, -12, 0, 12, 22.5, 25, 45, 60])
hdgs = np.array([-12, -6, -3, -1.5, 0, 1.5, 3, 6, 12])
tr_tab, _ = dots3DMP_create_trial_list(hdgs, mods, cohs, deltas,
                                       1, shuff=False)
tr_tab.columns = condlabels


# %% example raster and psth

sess, unit = 0, 0  # session, unit

events = data['Task'][sess].events
good_trs = events['goodtrial'].to_numpy(dtype='bool')
condlist = events[condlabels].loc[good_trs, :]

align_ev = 'stimOn'
align = events.loc[good_trs, align_ev].to_numpy(dtype='float64')


titles = ['Ves', 'Vis L', 'Vis H', 'Comb L', 'Comb H']

data['Task'][sess].units[unit].plot_raster(align, condlist,
                                     ['modality', 'coherenceInd'], 'heading',
                                     titles,
                                     trange=np.array([-2, 3]),
                                     binsize=0.05, sm_params={})

# %% trial firing rates

# get all unit firing rates across trials/conds, for tuning and task


binsize = 0.05


sm_params = {}
# sm_params = {'type': 'boxcar', 'binsize': binsize, 'width': 0.2}
# sm_params = {'type': 'gaussian', 'binsize': binsize,
#             'width': 0.5, 'sigma': 0.05}

align_ev = ['stimOn', 'stimOff']
trange = np.array([[-1, 1.25], [-0.3, 1]])

rates, tvecs, conds, _ = \
    zip(*data['Task'].apply(lambda x: x.calc_firing_rates(align_ev, trange,
                                                          binsize, sm_params,
                                                          condlabels)))

rates_cat, _ = FRutils.concat_aligned_rates(rates)


# %% cond avg in task, PSTH plotting

# cond_frs, cond_groups = [], []
# for f_in, cond in zip(rates_cat, conds):

#     # avg firing rate over time across units, each cond, per session
#     f_out, _, cg = FRutils.condition_averages(f_in, cond, cond_groups=tr_tab)
#     cond_frs.append(f_out)
#     cond_groups.append(cg)

# # stack 'em up. all units x conditions x time
# cond_frs_stacked = np.vstack(cond_frs)


# %%

