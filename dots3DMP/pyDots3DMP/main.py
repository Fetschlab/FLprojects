#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 08:49:37 2023

@author: stevenjerjian
"""

import numpy as np
import pickle

from dots3DMP_behavior import dots3DMP_create_trial_list

data_folder = '/Users/stevenjerjian/Desktop/FetschLab/Analysis/data/lucio_neuro_datasets'
filename = PurePath(data_folder, 'test_data.pkl')
with open(filename, 'rb') as file:
    this_df = pickle.load(file)

# use a very broad trange here, we can sub-slice it later as needed!
align = ['stimOn', 'saccOnset']
trange = np.array([[-2, 0.8], [-0.5, 2]])

condlabels = ['modality', 'coherence', 'heading', 'delta']

mods = np.array([1, 2, 3])
cohs = np.array([0.2, 0.7])
hdgs = np.array([-12, -6, -3, -1.5, 0, 1.5, 3, 6, 12])
deltas = np.array([0])
nreps = 1

trial_table, ntrials = \
    dots3DMP_create_trial_list(hdgs, mods, cohs, deltas,
                               nreps, shuff=False)
trial_table = np.stack(trial_table, axis=1)
u_conds = pd.DataFrame(trial_table[:, [1, 2, 0, 3]],
                       columns=condlabels)

par = 'Task'

binsize = 0.05
smbinsize = 0.2
sm_params = ('boxcar', int(smbinsize/binsize))

data = this_df[this_df[par].notna()][par]

# firing rate over time across units and trials, per session

# res = data.apply(get_aligned_rates, args=(align, trange, binsize,
#                                           sm_params, condlabels))

rates, avg_rates, ids, conds, tvecs, durs \
    = zip(*data.apply(get_aligned_rates, args=(align, trange, binsize,
                                               sm_params, condlabels)))

# avg firing rate over time across units, each cond, per session
cond_frs = []
ucond_all = []

for fr, cond in zip(rates, conds):
    frs = {}  # OrderedDict?
    # inner loop over different alignments
    for k, f in fr.items():
        frs[k], cg = condition_averages(f, cond, cond_groups=u_conds)
    cond_frs.append(frs)
    ucond_all.append(cg)

# avg rate within interval across units, each cond, per session
avg_cond_frs = []
ucond_all = []
for f, cond in zip(avg_rates, conds):
    f_out, cg = condition_averages(f, cond, cond_groups=u_conds)
    avg_cond_frs.append(f_out)
    ucond_all.append(cg)
    

# TODO plotting helpers
# TODO plotting and 
# for pseudopop, should check/force cond_groups to be the same across all units/sessions


#palette = sns.color_palette('PRGn', len(np.unique(condlist['heading'])))
