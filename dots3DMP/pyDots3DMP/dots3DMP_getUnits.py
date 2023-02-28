#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 11:41:29 2023

@author: stevenjerjian
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
from scipy.io import loadmat

from pathlib import Path 

# need a function to create unit dataframe, then a decorator to create PSTHs

datapath = '/Volumes/homes/fetschlab/data/lucio/lucio_neuro/20220308/lucio20220308_1/phy_WC/'
npy_files = glob.glob(f'{datapath}*.npy')


neuro_df = pd.DataFrame()

npysizes = []
for f in npy_files:
    name = Path(f).name
    name = name.split('.')[0]
    
    if name in ['amplitudes', 'spike_clusters',
                'spike_templates', 'spike_times']:
        dat = np.squeeze(np.load(f))
        
        dat_series = pd.Series(dat)
        
        neuro_df[name] = dat_series
    
    



        

