#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 11:41:29 2023

@author: stevenjerjian
"""

# need functions to get unit information for given session


import numpy as np
import pandas as pd
import glob
from pathlib import Path 
import os

# need a function to create unit dataframe, then a decorator to create PSTHs


def loadKSspikes(subject, date, set=1, datapath='/Volumes/homes/fetschlab/data/'):

    subfolder = f'{date}/{subject}{date}_{set}/'
    datapath = f'{datapath}{subject}/{subject}_neuro/'

    filepath = os.path.join(datapath,subfolder)
    npy_files = glob.glob(f'{filepath}*.npy')

    neuro_df = pd.DataFrame()

    npysizes = []
    
    for f in npy_files:
        name = Path(f).name
        name = name.split('.')[0]
    
        if name in ['amplitudes', 'spike_clusters',
                    'spike_templates', 'spike_times']:
            dat = np.squeeze(np.load(f))

            neuro_df[name] = pd.Series(dat)
    

# 
def getSessionInfo():
    pass




if __name__ == '__main__':
    pass
        

