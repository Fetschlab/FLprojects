#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 08:37:40 2023

@author: stevenjerjian
"""

import numpy as np
from datetime import date

from dataclasses import dataclass, field
from dots3DMP_FRutils import trial_psth, condition_averages

# %%

@dataclass
class Unit:

    spiketimes: np.ndarray = field(repr=False)
    amps: np.ndarray = field(repr=False)

    clus_id: int = field(default=0)
    rec_date: date = date.today().strftime("%Y%m%d")
    
    def isi(self):
        return np.diff(self.spiketimes)
    
    def isi_hist(self, binsize=0.02):
        x = np.arange(0, 2, binsize)
        return np.histogram(self.spiketimes, x)

    # TODO post_init on waveforms
    
    # TODO simple methods e.g. isi, autocorr, ifr
    # TODO inter-unit methods e.g. cross-corr in spiketimes, similarity
    # def ifr(self):
    # def spike_acorr(self, binsize=0.02, lag=)
    # def spike_xcorr(self, other, binsize=0.02, lag)

@dataclass
class ksUnit(Unit):
    """
    a unit extracted from Kilosort, should have a template with amplitude
    also cluster information, and lab/rig specific info
    """
    # TODO
    # wfs: np.ndarray = field(repr=False, default_factory=lambda: np.zeros(shape=int, dtype=np.float64))
    # template: np.ndarray = field(repr=False, default_factory=lambda: np.zeros(shape=int, dtype=np.float64))

    unique_id: int = field(init=False)
    temp_amp: float = np.nan

    # TODO force clus_group to be 0-3, clus_label to be UN, MU, SU, or noise
    clus_group: int = 0
    clus_label: str = ''
    channel: int = 0
    depth: int = field(default=0, metadata={'unit': 'mm'})
    rec_set: int = 1

    def __post_init__(self):
        self.unique_id = \
            int(f"{self.rec_date}{self.rec_set:02d}{self.clus_id:03d}")

    # TODO add contam pct?

    # TODO fix these up, should they be properties??
    @property
    def trial_psth(self, align_ev, **kwargs):
        return trial_psth(self.spiketimes, align_ev, **kwargs)

    @property
    def condition_psth(self, condlist, cond_groups=None):
        return condition_averages(self.frates, condlist, cond_groups)
        

    #def isi(self, binsize=0.01, max=0.2, plot=False): 
    # return isi hist, violations count
    #def corr(self, binsize=0.01, max=0.2, plot=False):
    #def ifr(self, binsize=0.1, plot=False)
    #def plot_waveforms(self):  
    #def wf_width(self):

    #def summary(self, binsize=0.01, max=0.2, plot=False):
    #  isi, corr, ifr, wf_width 2x2 subplots


@dataclass
class Population:
    """
    data from one recording set
        - metadata
        - list of unit instances, one per recorded unit
        - dictionary/pandas df of task events and times
    """

    rec_date: date = date.today().strftime("%Y%m%d")
    create_date: date = date.today().strftime("%Y%m%d")
    subject: str = ''
    session: str = ''
    rec_set: int = 1
    probe_num: int = 1
    device: str = ''
    
    pen_num: int = 1
    grid_xy: tuple[int] = field(default=(np.nan, np.nan))
    grid_type: str = ''
    area: str = ''

    mdi_depth: int = field(default=0, metadata={'unit': 'mm'})
    chs: list = field(default_factory=list, repr=False)
    sr: float = field(default=30000.0, metadata={'unit': 'Hz'})

    units: list = field(default_factory=list, repr=False)
    events: dict = field(default_factory=dict, repr=False)
