#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 08:37:40 2023

@author: stevenjerjian
"""

import numpy as np
import xarray as xr

import matplotlib.pyplot as plt
from datetime import date

from dataclasses import dataclass, field
from dots3DMP_FRutils import trial_psth, spks_eventplot

# %%

@dataclass
class Unit:

    spiketimes: np.ndarray = field(repr=False)
    amps: np.ndarray = field(repr=False)

    clus_id: int = field(default=0)
    rec_date: date = date.today().strftime("%Y%m%d")

    def __len__(self):
        return len(self.spiketimes)

    def isi(self):
        return np.diff(self.spiketimes)

    def isi_cv(self):
        y = self.isi()
        return np.std(y) / np.mean(y)

    def isi_hist(self, binsize=0.02, max_isi=2):
        x = np.arange(0, max_isi+1e-6, binsize)
        y, _ = np.histogram(self.isi(), x, density=True)
        return y, x[:-1]

    def ifr(self):
        return 1 / self.isi()

    def acf(self, binsize=0.02, maxlag=2):
        ...

    def ccf(self, other, binsize=0.02, maxlag=2):
        ...

    def mean_wf(self):
        ...

    # def summary(self, binsize=0.01, max=0.2, plot=False):
    #  isi, corr, ifr, wf_width 2x2 subplots


@dataclass
class ksUnit(Unit):
    """
    a unit extracted from Kilosort, should have a template with amplitude
    also cluster information, and lab/rig specific info
    """

    # TODO
    # wfs: np.ndarray = field(repr=False, default_factory=lambda:
    #    np.zeros(shape=int, dtype=np.float64))
    # template: np.ndarray = field(repr=False, default_factory=lambda:
    #    np.zeros(shape=int, dtype=np.float64))

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

    # TODO
    # add contam pct?
    # def isi(self, binsize=0.01, max=0.2, plot=False):
    # return isi hist, violations count
    # def corr(self, binsize=0.01, max=0.2, plot=False):
    # def ifr(self, binsize=0.1, plot=False)
    # def plot_waveforms(self):
    # def wf_width(self):

    def raster_plot(self, condlist, row, col,
                    align_ev, trange=np.array([-2, 3])):

        # condlist should be pandas df with conditions
        # align_ev should be np array of same length as condlist

        # good_trs = events['goodtrial'].to_numpy(dtype='bool')
        # condlist = events.loc[good_trs, :]
        # df = condlist[['modality', 'coherence', 'heading']]
        # align_ev = events.loc[good_trs, align_ev].to_numpy(dtype='float64')

        df = condlist
        _, _, df['spks'] = trial_psth(self.spiketimes,
                                      align_ev, trange)

        return spks_eventplot(df, row, col, condlist.iloc[:, -1])


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

    def calc_firing_rates(self, align=['stimOn'], trange=np.array([[-2, 3]]),
                          binsize=0.05, sm_params={},
                          condlabels=['modality', 'coherence', 'heading'],
                          return_Dataset=False):
        """

        Parameters
        ----------
        align : nested lists, optional
            DESCRIPTION. The default is ['stimOn'].
        trange : nested lists, optional
            DESCRIPTION. The default is np.array([[-2, 2]]).
        binsize : FLOAT, optional
            DESCRIPTION. The default is 0.05.
        sm_params : DICT, optional
            DESCRIPTION. The default is {}
        condlabels : LIST, optional
            DESCRIPTION. The default is ['modality', 'coherence', 'heading'].
        return_Dataset : BOOL, optional
            DESCRIPTION. The default is False.

        Raises
        ------
        ValueError
            If there are no 'good trials'

        Returns
        -------
        TYPE
            DESCRIPTION.

        """

        good_trs = self.events['goodtrial'].to_numpy(dtype='bool')
        condlist = self.events[condlabels].loc[good_trs, :]

        rates = []
        tvecs = []
        align_lst = []

        for al, t_r in zip(align, trange):

            if self.events.loc[good_trs, al].isna().all(axis=0).any():
                raise ValueError(al)

            align_ev = self.events.loc[good_trs, al].to_numpy(dtype='float64')

            # get spike counts and relevant t_vec for each unit - thanks chatGPT!
            # trial_psth in list comprehesnsion is going to generate a list of
            # tuples, the zip(*iter) syntax allows us to unpack the tuples into
            # separate variables
            spike_counts, t_vec, _ = \
                zip(*[(trial_psth(unit.spiketimes, align_ev, t_r,
                                  binsize, sm_params))
                      for unit in self.units])

            rates.append(np.asarray(spike_counts))
            tvecs.append(np.asarray(t_vec[0]))

            # previously wanted dict with key as alignment event, but al[0]
            # might not be unique!

            align_lst.append([al]*len(t_vec[0]))

        align_arr = np.concatenate(align_lst)

        if return_Dataset:
            # now construct a dataset
            arr = xr.DataArray(np.concatenate(rates, axis=2),
                               coords=[unitlabels,
                                       condlist.index,
                                       np.concatenate(tvecs)],
                               dims=['unit', 'trial', 'time'])

            cond_coords = {condlabels[c]: ('trial', condlist.iloc[:, c])
                           for c in range(len(condlist.columns))}

            ds = xr.Dataset({'firing_rate': arr},
                            coords={'align_event': ('time', align_arr),
                                    **cond_coords})

            ds.attrs['rec_date'] = self.rec_date
            ds.attrs['area'] = self.area

            return ds

        else:
            # return separate vars
            return rates, tvecs, condlist, align_lst

    def rel_event_times(self, align=['stimOn'], others=['stimOff']):

        good_trs = self.events['goodtrial'].to_numpy(dtype='bool')

        reltimes = {aev: [] for aev in align}
        for aev, oev in zip(align, others):

            if self.events.loc[good_trs, aev].isna().all(axis=0).any():
                raise ValueError(aev)

            align_ev = self.events.loc[good_trs, aev].to_numpy(dtype='float64')
            other_ev = self.events.loc[good_trs, oev].to_numpy(dtype='float64')
            reltimes[aev] = other_ev - align_ev[:, np.newaxis]
        
        return reltimes

# %% functions to run on an individual class instance 
# (in essence class methods, so make them so)
