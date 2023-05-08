#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 08:37:40 2023

@author: stevenjerjian
"""

import numpy as np
import xarray as xr

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from datetime import date

from dataclasses import dataclass, field
import dots3DMP_FRutils as FRutils

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

    # TODO
    # allow this to plot multiple alignments?
    def plot_raster(self, align, condlist, col, hue, titles,
                    trange=np.array([-2, 3]),
                    palette=None, hue_norm=(-12, 12),
                    binsize=0.05, sm_params={}):

        # TODO make color prefs extra kwargs e.g palette
        # TODO add psth plotting to this, as a second row in facetgrid

        # condlist should be pandas df with conditions
        # align_ev should be np array of same length as condlist

        df = condlist[col]
        df[hue] = condlist[hue]
        fr, x, spks = FRutils.trial_psth(self.spiketimes, align, trange,
                                         binsize=binsize, sm_params=sm_params)

        # if col is a list (i.e. 1 or more conditions), create a new grouping
        # column with unique condition groups. otherwise, just use that column
        if isinstance(col, list):
            ic, nC, cond_groups = FRutils.condition_index(df[col])
            df['grp'] = ic
        else:
            nC = len(np.unique(df[col]))
            df['grp'] = df[col]

        assert len(titles) == nC

        # df['fr'] = fr.tolist()
        # df['time'] = [x] * len(df.index)
        # fr_df = df.explode(['fr', 'time'])
        # fr_df = fr_df.reset_index()

        # fr_df[['fr', 'time']] = fr_df[['fr', 'time']].astype('float')
        df['spks'] = spks

        fig, axs = plt.subplots(nrows=2, ncols=nC, figsize=(20, 6))

        # set hue colormap
        uhue = np.unique(df[hue])
        if palette is None:
            palette = sns.diverging_palette(240, 10, n=100, as_cmap=True)

        if isinstance(hue_norm, tuple):
            hue_norm = mcolors.Normalize(vmin=hue_norm[0], vmax=hue_norm[1])

        for c, cond_df in df.groupby('grp'):

            # this time, create groupings based on hue, within cond_df
            ic, nC, cond_groups = FRutils.condition_index(cond_df[[hue]])

            # need to call argsort twice!
            # https://stackoverflow.com/questions/31910407/numpy-argsort-cant-see-whats-wrong !!!
            order = np.argsort(np.argsort(ic)).tolist()

            # TODO further sort within condition by user-specified input e.g. RT

            # get color for each trial, based on heading, convert to list
            colors = palette(hue_norm(cond_df[hue]).data.astype('float'))
            colors = np.split(colors, colors.shape[0], axis=0)

            # ==== raster plot ====
            # no need to re-order all the lists, just use order for lineoffsets
            ax = axs[0, c]
            ax.eventplot(cond_df['spks'], lineoffsets=order, color=colors)
            ax.set_ylim(-0.5, len(cond_df['spks'])+0.5)
            ax.invert_yaxis()
            ax.set_xlim(trange[0], trange[1])
            ax.set_title(titles[c])

            if c == 0:
                ax.set_ylabel('trial')
            else:
                ax.set_yticklabels([])

            # ==== PSTH plot ====
            
            # change this to use matplotlib...

            # c_df = fr_df.loc[fr_df['grp'] == c, :].dropna(axis=0)
            # hue_group, _, _ = FRutils.condition_index(c_df[[hue]], cond_groups)

            # # set colors for unique headings, with same mapping as raster
            # colors = palette(hue_norm(np.unique(c_df[hue])).data.astype('float'))
            # colors = np.split(colors, colors.shape[0], axis=0)

            # ax = axs[1, c]
            # sns.lineplot(x=c_df['time'], y=c_df['fr'],
            #              hue=hue_group, palette=colors,
            #              ax=ax, legend=False, errorbar=None)
            # ax.set_xlim(trange[0], trange[1])

            # if c == 0:
            #     ax.set_ylabel('spikes/sec')
            # else:
            #     ax.set_yticklabels([])
            #     ax.set_ylabel('')

        sm = plt.cm.ScalarMappable(cmap=palette, norm=hue_norm)
        sm.set_array([])
        # cbar_ax = fig.add_axes([0.1, 0.1, 0.05, 0.8])
        cbar = plt.colorbar(sm, ticks=list(uhue))
        cbar.set_label(hue)
        sns.despine()

        plt.show()

        return fig, axs


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

    def calc_firing_rates(self, align_ev='stimOn', trange=np.array([[-2, 3]]),
                          binsize=0.05, sm_params={},
                          condlabels=['modality', 'coherence', 'heading'],
                          return_Dataset=False):

        good_trs = self.events['goodtrial'].to_numpy(dtype='bool')
        condlist = self.events[condlabels].loc[good_trs, :]

        rates = []
        tvecs = []
        align_lst = []

        for al, t_r in zip(align_ev, trange):

            if self.events.loc[good_trs, al].isna().all(axis=0).any():
                raise ValueError(al)

            align = self.events.loc[good_trs, al].to_numpy(dtype='float64')

            # get spike counts and relevant t_vec for each unit
            # trial_psth in list comp is going to generate a list of tuples
            # the zip(*iter) syntax allows us to unpack the tuples into the
            # separate variables
            spike_counts, t_vec, _ = \
                zip(*[(FRutils.trial_psth(unit.spiketimes, align, t_r,
                                          binsize, sm_params))
                      for unit in self.units])

            rates.append(np.asarray(spike_counts))
            tvecs.append(np.asarray(t_vec[0]))

            # considered dict align_lst with key as alignment event, but al[0]
            # might not be unique!

            align_lst.append([al]*len(t_vec[0]))

        align_arr = np.concatenate(align_lst)

        unitlabels = np.array([u.clus_group for u in self.units])
        # unit_ids  = np.array([u.clus_id for u in popn.units])

        if return_Dataset:
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
