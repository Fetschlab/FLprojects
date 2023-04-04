# %% imports

import numpy as np
import pandas as pd
from pathlib import Path, PurePath
from datetime import date
import scipy.io as sio
from dataclasses import dataclass, field

import pdb


# %% define Neuron and Population classes


@dataclass
class Neuron:
    spiketimes: np.ndarray = field(repr=False)
    amps: np.ndarray = field(repr=False)

    #wfs: np.ndarray = field(repr=False, default_factory=lambda: np.zeros(shape=int, dtype=np.float64))
    #template: np.ndarray = field(repr=False, default_factory=lambda: np.zeros(shape=int, dtype=np.float64))

    temp_amp: float = field(default=0.0)

    nspks: int = field(init=False)

    clus_id: int = field(default=0)
    clus_group: int = field(default=0)
    clus_label: str = field(default='none')
    ch_depth: tuple[int] = field(default=(0, 1),
                                 metadata={'ord': ('contact number',
                                                   'contact_depth')})
    rec_date: date = date.today().strftime("%Y%m%d")
    set_num: int = 1

    def __post_init__(self):
        self.nspks = len(self.spiketimes)

    # static methods?

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
    rec_date: date
    create_date: date = date.today().strftime("%Y%m%d")
    subject: str = ''
    session: str = ''
    set_num: int = 1
    probe_num: int = 1
    pen_num: int = 1
    grid_xy: tuple[int] = field(default=(np.nan, np.nan))
    grid_type: str = 'odd_sym'
    area: str = ''

    mdi_depth: int = field(default=0, metadata={'unit': 'mm'})
    device: str = ''
    chs: list = field(default_factory=list, repr=False)
    sr: float = field(default=30000.0, metadata={'unit': 'Hz'})

    units: list = field(default_factory=list, repr=False)
    events: dict = field(default_factory=dict, repr=False)

    def __post_init__(self):
        self.nUnits = len(self.units)

# %% build Population class


def build_rec_popn(subject, rec_date, set_info,
                   groups2keep=['good', 'mua'],
                   data_folder='/Volumes/homes/fetschlab/data/'):
    """
    Population is a class instance containing all simultaneous recorded units
    (from one area/probe) and associated stimuli/behavioral events

    rec_info - metadata
    chs- chs for this population

    """
    session = f"{subject}{rec_date}_{set_info['set_num']}"
    filepath = PurePath(data_folder, session)

    # read in events for this set
    events_file = f"{subject}{rec_date}dots3DMPevents_{set_info['set_num']}.mat"
    events = sio.loadmat(PurePath(data_folder, 'rec_events', events_file),
                         simplify_cells=True)
    events = pd.DataFrame.from_dict(events['S'])

    # read in cluster groups (from manual curation)
    cgs = pd.read_csv(PurePath(filepath, 'cluster_group.tsv'), sep='\t')

    # read cluster info
    clus_info = pd.read_csv(PurePath(filepath, 'cluster_info.tsv'), sep='\t')

    ss = np.squeeze(np.load(PurePath(filepath, 'spike_times.npy')))
    sg = np.squeeze(np.load(PurePath(filepath, 'spike_clusters.npy')))
    # st = np.squeeze(np.load(PurePath(filepath, 'spike_templates.npy')))
    sa = np.squeeze(np.load(PurePath(filepath, 'amplitudes.npy')))

    # initialize the neural population
    rec_popn = Population(subject=subject.upper()[0], rec_date=rec_date,
                          session=session, chs=set_info['chs'], events=events)

    for key in set_info.keys():
        vars(rec_popn)[key] = set_info[key]

    # only go through clusters in this group of chs (i.e. one probe/area)
    these_clus_ids = clus_info.loc[clus_info['ch'].isin(set_info['chs']),
                                   'cluster_id']
    cgs = cgs.loc[cgs['cluster_id'].isin(these_clus_ids) &
                  cgs['group'].isin(groups2keep), :]

    for clus in cgs.itertuples():

        # get info for this cluster
        clus_id = clus.cluster_id
        unit_info = clus_info[clus_info['cluster_id'] == clus_id].to_dict(
            'records')[0]

        # create Neuron instance for this cluster
        unit = Neuron(spiketimes=ss[sg == clus_id]/rec_popn.sr,
                      amps=sa[sg == clus_id],
                      clus_id=clus_id, clus_label=clus.group,
                      clus_group=get_cluster_group(clus.group),
                      ch_depth=(unit_info['ch'], unit_info['depth']),
                      rec_date=rec_popn.rec_date, set_num=rec_popn.set_num)

        # ...and add it to the population
        rec_popn.units.append(unit)

    return rec_popn

# %% helpers


def get_cluster_group(clus_label, labels=['NaN', 'mua', 'good', 'noise']):
    """
    return numerical category for cluster type string label
    """

    try:
        return labels.index(clus_label)
    except ValueError:
        return None


# def get_set_info(sess_info, set_num, probe_num, probe_nChs):

#     rec_info = {}

#     rec_info['set_num'] = set_num
#     rec_info['probe_num'] = probe_num

#     ch_inds = np.logical_and(sess_info['chanlist'] >= probe_nChs[probe_num],
#                              sess_info['chanlist'] <= probe_nChs[probe_num+1])
#     rec_info['chs'] = sess_info['chanlist'][ch_inds]
#     rec_info['gridtype'] = sess_info['gridtype']

#     # paradigms = sess_info['par'][sess_info['rec_group'] == set_num]

#     if sess_info['numch'].size == 1:
#         rec_info['grid_xy'] = sess_info['gridxy']
#         rec_info['device'] = sess_info['probe_ID']
#         rec_info['mdi_depth'] = np.mean(sess_info['depths']
#                                         [sess_info['rec_group'] == set_num])
#     else:
#         gridxy = np.concatenate(sess_info['gridxy']).reshape(-1, 2)
#         rec_info['grid_xy'] = gridxy[probe_num]
#         rec_info['device'] = sess_info['probe_ID'][probe_num]
#         rec_info['mdi_depth'] = np.mean(sess_info['depths'][probe_num]
#                                         [sess_info['rec_group'] == set_num])

#     return rec_info

def get_set_info(sess_info):

    rec_info = {}

    rec_info['set_num'] = sess.REC_SET
    rec_info['probe_num'] = sess.PROBE_NUM
    rec_info['chs'] = np.arange(sess.MIN_CH, sess.MAX_CH+1)
    rec_info['area'] = sess.BRAIN_AREA
    rec_info['gt_depth_cm'] = sess.GT_DEPTH_CM
    rec_info['pen_num'] = sess.PEN_NO_TOTAL
    rec_info['mdi_depth'] = sess.MDI_DEPTH_UM
    rec_info['device'] = (sess.PROBE_TYPE, sess.PROBE_ID)
    rec_info['grid_xy'] = (sess.GRID_X, sess.GRID_Y)
    rec_info['grid_type'] = sess.GRID_TYPE

    return rec_info

# %%


if __name__ == '__main__':

    import pickle
    import pandas as pd
    from collections import OrderedDict
    from datetime import datetime

    # subject = input("Enter subject:")
    subject = 'lucio'
    datapath = '/Volumes/homes/fetschlab/data/'
    data_folder = Path(datapath, subject, f'{subject}_neuro/')

    # https://stackoverflow.com/questions/973473/getting-a-list-of-all-subdirectories-in-the-current-directory
    # rec_dates = [f.parts[-1] for f in data_folder.iterdir()
    #             if f.is_dir() and f.parts[-1][0:2] == '20']

    rec_info_file = '/Users/stevenjerjian/Desktop/FetschLab/Analysis/RecSessionInfo.xlsx'
    rec_info = pd.read_excel(rec_info_file, sheet_name=subject.lower())

    rec_info['DATE'] = rec_info['DATE'].apply(lambda x:
                                              x.date().strftime('%Y%m%d'))

    rec_data = OrderedDict()

    for index, sess in rec_info.iterrows():

        rec_date = sess['DATE']
        print(rec_date)

        rec_folder = PurePath(data_folder, str(rec_date))


        set_info = get_set_info(sess)

# this is all unnecessary now
        # read info file for this date
        # sess_info_file = f'{subject}{rec_date}dots3DMP_info.mat'
        # s = sio.loadmat(PurePath(rec_folder, sess_info_file),
        #                 simplify_cells=True)
        # sess_info = s['info']

        # if np.array(sess_info['chanlist']).size == 1:
        #    continue

        # if datetime.strptime(rec_date, '%Y%m%d') <= datetime.strptime(
        #       '20230330', '%Y%m%d'):
        #    sess_info['numch'] = np.array(sess_info['chanlist'].size)

        # nProbes = sess_info['gridxy'].ndim
        # nSets = np.max(sess_info['rec_group'])

        # probe_nChs = np.cumsum(sess_info['numch'])
        # probe_nChs = np.insert(probe_nChs, 0, 0)

        rec_popn = build_rec_popn(subject, rec_date, set_info,
                                  data_folder=rec_folder)

        # replace this with area reference?
        key = f'{rec_popn.session}_p{rec_popn.probe_num}'
        rec_data[key] = rec_popn

    with open('alldata.pkl', 'wb') as file:
        pickle.dump(rec_data, file)
