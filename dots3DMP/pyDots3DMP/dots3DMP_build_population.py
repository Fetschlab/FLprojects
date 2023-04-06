# %% imports

import numpy as np
import pandas as pd
from pathlib import Path, PurePath
from datetime import date
import scipy.io as sio
from dataclasses import dataclass, field
import pickle

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
    rec_set: int = 1

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

    rec_date: date = date.today().strftime("%Y%m%d")
    create_date: date = date.today().strftime("%Y%m%d")
    subject: str = ''
    session: str = ''
    rec_set: int = 1
    probe_num: int = 1
    pen_num: int = 1
    grid_xy: tuple[int] = field(default=(np.nan, np.nan))
    grid_type: str = ''
    area: str = ''

    mdi_depth: int = field(default=0, metadata={'unit': 'mm'})
    device: str = ''
    chs: list = field(default_factory=list, repr=False)
    sr: float = field(default=30000.0, metadata={'unit': 'Hz'})

    units: list = field(default_factory=list, repr=False)
    events: dict = field(default_factory=dict, repr=False)

# %% build Population class


def build_rec_popn(subject, rec_date, rec_info, data, data_folder):
    """
    Parameters
    ----------
    subject : TYPE
        DESCRIPTION.
    rec_date : TYPE
        DESCRIPTION.
    rec_info : TYPE
        DESCRIPTION.
    data : TYPE
        DESCRIPTION.
    data_folder : TYPE
        DESCRIPTION.

    Returns
    -------
    rec_popn : TYPE
        DESCRIPTION.

    """

    session = f"{subject}{rec_date}_{rec_info['rec_set']}"
    # filepath = PurePath(data_folder, session)

    # # read in cluster groups (from manual curation)
    # cgs = pd.read_csv(PurePath(filepath, 'cluster_group.tsv'), sep='\t')

    # # read cluster info
    # clus_info = pd.read_csv(PurePath(filepath, 'cluster_info.tsv'), sep='\t')

    # ss = np.squeeze(np.load(PurePath(filepath, 'spike_times.npy')))
    # sg = np.squeeze(np.load(PurePath(filepath, 'spike_clusters.npy')))
    # st = np.squeeze(np.load(PurePath(filepath, 'spike_templates.npy')))
    # sa = np.squeeze(np.load(PurePath(filepath, 'amplitudes.npy')))
    


    events = {**data['events'], **data['pldaps']}
    # events = pd.DataFrame({k: pd.Series(v) for k, v in events.items()})
    events = pd.DataFrame.from_dict(events, orient='index')
    events = events.T.convert_dtypes(infer_objects=True)
    events['heading'].loc[np.abs(events['heading']) < 0.01] = 0

    # initialize the neural population
    rec_popn = Population(subject=subject.upper()[0], rec_date=rec_date,
                          session=session, chs=rec_info['chs'], events=events)

    # add all the metadata
    for key in rec_info.keys():
        vars(rec_popn)[key] = rec_info[key]

    # add the units
    if isinstance(data['units']['cluster_id'], int):
        spk_times = np.array(data['units']['spiketimes'])
        unit = Neuron(spiketimes=spk_times, amps=np.empty(spk_times.shape),
                      clus_id=data['units']['cluster_id'],
                      clus_group=data['units']['cluster_type'],
                      ch_depth=(data['units']['ch'], data['units']['depth']),
                      rec_date=rec_popn.rec_date, rec_set=rec_popn.rec_set)
        rec_popn.units.append(unit)

    else:
        nUnits = data['units']['cluster_id'].size
        for u in range(nUnits):
            spk_times = np.array(data['units']['spiketimes'])
            unit = Neuron(spiketimes=spk_times,
                          amps=np.empty(spk_times.shape),
                          clus_id=data['units']['cluster_id'][u],
                          clus_group=data['units']['cluster_type'][u],
                          rec_date=rec_popn.rec_date, rec_set=rec_popn.rec_set)
        rec_popn.units.append(unit)

        # add these in later
        # ch_depth=(data['units']['ch'][u],
        #          data['units']['depth'][u]),

    # only go through clusters in this group of chs (i.e. one probe/area)
    # these_clus_ids = clus_info.loc[clus_info['ch'].isin(rec_info['chs']),
    #                                'cluster_id']
    # cgs = cgs.loc[cgs['cluster_id'].isin(these_clus_ids) &
    #               cgs['group'].isin(groups2keep), :]

    # for clus in cgs.itertuples():

    #     # get info for this cluster
    #     clus_id = clus.cluster_id
    #     unit_info = clus_info[clus_info['cluster_id'] == clus_id].to_dict(
    #         'records')[0]

    #     # create Neuron instance for this cluster
    #     unit = Neuron(spiketimes=ss[sg == clus_id]/rec_popn.sr,
    #                   amps=sa[sg == clus_id],
    #                   clus_id=clus_id, clus_label=clus.group,
    #                   clus_group=get_cluster_group(clus.group),
    #                   ch_depth=(unit_info['ch'], unit_info['depth']),
    #                   rec_date=rec_popn.rec_date, set_num=rec_popn.set_num)

    #     # ...and add it to the population
    #     rec_popn.units.append(unit)

    return rec_popn, events

# %%


def get_cluster_group(clus_label, labels=['unsorted', 'mua', 'good', 'noise']):
    """

    Parameters
    ----------
    clus_label : str
        0, 1, 2, or 3
    labels : list of strings, optional
        kilosort label. The default is ['unsorted', 'mua', 'good', 'noise'].

    Returns
    -------
    TYPE
        DESCRIPTION.

    """

    try:
        return labels.index(clus_label)
    except ValueError:
        return None


# %%


if __name__ == '__main__':

    # TODO - clean this up...
    # subject = input("Enter subject:")
    subject = 'lucio'
    datapath = '/Volumes/homes/fetschlab/data/'
    data_folder = Path(datapath, subject, f'{subject}_neuro/')

    mat_data_file = 'lucio_20220512-20230331_neuralData.mat'
    mat_folder = '/Users/stevenjerjian/Desktop/FetschLab/Analysis/data/lucio_neuro_datasets'

    m = sio.loadmat(PurePath(mat_folder, mat_data_file),
                    simplify_cells=True)
    data = m['dataStruct']

    # https://stackoverflow.com/questions/973473/getting-a-list-of-all-subdirectories-in-the-current-directory
    # rec_dates = [f.parts[-1] for f in data_folder.iterdir()
    #             if f.is_dir() and f.parts[-1][0:2] == '20']

    rec_info_file = '/Users/stevenjerjian/Desktop/FetschLab/Analysis/RecSessionInfo.xlsx'
    rec_info = pd.read_excel(rec_info_file, sheet_name=subject.lower())
    rec_info = rec_info.convert_dtypes(infer_objects=True)
    rec_info.columns = rec_info.columns.str.lower()

    rec_info['date'] = rec_info['date'].apply(lambda x:
                                              x.date().strftime('%Y%m%d'))

    rec_info.rename(columns={'mdi_depth_um': 'mdi_depth',
                             'gt_depth_cm': 'gt_depth',
                             'pen_no_total': 'pen_num',
                             'brain_area': 'area'}, inplace=True)

    rec_df = pd.DataFrame.from_dict(rec_info)
    rec_df['data'] = pd.NA
    
    # TODO - decorator to create a column in rec_df for each paradigm
    par = 'dots3DMP'

    for index, sess in enumerate(data):

        # rec_date = sess['DATE']
        rec_date = str(sess['date'])
        rec_set = sess['rec_set']
        
        print(rec_date, rec_set)

        rec_folder = PurePath(data_folder, rec_date)

        idx = (rec_info['date'] == rec_date) & (rec_info['rec_set'] == set_num)

        
        data = sess['data'][par]




        if 'units' not in data or not any(idx):
            continue

        rec_sess_info = rec_info.loc[idx, :].to_dict('records')[0]
        rec_sess_info['chs'] = np.arange(rec_sess_info['min_ch'],
                                         rec_sess_info['max_ch']+1)
        rec_sess_info['grid_xy'] = (rec_sess_info['grid_x'],
                                    rec_sess_info['grid_y'])

        rec_popn, events = build_rec_popn(subject, rec_date, rec_sess_info,
                                          data, data_folder=rec_folder)

        rec_df['data'][idx] = rec_popn

    filename = PurePath(mat_folder, 'test_data.pkl')
    with open(filename, 'wb') as file:
        pickle.dump(rec_df, file)
        
    # with open('test_data.pkl', 'rb') as file:
    #    this_df = pickle.load(file)