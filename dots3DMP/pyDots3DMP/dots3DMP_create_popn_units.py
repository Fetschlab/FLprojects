# %% imports

import numpy as np
import pandas as pd
from pathlib import Path, PurePath
from datetime import date
import scipy.io as sio
from dataclasses import dataclass, field

# %% define Neuron and Population classes

@dataclass
class Neuron:
    spiketimes: np.ndarray = field(repr=False)
    amps: np.ndarray  = field(repr=False)

    #wfs: np.ndarray = field(repr=False, default_factory=lambda: np.zeros(shape=int, dtype=np.float64))
    #template: np.ndarray = field(repr=False, default_factory=lambda: np.zeros(shape=int, dtype=np.float64))

    temp_amp: float = field(default=0.0)

    nspks: int = field(init=False)

    unit_id: int = field(default=0)
    clus_group: int = field(default=0)
    clus_label: str = field(default='none')
    ch_depth: tuple[int] = field(default=(0,1), 
                                 metadata={'ord':('contact number','contact_depth')})
    
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
    subject: str = 'test'
    set_num: int = 1
    pen_no: int = 1
    grid_xy: tuple[int] = field(default=(np.nan,np.nan))
    grid_type: str = 'odd_sym'

    mdi_depth: int = field(default=0,metadata={'unit':'mm'})
    area: str = 'MSTd'

    device: str = ''
    numchs: int = 32
    
    sr: float = field(default=30000.0, metadata={'unit': 'Hz'})

    units: list = field(default_factory=list, repr=False)
    events: dict = field(default_factory=dict, repr=False)

    def __post_init__(self):
        self.nUnits = len(self.units)

# %%

def build_rec_popn(subject, rec_date, set_num=1, probe_num=1):
    """
    Population is a class instance containing all simulataneously recorded neurons 
    (from one area/probe) and associated stimuli/behavioral events

    """
    datapath='/Volumes/homes/fetschlab/data/'

    session = f'{subject}{rec_date}_{set_num}'
    data_folder = PurePath(datapath, subject, f'{subject}_neuro/', str(rec_date))
    filepath = PurePath(data_folder,session)

    # read info file
    # need to decide how to handle rec_sets with more than one area, 
    # and assign the area information / split the chs and units
    sess_info_file = f'{subject}{rec_date}dots3DMP_info.mat'  
    s = sio.loadmat(PurePath(data_folder,sess_info_file))
    sess_info = s['info']

    # read in events for this set
    set_events_file = f'{subject}{rec_date}dots3DMPevents_{set_num}.mat'
    set_events_path = '/Users/stevenjerjian/Desktop/FetschLab/Analysis/data/lucio_neuro/rec_events'
    events = sio.loadmat(PurePath(set_events_path,set_events_file), simplify_cells=True)
    events = pd.DataFrame.from_dict(events['S'])
                          
    # read in cluster groups (from manual curation)
    cgs = pd.read_csv(PurePath(filepath,'cluster_group.tsv'),sep='\t')

    # read cluster info
    clus_info = pd.read_csv(PurePath(filepath,'cluster_info.tsv'),sep='\t')

    ss = np.squeeze(np.load(PurePath(filepath,'spike_times.npy')))
    sg = np.squeeze(np.load(PurePath(filepath,'spike_clusters.npy')))
    st = np.squeeze(np.load(PurePath(filepath,'spike_templates.npy')))
    sa = np.squeeze(np.load(PurePath(filepath,'amplitudes.npy')))

    # initialize the neural population
    rec_popn = Population(subject=subject,rec_date=rec_date, set_num=set_num, pen_no=sess_info['pen']
                           events=events)
    
    for clus in range(len(cgs)):

        # get info for this cluster
        unit_id = cgs['cluster_id'][clus]
        unit_info = clus_info[clus_info['cluster_id']==unit_id].to_dict('records')[0]
        label = unit_info['group']

        # create Neuron instance for this cluster
        unit = Neuron(spiketimes=ss[sg==unit_id]/rec_set.sr, 
                      amps=sa[sg==unit_id],
                      unit_id=unit_id, clus_label=label,
                      clus_group=get_cluster_group(label),
                      ch_depth=(unit_info['ch'], unit_info['depth']),
                      rec_date=rec_set.rec_date,set_num=rec_set.set_num)
        
        # add it to the population
        rec_popn.units.append(unit)

    return rec_set


def get_cluster_group(clus_label,labels=['NaN','mua','good','noise']):
    """
    return numerical category for cluster type string label
    """

    try:
        return labels.index(clus_label)
    except ValueError:
        return None





# %%

if __name__ == '__main__':


    rec_set = build_rec_popn('lucio', 20220512, set_num=1)


    