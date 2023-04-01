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

    clus_id: int = field(default=0)
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
    create_date: date = date.today().strftime("%Y%m%d")
    subject: str = 'test'
    set_num: int = 1
    pen_no: int = 1
    grid_xy: tuple[int] = field(default=(np.nan,np.nan))
    grid_type: str = 'odd_sym'

    mdi_depth: int = field(default=0,metadata={'unit':'mm'})
    device: str = ''
    chs: list = field(default_factory=list, repr=False)
    
    sr: float = field(default=30000.0, metadata={'unit': 'Hz'})

    units: list = field(default_factory=list, repr=False)
    events: dict = field(default_factory=dict, repr=False)

    def __post_init__(self):
        self.nUnits = len(self.units)

# %%

def build_rec_popn(subject, rec_date, sess_info, set_num=1, probe_num=1,
                   groups2keep = ['good','mua'], 
                   datapath = '/Volumes/homes/fetschlab/data/'):
    """
    Population is a class instance containing all simulataneously recorded neurons 
    (from one area/probe) and associated stimuli/behavioral events

    """
    session = f'{subject}{rec_date}_{set_num}'
    filepath = PurePath(data_folder,session)

    # read in events for this set
    set_events_file = f'{subject}{rec_date}dots3DMPevents_{set_num}.mat'
    set_events_path = '/Users/stevenjerjian/Desktop/FetschLab/Analysis/data/lucio_neuro/rec_events'
    events = sio.loadmat(PurePath(set_events_path,set_events_file), simplify_cells=True)
    events = pd.DataFrame.from_dict(events['S'])

    mdi_depth = np.mean(sess_info['depths'][sess_info['rec_group']==set_num])
    paradigms = sess_info['par'][sess_info['rec_group']==set_num]
    
    if sess_info['gridxy'].ndim == 1:
            gridxy = sess_info['gridxy']
            gridtype = sess_info['gridtype']
            device = sess_info['probe_ID']
            chs = sess_info['chanlist']
    else:
        gridxy = sess_info['gridxy'][probe_num]
        gridtype = sess_info['gridtype'][probe_num]
        device = sess_info['probeID'][probe_num]
        chs = sess_info['chanlist'][probe_num]        
                          
    # read in cluster groups (from manual curation)
    cgs = pd.read_csv(PurePath(filepath,'cluster_group.tsv'),sep='\t')

    # read cluster info
    clus_info = pd.read_csv(PurePath(filepath,'cluster_info.tsv'),sep='\t')

    ss = np.squeeze(np.load(PurePath(filepath,'spike_times.npy')))
    sg = np.squeeze(np.load(PurePath(filepath,'spike_clusters.npy')))
    st = np.squeeze(np.load(PurePath(filepath,'spike_templates.npy')))
    sa = np.squeeze(np.load(PurePath(filepath,'amplitudes.npy')))

    # initialize the neural population
    rec_popn = Population(subject=subject,rec_date=rec_date, set_num=set_num, pen_no=sess_info['pen'],
                          grid_xy=gridxy, grid_type=gridtype, device=device, chs=chs, events=events)
    

    # only go through clusters in this group of chs (i.e. one probe/area)
    these_clus_ids = clus_info.loc[clus_info['ch'].isin(chs),'cluster_id']
    cgs = cgs.loc[cgs['cluster_id'].isin(these_clus_ids) & cgs['group'].isin(groups2keep),:]
    
    for clus in cgs.itertuples():

        # get info for this cluster
        clus_id = clus.cluster_id
        unit_info = clus_info[clus_info['cluster_id']==clus_id].to_dict('records')[0]

        # create Neuron instance for this cluster
        unit = Neuron(spiketimes=ss[sg==clus_id]/rec_popn.sr, 
                      amps=sa[sg==clus_id],
                      clus_id=clus_id, clus_label=clus.group,
                      clus_group=get_cluster_group(clus.group),
                      ch_depth=(unit_info['ch'], unit_info['depth']),
                      rec_date=rec_popn.rec_date,set_num=rec_popn.set_num)
        
        # ...and add it to the population
        rec_popn.units.append(unit)

    return rec_popn

# %% helpers

def get_cluster_group(clus_label, labels=['NaN','mua','good','noise']):
    """
    return numerical category for cluster type string label
    """

    try:
        return labels.index(clus_label)
    except ValueError:
        return None


def trial_psth(spiketimes, align_ev, trange, binsize=0.05, all_trials=False, normalize=False):
    pass

    nTr = align_ev.shape[0]

    if nTr==0:
        return
    else:
        
        if align_ev.ndim==2:
            align_ev = np.sort(align_ev,axis=1)
            ev_order = np.argsort(align_ev,axis=1)
            which_ev = ev_order[0] # align to this event
        else:
            align_ev[:,1] = align_ev[:,0]
            which_ev = 0


        tr_starts = align_ev[:,0] + trange[0]
        tr_ends = align_ev[:,1] + trange[1]

        durs = tr_ends - tr_starts

        itr_start = 1
        itr_end = nTr

        if not all_trials and spiketimes.any():
            itr_start = np.argmin(np.abs(tr_starts - spiketimes[0]))
            itr_end = np.argmin(np.abs(tr_ends - spiketimes[-1]))
        
        







# %%

if __name__ == '__main__':

    #subject = input("Enter subject:")

    from dateutil.parser import parse

    subject = 'lucio'
    datapath='/Volumes/homes/fetschlab/data/'
    data_folder = Path(datapath, subject, f'{subject}_neuro/')

    # https://stackoverflow.com/questions/973473/getting-a-list-of-all-subdirectories-in-the-current-directory
    rec_dates = [f.parts[-1] for f in data_folder.iterdir() 
        if f.is_dir() and f.parts[-1][0:2]=='20']

    for rec_date in rec_dates:
        rec_folder = PurePath(datapath, subject, f'{subject}_neuro/', str(rec_date))

       # read info file for this date
        sess_info_file = f'{subject}{rec_date}dots3DMP_info.mat'  
        s = sio.loadmat(PurePath(rec_folder,sess_info_file), simplify_cells=True)
        sess_info = s['info']

        if rec_date <= 20230330:
            sess_info['numch'] = sess_info['chanlist'][-1]

        nProbes = sess_info['gridxy'].ndim
        nSets = np.max(sess_info['rec_group'])

        for s in range(nSets):
            for p in range(nProbes):
                rec_popn = build_rec_popn('lucio', rec_date, sess_info,
                                         set_num=s+1, probe_num=p)


    