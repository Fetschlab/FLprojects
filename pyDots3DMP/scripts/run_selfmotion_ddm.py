# %% ----------------------------------------------------------------
from datetime import datetime
import logging

import numpy as np

# temp hack to import from behavior module...should setup pydots3dmp as a package
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from behavior.selfmotionddm import SelfMotionDDM
from behavior.utils import data_cleanup

# %% ----------------------------------------------------------------
# Set up loggers

def setup_logger(log_file_path=None): 
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)

    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
        handler.close()

    # # Create handlers
    c_handler = logging.StreamHandler()
    c_handler.setLevel(logging.INFO) 
    c_format = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    c_handler.setFormatter(c_format)
    logger.addHandler(c_handler) 

    if log_file_path:
        f_handler = logging.FileHandler(log_file_path)
        f_handler.setLevel(logging.INFO) 
        f_format = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        f_handler.setFormatter(f_format)
        logger.addHandler(f_handler)

    return logger


# %% ----------------

# co-ordinate vectors for accumulator process
grid_vec = np.arange(-3, 0, 0.01)
time_vec = np.arange(0, 2, 0.05)

# NOTE on kmult:
# three-element list --> independent sensitivites for ves, low coh vis and high coh vis
# two-element list --> independent sensitivites for ves and vis. low coh and high coh proportional based on values
# one-element --> low coh and high coh proportional based on values, mean(low, high) = ves
#
# all other init_params values are lists, but can be:
# single element - same value for all modalities
# three elements - independent value per modality    

init_params = {
    'kmult': [1, 1],            # ves, vis sensitivites
    'bound': [1, 1, 1],         # ves, vis, comb bounds
    'non_dec_time': [0.3],      # non-decision time (secs)
    'wager_thr': [1, 1, 1],     # log odds threshold for high bets
    'wager_alpha': [0.05],      # base rate of high bets
}

# initialize DDM object for fitting
ddm = SelfMotionDDM(
    grid_vec=grid_vec,
    tvec=time_vec,
    **init_params, 
    stim_scaling=True,      # scale ves/vis according to acc/vel signals
    return_wager=True       # whether to compute pdfs and log odds maps, and return wagers
    )


# Get data for fitting
datafilepath = "/Users/stevenjerjian/Desktop/Academia/FetschLab/PLDAPS_data/dataStructs/lucio_20220512-20230606.csv"
data = data_cleanup(datafilepath)  # this is a hack function to quickly load and clean the data, could be improved/generalized with options

# remove one-target wager trials
data = data.loc[data['oneTargConf'] == 0, :]  

# split data into X (conditions) and y (outcomes)
X_all = data[['modality', 'coherence', 'delta', 'heading']]
y_all = data[['choice', 'PDW', 'RT']].astype({'choice': 'int', 'PDW': 'int'})

# fit only delta=0 trials
delta0 = X_all['delta'] == 0

X = X_all[delta0].reset_index(drop=True)
y = y_all[delta0].reset_index(drop=True)
y['RT'] += 0.3 # to reverse the downshift from motion platform latency

# Run the fitting
ddm.fit(X, y)

# ddm.simulate(X, n_samples=10, sample_dvs=True, seed=1)

y_pred, y_pred_samp = ddm.predict(X, n_samples=1, cache_accumulators=True, seed=1)
