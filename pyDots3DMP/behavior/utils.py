
import numpy as np
import pandas as pd
from pathlib import Path, PurePath
from typing import Literal, Optional, Union

from itertools import product

# %% basic stats functions
def prop_se(x):
    """standard error of proportion"""
    return np.sqrt((np.mean(x)*(1-np.mean(x))) / len(x))


def cont_se(x):
    """standard error of continuous variable"""
    return np.std(x) / np.sqrt(len(x))


def gaus(x, ampl, mu, sigma, bsln):
    return ampl * np.exp(-(x-mu)**2 / (2*sigma**2)) + bsln


def prop_se_minmax(x):
    se = prop_se(x)
    return (np.mean(x)-se, np.mean(x)+se)


def log_lik_bin(y, y_hat):
    # return np.sum(y * np.log(y_hat) + (1 - y) * np.log(y_hat))
    return np.sum(np.log(y_hat[y==1])) + np.sum(np.log(1 - y_hat[y==0]))

def log_lik_cont(y_hat):
    return np.sum(np.log(y_hat))

def margconds_from_intersection(prob_ab, prob_a):
    """
    :param prob_ab: joint probability of A and B
    :param prob_a: marginal probability of A
    :return: 
        a_given_b: conditional probability of A given B
        prob_b: marginal probability of B
    """

    prob_a = prob_a.reshape(-1, 1)  # make it 2-D, for element-wise and matrix mults below
    b_given_a = (prob_ab / np.sum(prob_ab)) / prob_a
    prob_b = prob_a.T @ b_given_a
    prob_b[prob_b==0] = np.finfo(np.float64).tiny
    a_given_b = b_given_a * prob_a / prob_b

    # assert np.allclose(prob_b, 1.0) and np.allclose(np.sum(a_given_b, axis=0), [1.0, 1.0])
    # assert np.allclose(prob_a, 1.0) and np.allclose(prob_ab, 1.0)
    
    return a_given_b, prob_b.flatten()
    
    # TODO attempt to do for multiple headings at once...work in progress
    # proba_full = np.expand_dims(prob_a, axis=1)
    # b_given_a = (prob_ab / np.sum(prob_ab, axis=(0, 1))) / proba_full
    # prob_b = np.zeros_like(prob_a)
    # for d in range(prob_b.shape[1]):
    #     prob_b[:, d] = prob_a[:, d].T @ b_given_a[:, :, d]
    
    # a_given_b = b_given_a * proba_full / np.expand_dims(prob_b, axis=0)
    # return a_given_b, prob_b


def intersection_from_margconds(a_given_b, prob_a, prob_b):
    """
    Recover intersection of a and b using conditionals and marginals, according to Bayes theorem
    Essentially the inverse of margconds_from_intersection, and the two can be used together e.g.
    to update the intersections after adding a base_rate to prob_b

    :param a_given_b:
    :param prob_a:
    :param prob_b:
    :return:
        prob_ab: joint probability of A and B
        b_given_a
    """

    prob_ab = a_given_b * prob_b
    b_given_a = prob_ab / prob_a

    return prob_ab, b_given_a

# %% dataframe cleaning functions

def drop_breakfix(df, columns: Union[list, str]="choice") -> pd.DataFrame:
    """drop trials with breakfixes"""
    return df.dropna(subset=columns, axis=0)


def drop_one_target_wager(df, columns: Union[list, str]="oneTargChoice") -> pd.DataFrame:
    """drop one-target wager trials"""
    if isinstance(columns, str):
        columns = [columns]
    return df.loc[(df[columns] == 0).all(axis=1), :]


def drop_columns(df, columns) -> pd.DataFrame:
    """drop specified columns from dataframe"""
    df = df.loc[:, ~df.columns.str.startswith('Unnamed')] # drop any unnamed columns first
    return df.drop(columns, axis=1)


def zero_one_choice(df) -> pd.DataFrame:
    """recode choice column to be 0 and 1 instead of 1 and 2"""
    if np.max(df['choice']) == 2:
        df['choice'] -= 1
    return df


def drop_outlierRTs_bygroup(func):
    def wrapper(df, grp_cols=None, *args, **kwargs):
        return df.groupby(grp_cols).apply(func, *args, **kwargs).reset_index(drop=True)
    return wrapper


#@drop_outlierRTs_bygroup
def drop_outlierRTs(
    df,
    rt_range: Union[float, tuple] = (0, np.inf),
    metric: Literal["precomputed", "stdev", "percentile"] = "precomputed"
    ) -> pd.DataFrame:
    """drop trials with outlier RTs, based on specified metric

    Args:
        df (pd.DataFrame): behavioral data
        rt_range (tuple, optional): range of RTs to keep. Defaults to (0, np.inf).
        metric: str, optional): method to determine RT range. Options are:
            "precomputed": use rt_range as is (min, max)
            "stdev": use mean +/- rt_range * std
            "percentile": use rt_range as percentiles (min, max)
    """

    if metric == "stdev":
        assert isinstance(rt_range, (int, float)), "rt_range must be a single float for 'stdev' metric"
    else:
        assert isinstance(rt_range, tuple) and len(rt_range) == 2, "rt_range must be a tuple (min, max)"

    if metric == "precomputed":
        min_rt, max_rt = rt_range

    elif metric == "stdev":
        rt_std = df['RT'].std()
        rt_mu = df['RT'].mean()
        min_rt = np.max(rt_mu - rt_range*rt_std, 0)
        max_rt = rt_mu + rt_range*rt_std

    elif metric == "percentile":
        assert 0 <= rt_range[0] < rt_range[1] <= 100, "percentile range must be between 0 and 100"
        prc_rt = np.percentile(df['RT'], rt_range)
        min_rt, max_rt = prc_rt[0], prc_rt[1]

    return df.loc[(df['RT'] > min_rt) & (df['RT'] <= max_rt), :]


def bin_conditions(
    df,
    bin_ranges,
    bin_labels: Optional[dict] = None
    ) -> pd.DataFrame:
    """
    group multiple condition values into bins

    this is useful for consolidating the number of unique stimulus conditions 
    e.g. 

    bins = {
        'coherence': [0, 0.5, 1],
        'heading': [-14, -8, -4, -2, -1, 1, 2, 4, 8, 14],
        'delta': [-5, -2, 2, 5],
    }
    labels = {
        'coherence': [0.2, 0.7],
        'heading': [-12, -6, -3, -1.5, 0, 1.5, 3, 6, 12],
        'delta': [-3, 0, 3],
    }

    will group coherence values into low (0.2) and high (0.7), heading values into 9 bins,
    and delta values into 3 bins.
    Any coherence between 0 and 0.5 will be labeled as 0.2, and between 0.5 and 1 as 0.7, etc.
    Similarly any heading between -14 and -8 will be labeled as -12, etc.
    
    Args:
        df (pd.DataFrame): behavioral data
        bin_ranges (dict): dictionary specifying bin edges for each condition
        bin_labels (dict): dictionary specifying labels for each bin

    if bin_labels is None, the midpoint of each bin will be used as the label.
    Returns:
        pd.DataFrame: dataframe with binned condition columns
    """
    
    if bin_labels is None:
        bin_labels = {}
        for col, bins in bin_ranges.items():
            labels = []
            for i in range(len(bins)-1):
                labels.append((bins[i] + bins[i+1]) / 2)
            bin_labels[col] = labels
            
    for col, bins in bin_ranges.items():
        df[col] = pd.cut(df[col], bins=bins, labels=bin_labels[col])
        df[col] = df[col].astype('float')
        
    return df


def data_cleanup(
    filepath: Union[Path, str],
    drop_cols=None,
    save_file: bool = False
    ) -> pd.DataFrame:

    # TODO add kwargs for drop and binning parameters below, currently hardcoded...
    # TODO add print statements to explain, allow user-inputs to specify what functions to use?

    folder = "/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/dataStructs/"
   
    clean_filename = Path(filepath).stem + "_clean.csv"

    bhv_df = pd.read_csv(filepath)

    bins = {
        'coherence': [0, 0.5, 1],
        'heading': [-14, -8, -4, -2, -1, 1, 2, 4, 8, 14],
        'delta': [-5, -2, 2, 5],
    }
    labels = {
        'coherence': [0.2, 0.7],
        'heading': [-12, -6, -3, -1.5, 0, 1.5, 3, 6, 12],
        'delta': [-3, 0, 3],
    }

    if drop_cols is None:
        drop_cols = ["TargMissed", "confRT", "insertTrial", "filename", "subj", "trialNum",
                        "amountRewardLowConfOffered", "amountRewardHighConfOffered", "reward"]

    bhv_df_clean = (bhv_df
                    .pipe(drop_breakfix, columns=['choice', 'PDW'])
                    .pipe(drop_one_target_wager)
                    .pipe(zero_one_choice)
                    .pipe(drop_columns, columns=drop_cols)
                    .pipe(drop_outlierRTs, rt_range=(0.25, 2))
                    .pipe(bin_conditions, bin_ranges=bins, bin_labels=labels)
                    )

    # force all ves to be low coherence, by convention
    bhv_df_clean.loc[bhv_df_clean['modality'] == 1, 'coherence'] = bhv_df_clean['coherence'].min()

    # convert modality column to category type
    bhv_df_clean['modality'] = bhv_df_clean['modality'].astype('category')

    if save_file:
        bhv_df_clean.to_csv(PurePath(folder, clean_filename))
    
    return bhv_df_clean


def recode_onetarget_wagers(
    df: pd.DataFrame,
    remove_one_targ: bool = True
    ) -> pd.DataFrame:
    """Recode or remove one-target wager trials"""

    if remove_one_targ:
        # create new df with one-target PDW trials removed
        df = df.loc[df['oneTargConf'] == 0, :]

    else:
        # create a single category type column for PDW, with oneTarg trials coded as "2"
        df["PDW_1targ"] = df['PDW']
        df.loc[df['oneTargConf'] == 1, 'PDW_1targ'] = 2
        df["PDW_1targ"] = df['PDW_1targ'].astype("category")

    return df


def dots3DMP_create_trial_list(
    hdgs: list,
    mods: list,
    cohs: list,
    deltas: list,
    nreps: int = 1,
    shuff: bool = True
    ) -> pd.DataFrame:
    """Create a trial list of stimulus conditions for dots3DMP task"""

    if isinstance(shuff, int):
        np.random.seed(shuff) 

    num_hdg_groups = any([1 in mods]) + any([2 in mods]) * len(cohs) + \
        any([3 in mods]) * len(cohs) * len(deltas)
    hdg = np.tile(hdgs, num_hdg_groups)

    coh = np.empty_like(hdg, dtype=float)
    modality = np.empty_like(hdg, dtype=float)
    delta = np.empty_like(hdg, dtype=float)

    if 1 in mods:
        coh[:len(hdgs)] = cohs[0]
        delta[:len(hdgs)] = 0
        modality[:len(hdgs)] = 1
        last = len(hdgs)
    else:
        last = 0

    # visual has to loop over cohs
    if 2 in mods:
        for c in range(len(cohs)):
            these = slice(last + c*len(hdgs), last + c*len(hdgs) + len(hdgs))
            coh[these] = cohs[c]
            delta[these] = 0
            modality[these] = 2
        last = these.stop

    # combined has to loop over cohs and deltas
    if 3 in mods:
        for c in range(len(cohs)):
            for d in range(len(deltas)):
                here = last + c*len(hdgs)*len(deltas) + d*len(hdgs)
                these = slice(here, here + len(hdgs))
                coh[these] = cohs[c]
                delta[these] = deltas[d]
                modality[these] = 3

    # Now replicate times nreps and shuffle (or not):
    condlist = np.column_stack((modality, coh, delta, hdg))
    trial_table = np.tile(condlist, (nreps, 1))
    ntrials = len(trial_table)

    if shuff:
        if isinstance(shuff, int):
            trial_table = trial_table[np.random.default_rng(shuff).permutation(ntrials)]
        else:
            trial_table = trial_table[np.random.permutation(ntrials)]

    trial_table = pd.DataFrame(trial_table, columns=['modality', 'coherence', 'delta', 'heading'])

    return trial_table


def add_trial_outcomes(
    trial_table: pd.DataFrame,
    outcomes: Optional[dict] = None
    ) -> pd.DataFrame:
    """
    Replicate a trial table of stimulus conditions N times to include possible trial outcomes
    e.g. choice and wager

    Args:
        trial_table (pd.DataFrame): initial trial conditions list
        outcomes (dict, optional): dictionary specifying possible outcomes for each outcome type.
            Defaults to {'choice': [0, 1], 'PDW': [0, 1], 'oneTargConf': [0]}.

    Returns:
        pd.DataFrame: final trial conditions list, with additional columns for trial outcomes
    """
    if outcomes is None:
        outcomes = {'choice': [0, 1], 'PDW': [0, 1], 'oneTargConf': [0]}

    combinations = list(product(*outcomes.values()))
    df = pd.DataFrame(combinations, columns=outcomes.keys())
    
    new_trial_table = trial_table.loc[trial_table.index.repeat(len(df))].reset_index(drop=True)
    df_rep = pd.concat([df] * len(trial_table), ignore_index=True)
    new_trial_table = pd.concat([new_trial_table, df_rep], axis=1)

    return new_trial_table


def dots3DMP_create_conditions(
    conds_dict: dict[str, list],
    cond_labels: Optional[list[str]] = None
    ) -> pd.DataFrame:
    """
    Create a trial list of stimulus and response conditions for dots3DMP task
    """
    
    # TODO allow user to set stim and res_keys
    stim_keys = ['mods', 'cohs', 'deltas', 'hdgs']    
    res_keys = ['choice', 'correct', 'PDW', 'oneTargConf']
    
    ss_dict = {key: conds_dict.get(key, [0]) for key in stim_keys}
    conds = dots3DMP_create_trial_list(**ss_dict, nreps=1, shuff=False)
                
    rr_dict = {key: conds_dict.get(key, [0]) for key in res_keys}
    conds = conds.pipe(add_trial_outcomes, rr_dict)
    
    indices = [
        idx for idx, el in enumerate(stim_keys+res_keys) if el in conds_dict
               ]
    
    conds = conds[conds.columns[indices]]
    
    if cond_labels:
        conds.columns = cond_labels

    return conds
    
    
    
    
                               
    
    
        
                               
    
    