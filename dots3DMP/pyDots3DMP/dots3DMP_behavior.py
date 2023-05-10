# %%

# TODO
# functions for:
# 1. extracting behavior from neural recording dataset
# 2. RT quantiles, plot vs confidence/accuracy
# 3. correct vs error metrics
# 4. decorator to loop any function over a grouping variable (i.e. day/block)
# 5. 

# regression analyses
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# import seaborn.objects as so
import glob

import statsmodels.api as sm
from scipy.optimize import curve_fit
# from sklearn.linear_model import LogisticRegression

# %% helper functions


def prop_se(x):
    return np.sqrt((np.mean(x)*(1-np.mean(x))) / len(x))


def cont_se(x):
    return np.std(x) / np.sqrt(len(x))


def gaus(x, ampl, mu, sigma, bsln):
    return ampl * np.exp(-(x-mu)**2 / (2*sigma**2)) + bsln


def clean_behavior_df(df):

    # clean up dataframe columns and types
    df = df.loc[:, ~df.columns.str.startswith('Unnamed')]

    # convert datetime to actual pandas datetime
    df['datetime'] = pd.to_datetime(
        df.loc[:, 'datetime'], infer_datetime_format=True)

    df['modality'] = df.loc[:, 'modality'].astype('category')

    if np.max(df['choice']) == 2:
        # to make 0..1 again, like PDW - more convenient for calculations
        df['choice'] -= 1

    # remove one-target choice
    df = df.loc[df['oneTargChoice'] == 0]

    # note we leave oneTargConf in, as this might be useful for some analyses
    # and can be removed easily later if desired
    
    # remove outlier RTs
    # hard limits, or based on things beyond percentile/SDs from mean range
    

    return df

# %%


def behavior_means(df, conftask=2, RTtask=True,
                   by_delta=False, splitPDW=False):
    """

    Parameters
    ----------
    df : TYPE
        DESCRIPTION.
    conftask : TYPE, optional
        DESCRIPTION. The default is 2.
    RTtask : TYPE, optional
        DESCRIPTION. The default is True.
    by_delta : TYPE, optional
        DESCRIPTION. The default is False.
    splitPDW : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    pRight : TYPE
        DESCRIPTION.
    pHigh : TYPE
        DESCRIPTION.
    meanRT : TYPE
        DESCRIPTION.

    """

    grp_list = ['modality', 'coherence', 'heading']

    if by_delta is True:
        grp_list.insert(-1, 'delta')

    if splitPDW is True and conftask > 0:

        if conftask == 1:
            df['PDW'] = df['conf'] > 0.5

        grp_list.append('PDW')

    # remove cue conflict trials if by_delta is False
    if by_delta is False:
        df = df.loc[df['delta'] == 0]

    # for pHigh calculations, remove one-target confidence trials
    df_noOneTarg = df.loc[df['oneTargConf'] == 0]

    # pHigh calculation
    pHigh = None
    if conftask == 1:
        pHigh = df_noOneTarg.groupby(by=grp_list)['conf'].agg(
            [np.mean, cont_se]).dropna(axis=0).reset_index()
    elif conftask == 2:
        pHigh = df_noOneTarg.groupby(by=grp_list)['PDW'].agg(
            [np.mean, prop_se]).dropna(axis=0).reset_index()

    # pRight calculation
    pRight = df.groupby(by=grp_list)['choice'].agg(
        [np.mean, prop_se]).dropna(axis=0).reset_index()

    # this effectively achieves the same thing
    # pRight = pd.pivot_table(df, values='choice',
    #       index=['modality','coherence'], columns='heading',
    #       aggfunc=(np.mean,prop_se)).stack().reset_index()

    # meanRT calculation
    meanRT = None
    if RTtask is True:
        meanRT = df.groupby(by=grp_list)['RT'].agg(
            [np.mean, cont_se]).dropna(axis=0).reset_index()

    return pRight, pHigh, meanRT

# %% descriptive fit to basic curves (none, logistic, or gaussian)


def behavior_fit(df, fitType='logistic', by=None, numhdgs=200):

    # TODO implement by parameter
    # if by == 'correct', fit PDW and RT separately for correct/error trials
    # if by == 'PDW', fit choice and RT separatley for different wagers
    #                                        (high, low, one-targ)

    if np.max(df['choice']) == 2:
        df['choice'] -= 1

    hdgs = np.unique(df['heading']).reshape(-1, 1)
    xhdgs = np.linspace(np.min(hdgs), np.max(hdgs), numhdgs).reshape(-1, 1)

    outlist = ['pRight', 'PDW', 'RT']
    modnames = ['ves', 'vis', 'comb']
    attr_dict = {'intercept_': [], 'coef_': [], 'prediction': []}
    fit_results = dict(zip(outlist, [dict(zip(modnames,
                                              [attr_dict for m in modnames]))
                                     for outcome in outlist]))


    # always do fits separately for each mod/coh
    for modality in pRight['modality'].unique():
        for c, coh in enumerate(df['coherence'].unique()):

            if modality == 1:
                X = df.loc[df['modality'] == modality,
                           ['heading', 'choice', 'PDW', 'RT']]
            else:
                X = df.loc[(df['modality'] == modality) &
                           (df['coherence'] == coh),
                           ['heading', 'choice', 'PDW', 'RT']]

            if fitType == 'logistic':

                # using scikit-learn package
                # logreg = LogisticRegression().fit(
                #               X['heading'].to_numpy().reshape(-1, 1),
                #               X['choice'].to_numpy())
                # yhat = logreg.predict_proba(xhdgs)[:, 1]
                # coef_, intercept_ = logreg.coef_, logreg.intercept_

                logreg = sm.Logit(X['choice'], sm.add_constant(X['heading'])).fit()
                yhat = logreg.predict(sm.add_constant(xhdgs))
                coef_, intercept_ = logreg.params['heading'], logreg.params['const']

            elif fitType == 'gaussian':

                probreg = sm.Probit(X['choice'], sm.add_constant(X['heading'])).fit()
                yhat = probreg.predict(sm.add_constant(xhdgs))
                coef_, intercept_ = probreg.params['heading'], probreg.params['const']

                # fits for PDW and RT
                popt_PDW, _ = curve_fit(gaus, xdata=X['heading'],
                                        ydata=1-X['PDW'], p0=[0.1, 0, 3, 0.5])

                popt_RT, _ = curve_fit(gaus, xdata=X['heading'],
                                       ydata=X['RT'], p0=[0.1, 0, 3, 0.5])

                fitPDW = 1 - gaus(xhdgs, *popt_PDW)
                fitRT = gaus(xhdgs, *popt_RT)

                fit_results['PDW'][modnames[modality-1]]['prediction'].append(fitPDW)
                fit_results['PDW'][modnames[modality-1]]['coef_'].append(popt_PDW)

                fit_results['RT'][modnames[modality-1]]['prediction'].append(fitRT)
                fit_results['RT'][modnames[modality-1]]['coef_'].append(popt_RT)

            fit_results['pRight'][modnames[modality-1]]['prediction'].append(yhat)
            fit_results['pRight'][modnames[modality-1]]['intercept_'].append(intercept_)
            fit_results['pRight'][modnames[modality-1]]['coef_'].append(coef_)

    return xhdgs, fit_results

# %%


def cue_weighting(fit_results):

    wves_emp = wves_pred = None

    for res in fit_results.keys():
        ...
        # TODO calculate wves pred and wves emp for each of pRight, PDW, RT
        # TODO first fix the behavior_fit function to do by_delta as well...
    return wves_emp, wves_pred

# %%


def plot_behavior_means(pRight, pHigh, meanRT):
    """


    Parameters
    ----------
    pRight : TYPE
        DESCRIPTION.
    pHigh : TYPE
        DESCRIPTION.
    meanRT : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """

    # (
    #     so.Plot(pRight, x='heading',y='mean', color='modality')
    #     .facet('coherence')
    #     .add(so.Line())
    # )

    condcols = ['k', 'r', 'b']
    fig, ax = plt.subplots(3, 2, figsize=(8, 12))

    for modality in pRight['modality'].unique():

        for c, coh in enumerate(pRight['coherence'].unique()):

            if modality == 1:
                choice_mc = pRight[pRight['modality'] == modality]
                pdw_mc = pHigh[pHigh['modality'] == modality]
                rt_mc = meanRT[meanRT['modality'] == modality]
            else:
                choice_mc = pRight[(pRight['modality'] == modality) & (
                    pRight['coherence'] == coh)]
                pdw_mc = pHigh[(pHigh['modality'] == modality)
                               & (pHigh['coherence'] == coh)]
                rt_mc = meanRT[(meanRT['modality'] == modality)
                               & (meanRT['coherence'] == coh)]

            plt.sca(ax[0][c])
            ax[0][c].errorbar(choice_mc['heading'], choice_mc['mean'],
                              choice_mc['prop_se'], ecolor=condcols[modality-1])
            ax[0][c].set_xticks(choice_mc['heading'])

            plt.sca(ax[0][c])
            ax[1][c].errorbar(choice_mc['heading'],
                              pdw_mc['mean'], pdw_mc['prop_se'])
            ax[1][c].set_xticks(choice_mc['heading'])

            plt.sca(ax[0][c])
            ax[2][c].errorbar(choice_mc['heading'],
                              rt_mc['mean'], rt_mc['cont_se'])
            ax[2][c].set_xticks(choice_mc['heading'])
    plt.show()

# %%


def RTquantiles(df, nq=5, cond_columns=['modality', 'coherence', 'heading'],
                depvar='PDW'):

    df['RTq'] = df.groupby(cond_columns)['RT'].transform(
        lambda x: pd.qcut(x, nq, labels=False))

    # find a way to return the RT (mean/median) of each quantile
    quantile_vals = df.groupby(cond_columns)['RT'].transform(
        lambda x: pd.qcut(x, nq))

    result = df.groupby(cond_columns + ['RTq'])[depvar].mean().reset_index()

    return result




# %% define trial list


def dots3DMP_create_trial_list(hdgs, mods, cohs, deltas, nreps, shuff=True):

    # if shuff:
    #     np.random.seed(42)  # for reproducibility

    num_hdg_groups = any([1 in mods]) + any([2 in mods]) * len(cohs) + \
        any([3 in mods]) * len(cohs) * len(deltas)
    hdg = np.tile(hdgs, num_hdg_groups)

    coh = delta = modality = np.empty_like(hdg)

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
    condlist = np.column_stack((hdg, modality, coh, delta))
    trial_table = np.tile(condlist, (nreps, 1))
    ntrials = len(trial_table)

    if shuff:
        trial_table = trial_table[np.random.permutation(ntrials)]

    # TODO check this
    trial_table = np.stack(trial_table.T, axis=1)
    trial_table = pd.DataFrame(trial_table[:, [1, 2, 0, 3]],
                               columns=['modality', 'coherence',
                                        'heading', 'delta'])

    return trial_table, ntrials


# %% main

if __name__ == "__main__":

    subject = 'lucio'
    folder = '/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/dataStructs/'
    files = glob.glob(f'{folder}*.csv')
    df = pd.read_csv(files[0])
    # TODO probably just specify the file...
    
    df = clean_behavior_df(df)
    
    pRight, pHigh, meanRT = behavior_means(
        df, conftask=2, RTtask=True, by_delta=False, splitPDW=False)

    # plot_behavior_means(pRight, pHigh, meanRT)
    xhdgs, fit_results = behavior_fit(df, fitType='gaussian')

    # mods = np.array([1, 2, 3])
    # cohs = np.array([1, 2])
    # hdgs = np.array([-12, -6, -3, -1.5, 0, 1.5, 3, 6, 12])
    # deltas = np.array([-3, 0, 3])
    # nreps = 1

    # trial_table, ntrials = \
    #     dots3DMP_create_trial_list(hdgs, mods, cohs, deltas,
    #                                nreps, shuff=False)