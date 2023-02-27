import numpy as np
import pandas as pd


def prop_se(x):
    return np.sqrt((np.mean(x)*(1-np.mean(x))) / len(x))


def cont_se(x):
    return np.std(x) / np.sqrt(len(x))


def behavior_means(df, conftask=2, RTtask=True, by_delta=False, splitPDW=False):

    grp_list = ['modality', 'coherence', 'heading']

    if by_delta is True:
        grp_list.insert(-1, 'delta')

    if splitPDW is True and conftask > 0:

        if conftask == 1:
            df['PDW'] = df['conf'] > 0.5

        grp_list.append('PDW')

    # clean up dataframe columns and types
    df = df.loc[:, ~df.columns.str.startswith('Unnamed')]
    df['datetime'] = pd.to_datetime(
        df.loc[:, 'datetime'], infer_datetime_format=True)

    # leave choice, PDW and correct as int for convenience...
    df[['oneTargChoice', 'oneTargConf']] = df.loc[:, [
        'oneTargChoice', 'oneTargConf']].astype(bool)

    df['modality'] = df.loc[:, 'modality'].astype('category')

    if np.max(df['choice']) == 2:
        # to make 0..1 again, like PDW - more convenient for mean calculations
        df['choice'] -= 1

    # remove one-target choice
    df = df.loc[df['oneTargChoice'] == False]

    # remove cue conflict trials if by_delta is False
    if by_delta is False:
        df = df.loc[df['delta'] == 0]

    # for pHigh calculations, remove one-target confidence trials
    df_noOneTarg = df.loc[df['oneTargConf'] == False]

    if conftask == 1:
        pHigh = df_noOneTarg.groupby(by=grp_list)['conf'].agg(
            [np.mean, cont_se]).dropna(axis=0).reset_index()
    elif conftask == 2:
        pHigh = df_noOneTarg.groupby(by=grp_list)['PDW'].agg(
            [np.mean, prop_se]).dropna(axis=0).reset_index()
    else:
        pHigh = None

    pRight = df.groupby(by=grp_list)['choice'].agg(
        [np.mean, prop_se]).dropna(axis=0).reset_index()

    # this effectively achieves the same thing
    # pRight = pd.pivot_table(df, values='choice', index=['modality','coherence'], columns='heading',
    #     aggfunc=(np.mean,prop_se)).stack().reset_index()

    if RTtask is True:
        meanRT = df.groupby(by=grp_list)['RT'].agg(
            [np.mean, cont_se]).dropna(axis=0).reset_index()
    else:
        meanRT = None

    return pRight, pHigh, meanRT


# descriptive fit (none, logistic, or gaussian)
def behavior_fit(df, fitType='logistic', numhdgs=200):

    import statsmodels.api as sm
    #from sklearn.linear_model import LogisticRegression

    if np.max(df['choice']) == 2:
        df['choice'] -= 1

    hdgs = np.unique(df['heading']).reshape(-1, 1)
    xhdgs = np.linspace(np.min(hdgs), np.max(hdgs), numhdgs).reshape(-1, 1)

    outlist  = ['pRight','pHigh','meanRT']
    modnames = ['ves', 'vis', 'comb']
    attr_dict = {'intercept_': [], 'coef_': []}
    fit_results = dict(zip(outlist, [dict(zip(modnames, [attr_dict for m in modnames])) for outcome in outlist]))

    for modality in pRight['modality'].unique():
        for c, coh in enumerate(df['coherence'].unique()):
            
            if modality == 1:
                X = df.loc[df['modality'] == modality,
                                   ['heading', 'choice','PDW','RT']]
            else:
                X = df.loc[(df['modality'] == modality) & (df['coherence'] == coh),
                                    ['heading', 'choice','PDW','RT']]


            if fitType == 'logistic':

                pHigh_fit, meanRT_fit = None, None

                # using scikit-learn package
                #logreg = LogisticRegression().fit(
                #               choice_mc['heading'].to_numpy().reshape(-1, 1),
                #               choice_mc['choice'].to_numpy())
                #pRight_fit = logreg2.predict_proba(xhdgs)[:, 1]
                #coef_, intercept_ = logreg.coef_, logreg.intercept_

                logreg = sm.Logit(X['choice'], sm.add_constant(X['heading'])).fit()
                pRight_fit = logreg.predict(sm.add_constant(xhdgs))
                coef_, intercept_ = logreg.params['heading'], logreg.params['const']

                # store logreg attributes
                fit_results['pRight'][modnames[modality-1]
                                   ]['intercept_'].append(intercept_)
                fit_results['pRight'][modnames[modality-1]
                                   ]['coef_'].append(coef_)

            elif fitType == 'gaussian':
                
                probreg = sm.Probit(X['choice'], sm.add_constant(X['heading'])).fit()
                pRight_fit = probreg.predict(sm.add_constant(xhdgs))
                coef_, intercept_ = probreg.params['heading'], probreg.params['const']

                # store logreg attributes
                fit_results['pRight'][modnames[modality-1]
                                   ]['intercept_'].append(intercept_)
                fit_results['pRight'][modnames[modality-1]
                                   ]['coef_'].append(coef_)
                
                probreg = sm.Probit(1-X['PDW'], sm.add_constant(X['heading'])).fit()
                pRight_fit = probreg.predict(sm.add_constant(xhdgs))
                
                # store logreg attributes
                fit_results['pRight'][modnames[modality-1]
                                   ]['intercept_'].append(intercept_)
                fit_results['pRight'][modnames[modality-1]
                                   ]['coef_'].append(coef_)
                
                


    return xhdgs, pRight_fit, pHigh_fit, meanRT_fit


def plot_behavior_means(pRight, pHigh, meanRT):

    import matplotlib.pyplot as plt

    condcols = ['k', 'r', 'b']
    fig, ax = plt.subplots(3, 2, figsize=(8, 12))

    def plot_fcn(data, mods, cohs): # to apply a function instead of loop
        pass

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


if __name__ == "__main__":

    import glob
    import pandas as pd

    subject = 'lucio'
    folder = '/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/dataStructs/'
    files = glob.glob(f'{folder}*.csv')
    df = pd.read_csv(files[0])

    pRight, pHigh, meanRT = behavior_means(
        df, conftask=2, RTtask=True, by_delta=False, splitPDW=False)
    plot_behavior_means(pRight, pHigh, meanRT)
