%% dots3DMP behavioral plots and analyses
% SJ 10/2021

% for SfN Poster

% datasets
% lucio - RT, Mar to Aug 2021
% human - non-RT up to Feb 2020
% human - RT, Mar 2020 to Oct 2021
% simulated data, different models

clear; clc; close all

% check that plotting and labelling works correctly
% compare RT profiles in ves and comb split by heading and coherence

%% select subject, load the data

datapath = '/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/dataStructs';
subject = 'zarya';
export_figs = 0;

fig_folder = '/Users/stevenjerjian/Desktop/FetschLab/Analysis/figs';

conftask = 2;
RTtask   = 1;

data = dots3DMP_loadBehaviorData(subject,datapath,conftask,RTtask);
RTlims = [0.25 2.25];

if ~isfield(data,'oneTargConf')
    data.oneTargConf = false(size(data.heading));
end

fnames = fieldnames(data);


if strcmp(subject,'human')
    if RTtask
        % not enough good data, so let's just remove for now?
        removethese = data.heading==0;
        for f=1:length(fnames)
            data.(fnames{f})(removethese) = [];
        end
    end
    removethese = strcmp(data.subj,'SBG');
    for f=1:length(fnames)
        data.(fnames{f})(removethese) = [];
    end
end


if RTtask
    removethese = data.RT < RTlims(1) | data.RT > RTlims(2);

    fprintf('%d trials outside of RT lims...removing\n',sum(removethese));
    
    for f=1:length(fnames)
        data.(fnames{f})(removethese) = [];
    end
end

mods   = [1 2 3]; %unique(data.modality); 
cohs   = unique(data.coherence); 
deltas = unique(data.delta);
% deltas = [-3 3];
deltas = 0;
hdgs   = unique(data.heading);

%% Zarya, remove older data

removethese = data.date >= 20230901;
for f=1:length(fnames)
    data.(fnames{f})(removethese) = [];
end

%% basic parsing and plot of logistic fits

% means per condition, logistic fits
parsedData = dots3DMP_parseData(data,mods,cohs,deltas,hdgs,conftask,RTtask); 
dots3DMP_plots(parsedData,mods,cohs,deltas,hdgs,conftask,RTtask)

%% gaussian fits and plots
gfit = dots3DMP_fit_cgauss(data,mods,cohs,deltas,conftask,RTtask); 

% separate subplots for each coh, with all mods on same subplot
gfits_fig = dots3DMP_plots_cgauss_byCoh(gfit,parsedData,mods,cohs,deltas,hdgs,conftask,RTtask);

if export_figs
    figname = sprintf('%s_%s_gfit_behavior.pdf',subject,date);
%     savefig(fullfile(figfolder,figname));
    exportgraphics(gfits_fig,fullfile(figfolder,figname),'Resolution',300); 
end

% or separate subplots for each mod/delta, and all cohs on same subplot -
% needs work to look nice if it's going to be used publicly
% dots3DMP_plots_cgauss_byModDelta(gfit,parsedData,mods,cohs,deltas,hdgs,conftask,RTtask)

%% psychophysical cue weights
% assume wvis always = 1-wves

wves = dots3DMP_cueWeights(gfit,cohs,deltas,conftask);

% bootstrapping for error bars
rng(28); % for reproducibility
nboots = 100;
[gfitBoot,wvesBoot] = dots3DMP_cgauss_bootstrap_func(data,gfit,mods,cohs,deltas,nboots,conftask,RTtask);

% plot the weights
dots3DMP_plotCueWeights(wves,wvesBoot,cohs,conftask)

%% Confidence as function of decision time (RT quantiles)

if ~isfield(data,'correct')
    data.correct = (data.heading>0 & data.choice==2) | (data.heading<0 & data.choice==1) | (data.heading==0 & rand<0.5);
end

% third argument specifies which trials to use 
% -1: all, 0: error trials only, 1: correct trials only
% -1 will also lead to plotting of p(correct) as function of RT quantiles
% 0 will only use weaker stimuli (set to bottom 3), since stronger stimuli
% produce few errors 


% dots3DMP_RTquantiles(data,conftask,-1)

% newer version, SJJ late October 2021
% plotOption == 0 - plot errors/low bet only
% plotOption == 1 - plot correct/high bet only
% plotOption == 2 - plot correct/error or high/low bet separately
% plotOption == -1 - plot all trials
RTquant_fh = dots3DMP_RTquantiles(data,conftask,4,-1);

if export_figs
    exportgraphics(RTquant_fh(1),'RTquantiles.pdf','Resolution',300,'append',1);
    exportgraphics(RTquant_fh(2),'RTquantiles.pdf','Resolution',300,'append',1);
end


%% PDW and RT for correct vs incorrect trials

dots3DMP_CorrectVsErrorCurves(data,conftask,RTtask,1)

%% Psychometric curves split by PDW

% Choice curves for low bet should be shallower - indicative that PDW is
% predictive of accuracy, even within a stimulus condition!

cohs = 0.2;
parsedData_byConf = dots3DMP_parseData_byConf(data,mods,cohs,deltas,hdgs,conftask,RTtask); 
gfit_byConf       = dots3DMP_fit_cgauss_byConf(data,mods,cohs,deltas,conftask,RTtask);

% plot it
fh = dots3DMP_plots_cgauss_byConf(gfit_byConf,parsedData_byConf,mods,cohs,deltas,hdgs,conftask,RTtask);

%% split by reward offered
rewRatio = data.amountRewardHighConfOffered ./ data.amountRewardLowConfOffered;
nbins = 4;
confQ = [0 quantile(rewRatio,nbins-1) inf];
confGroup = discretize(rewRatio, confQ); 

splitPDW = 0;
removeOneTarg = 0;

clrseqs = {'Greys','Reds','Blues'};

parsedData_conf = dots3DMP_parseData_multiConf(data,mods,cohs,deltas,hdgs,confGroup,conftask,RTtask,removeOneTarg,splitPDW); % don't remove 1-targets, and don't split by hi/lo, so we can plot P(high bet) as function of reward ratio
nConfGroups = length(parsedData_conf.confGroups);

clear splitcols
for m = 1:length(mods)
    if splitPDW
        splitcols{mods(m)} = cbrewer('seq',clrseqs{m},ceil(nConfGroups/2));
    else
        splitcols{mods(m)} = cbrewer('seq',clrseqs{m},nConfGroups);
    end
end

% confQ_means = accumarray(confGroup,rewRatio,[nbins 1],@mean);

for b = 1:nbins
    parsedData_conf.confGroupLabels{b} = num2str(mean(rewRatio(confGroup==b)));
end
dots3DMP_plots_multiConf(parsedData_conf,mods,cohs,deltas,hdgs,conftask,RTtask,splitPDW,splitcols)


%% split high,low,1-targ (this should supersede above '_byConf')
% % this only works if the dataset still includes oneTargConf trials!
confGroup = double(data.PDW)+1;
confGroup(logical(data.oneTargConf))= 3;
splitPDW = 1;
removeOneTarg = 0;

parsedData_conf = dots3DMP_parseData_multiConf(data,mods,cohs,deltas,hdgs,confGroup,conftask,RTtask,removeOneTarg,splitPDW); % don't remove 1-targets, and don't split by hi/lo, so we can plot P(high bet) as function of reward ratio
nConfGroups = length(parsedData_conf.confGroups);

for m = 1:length(mods)
        splitcols{mods(m)} = cbrewer('qual','Dark2',nConfGroups);
end

parsedData_conf.confGroupLabels = {'Low','High','1-targ'};
dots3DMP_plots_multiConf(parsedData_conf,mods,cohs,deltas,hdgs,conftask,RTtask,splitPDW,splitcols)


%% relationship between confidence and cue weights

% calculate weights separately for low and high coh?
% CF notes - this probably doesn't make sense to try and do however - 
% we don't know mapping between high conf in single cues to high conf in comb, or if one even exists!

% anyway, it's easy to do technically, just use gfit_byConf from above
    
% wves_byConf = dots3DMP_cueWeights(gfit_byConf,cohs,deltas,conftask,1);

% or we can do a piecewise way (part 2 of this function)
% this function actually runs three separate analyses relating confidence
% and cue conflict

% 1. compare shifts of choice mu vs conf mu in non-conflict and conflict conditions 
% 2. piecewise comparison of PRight for fixed absolute heading, different conflicts as function of coh
% 3. average confidence in conflict vs no conflict (low headings only)

dots3DMP_ConfDelta(data,gfit,cohs,deltas,hdgs,conftask,RTtask,3)


%% MODELLING

% options.errfun = 'dots3DMP_fit_2Dacc_err_sepbounds_noSim';
options.errfun = 'dots3DMP_fit_2Dacc_err_noSim';
options.fitMethod = 'fms'; %'fms','fmsbnd','global','multi','pattern','bads'
options.whichFit = 3; % 0 - choice only, 1 - choice + RT, 2 - choice + conf, 3 - ALL
% options.paramNames = {'kves','kvisLo','kvisHi','BVes','BVis','BComb','TndVe','TndVi','TndCo','T-Conf','theta','cLVes','cLVis','cLComb'};
options.sepbounds = 0;

fixAll = 0;

% initial guess (or hand-tuned params)
% kmult   = 0.3;
% kvis    = kmult.*cohs';
% kves    = mean(kvis);
% kves    = 46;
% B       = 0.5;
% Tnd     = [0.4 0.5 0.4];
% Ttc     = 0; % time to confidence, ignored for now...

% hand-tuning 10/25/2021 196 runs
% kves = 0.237;
% kvis = [0.103 0.33];
% B    = [0.330 0.383671 0.342];
% Tnd  = [0.48 0.62 0.55];
% Ttc  = 0;
% cL   = [0.07 0.22 0.085]; % confLapse rate, lapse rate of high bets
% fixed = [0 0 0 0 0 0 0 0 0 1];
% options.paramNames = {'kves','kvisLo','kvisHi','BVes','BVis','BComb','TndVe','TndVi','TndCo','T-Conf','theta','cLVes','cLVis','cLComb'};

% one bound, lucio RT
% kves = 0.2213;
% kvis = [0.1504 0.33];
% B    = 0.33;
% Tnd  = [0.48 0.68 0.56];
% Ttc  = 0;
% cL   = [0.07 0.38 0.09]; % confLapse rate, lapse rate of high bets/high conf
% fixed = [0 0 0 0 0 0 0 1];
% options.paramNames = {'kves','kvisLo','kvisHi','B','TndVe','TndVi','TndCo','T-Conf','theta','cLVes','cLVis','cLComb'};

% one bound, human (RT)
% kves = 0.15;
% kvis = [0.05 0.2];
% B    = 1.2;
% Tnd  = [0.75 0.75 0.75];
% Ttc  = 0;
% cL   = [0.25 0.25 0.25]; % confLapse rate, lapse rate of high bets, or in human case, random

% lucio
kves = 0.23;
kvis = [0.15 0.32];

B = 0.33;
Tnd  = [0.49 0.69 0.56];
Ttc  = 0;
cL = [0.07 0.38 0.09];
theta = 0.07;



%%
% options.paramNames = {'kves','kvisLo','kvisHi','B','TndVe','TndVi','TndCo','T-Conf','cLVes','cLVis','cLComb'};
options.paramNames = {'kves','kvisLo','kvisHi','B','TndVe','TndVi','TndCo','T-Conf','theta','cLVes','cLVis','cLComb'};

guess   = [kves kvis B Tnd Ttc];

if RTtask
    fixed = [0 0 0 0 0 0 0 1]; % RT
else
    fixed = [0 0 0 0 1 1 1 1]; % non-RT
end


if conftask==1
    guess = [guess cL];
    fixed = [fixed 0 0 0];
elseif conftask==2 % PDW
    theta = 0.07;
    guess = [guess theta cL];
    fixed = [fixed 0 0 0 0];
end

options.lowerbound = [0.05 0.05 0.2 0.2 0.2 0.2 0.2 0 0.1 0.1 0.1];
options.upperbound = [2 2 2 2 1 1 1 1 0.3 0.3 0.3];


% ************************************
% set all fixed to 1 for hand-tuning, or 0 for full fit
if fixAll, fixed(:)=1; end
% ************************************

% plot error trajectory (prob doesn't work with parallel fit methods)
options.ploterr  = 1;
options.dummyRun = 0;
options.RTtask   = RTtask;
options.conftask = conftask; % 1 - sacc endpoint, 2 - PDW

if options.ploterr, options.fh = 400; end

options.runInterpFit = 1;
[X, err_final, fit, fitInterp] = dots3DMP_fitDDM(data,options,guess,fixed);
dots3DMP_plots_fit_byCoh(data,fitInterp,conftask,RTtask,0);


