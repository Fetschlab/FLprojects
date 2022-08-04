%% dots3DMP behavioral plots and analyses
% SJ 10/2021

% for SfN Poster

% datasets
% lucio - RT, Mar to Aug 2021
% human - non-RT up to Feb 2020
% human - RT, Mar 2020 to Oct 2021
% simulated data, different models

clear; clc; close all
cd /Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/dataStructs
addpath(genpath('/Users/stevenjerjian/Desktop/FetschLab/Analysis/codes/'))

%% select subject, load the data

subject = 'lucio';

switch subject
    
    case 'lucio'
%         load('lucio_20210315-20210805_clean.mat') % recent lucio data, PDW + RT
        load('lucio_20211101-20220602_clean.mat') % recent lucio data, PDW + RT

        conftask = 2; % 1=colorbars, 2=PDW
        RTtask   = 1;
        
        RTlims = [0.25 2.25];

    case 'human'

        conftask = 1;
        RTtask   = 1; % change this to select RT or non-RT dataset
       
        if ~RTtask
            load('human_20190625-20191231_nonRT_clean.mat');% human non-RT, SfN 2021
            
        else
%             load('human_20200213-20210922_RT_clean.mat') % human RT, SfN 2021
%             load('human_20200213-20220113_RT_clean_Jan2022.mat') % human RT, Jan 2022
            load('human_20200213-20220317_RT_clean_Mar2022.mat') % human RT, Jan 2022

            RTlims = [0.25 2.5];
        end
        
        
    case 'simul' % load simulated data

        load('2DAccSim_conftask2_8100trs.mat')
        % conftask & RTtask should already be saved in file
end

fnames = fieldnames(data);

if RTtask
    removethese = data.RT < RTlims(1) | data.RT > RTlims(2);
    
    for f=1:length(fnames)
        data.(fnames{f})(removethese) = [];
    end
end

if strcmp(subject,'human') && RTtask
    % not enough good data, so let's just remove for now?
    removethese = data.heading==0;
    for f=1:length(fnames)
        data.(fnames{f})(removethese) = [];
    end
end

mods   = unique(data.modality); 
cohs   = unique(data.coherence); 
deltas = unique(data.delta);
% deltas = [-3 3];
hdgs   = unique(data.heading);

%% basic parsing and descriptive gaussian fits 

% means per condition, logistic fits
parsedData = dots3DMP_parseData(data,mods,cohs,deltas,hdgs,conftask,RTtask); 

% gaussian fits
gfit = dots3DMP_fit_cgauss(data,mods,cohs,deltas,conftask,RTtask); 

%% summary plots
% logistic fit plots
% dots3DMP_plots(parsedData,mods,cohs,deltas,hdgs,conftask,RTtask)

% separate subplots for each coh, with all mods on same subplot
dots3DMP_plots_cgauss_byCoh(gfit,parsedData,mods,cohs,deltas,hdgs,conftask,RTtask)

% or...separate subplots for each mod/delta, and all cohs on same subplot
% n.b. - this one needs tidying to look nice
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

% newer version, SJJ late October 2021
% plotOption == 0 - plot errors/low bet only
% plotOption == 1 - plot correct/high bet only
% plotOption == 2 - plot correct/error or high/low bet separately
% plotOption == -1 - plot all trials
dots3DMP_RTquantiles(data,conftask,1)

%% PDW and RT for correct vs incorrect trials

dots3DMP_CorrectVsErrorCurves(data,conftask,RTtask,1)

%% Psychometric curves split by PDW

% Choice curves for low bet should be shallower - indicative that PDW is
% predictive of accuracy, even within a stimulus condition!

parsedData_byConf = dots3DMP_parseData_byConf(data,mods,cohs,deltas,hdgs,conftask,RTtask); 
% gfit_byConf       = dots3DMP_fit_cgauss_byConf(data,mods,cohs,deltas,conftask,RTtask);

% plot it
dots3DMP_plots_byConf(parsedData_byConf,mods,cohs,deltas,hdgs,conftask,RTtask);
% dots3DMP_plots_cgauss_byConf(gfit_byConf,parsedData_byConf,mods,cohs,deltas,hdgs,conftask,RTtask)

%% Psychometric curves split by (previous) PDW

% doesn't seem like nBack 1 has any influence on accuracy/RT on subsequent
% trial (within or across mod)...
% something weird with PDW though, need to check code logic
nBack = 1; withinMod = 0;
parsedData = dots3DMP_parseDataPrevConf(data,mods,cohs,deltas,hdgs,conftask,RTtask,nBack,withinMod);
dots3DMP_plots_byConf(parsedData,mods,cohs,deltas,hdgs,conftask,RTtask);

%% split curves for PDW, split  by relative reward for high bet

rewRatio = data.amountRewardHighConfOffered ./ data.amountRewardLowConfOffered;
nbins = 4;
confQ = [0 quantile(rewRatio,nbins-1) inf];
confGroup = discretize(rewRatio, confQ); 
splitPDW  = 2; isPDW = 0;
% splitPDW = 0 - don't split trials at all (show all, together), 1 - show only low bets, 2 -
% show only high bets, 3 - show only 1-target, 4 - show all, separately

% just plot behavioral outcomes split by PDW, and
% make 1-target trials a separate category (i.e. line)
% can remove them totally by setting removeOneTarg=0;
% confGroup = double(data.PDW)+1;
% confGroup(logical(data.oneTargConf))= 3;
% splitPDW = 0; % this is redundant in this case
% isPDW = 1;

removeOneTarg = 0;
parsedData = dots3DMP_parseData_multiConf(data,mods,cohs,deltas,hdgs,confGroup,conftask,RTtask,removeOneTarg,splitPDW); % don't remove 1-targets, and don't split by hi/lo, so we can plot P(high bet) as function of reward ratio
dots3DMP_plots_multiConf(parsedData,mods,cohs,deltas,hdgs,conftask,RTtask,splitPDW,isPDW)


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


