% dots3DMP_fitDDM_wrapper.m

% Generalized wrapper script for dots3DMP behavior fitting
% SJ 12-2022

% TO DO

% 


clear; close all;
addpath(genpath('/Users/stevenjerjian/Desktop/FetschLab/Analysis/codes/FLprojects/'))
% addpath(genpath('/Users/stevenjerjian/Desktop/FetschLab/Analysis/codes/third-party-codes/WolpertMOI_collapse'))
addpath(genpath('/Users/stevenjerjian/Desktop/FetschLab/Analysis/codes/third-party-codes/WolpertMOI'))

%% try fitting simulated data to recover the generative parameters

cd /Users/stevenjerjian/Desktop/FetschLab/Analysis/data/dots3DMP_DDM
load tempsim.mat

%% or real data

% load lucio_20220301-20221006_clean

%% some bookkeeping, then parse data, and plot if desired

if ~exist('allowNonHB','var'); allowNonHB=0; end

% ignore unabsorbed probability, ie no need for weighted sum with max_dur
% in RT calculation, e.g. when sim excluded trials that failed to hit bound
options.ignoreUnabs = ~allowNonHB;

options.RTtask   = 1; 
options.conftask = 2; % 1=continuous/rating, 2=PDW

% parse trial data into aggregated and other support vars
mods   = unique(data.modality);
cohs   = unique(data.coherence);
deltas = unique(data.delta);
hdgs   = unique(data.heading);

RTCorrOnly = 0;
if ~exist('parsedData','var')  % e.g., if simulation was run
    parsedData = dots3DMP_parseData(data,mods,cohs,deltas,hdgs,options.conftask,options.RTtask);
end

% **** 
% optional [data will be plotted below regardless, along with the fits]
% forTalk = 0;
% dots3DMP_plots(parsedData,mods,cohs,deltas,hdgs,options.conftask,options.RTtask)

% convert choice (back) to 0...1, for fitting
if max(data.choice(~isnan(data.choice)))==2
    data.choice = logical(data.choice-1);    
end


%% now the fitting itself

%****** first select which model to fit ********
modelID=1; options.errfcn = @dots3DMP_errfcn_DDM_2D_wConf_noMC; % 2D DDM aka anticorrelated race, for RT+conf [Kiani 14 / van den Berg 16 (uses Wolpert's images_dtb_2d (method of images, from Moreno-Bote 2010))]
%***********************************************

% options.nreps  = 100;
% options.confModel = 'evidence+time';

% SJ 10/2021, no longer doing model fits via Monte Carlo
options.runInterpFit = 0; 

options.fitMethod = 'fms'; %'fms','global','multi','pattern','bads'

guess = [origParams.kmult, origParams.B, origParams.theta, origParams.alpha, origParams.TndMean/1000];
% guess = guess.*(0.9+0.2.*rand(size(guess)));
% guess = [20, 1.5, origParams.theta, origParams.alpha, origParams.TndMean/1000];

fixed = zeros(1,length(guess));

% ************************************
% set all fixed to 1 for hand-tuning, or 0 for full fit
fixed(:)=1;
% ************************************

% fixed = [0 0 0 0 0 1 1 1 1];


options.plot      = 0; % plot confidence maps
options.feedback  = 2; % plot error trajectory (prob doesn't work with parallel fit methods)
options.useVelAcc = 0;
options.whichFit  = {'choice','RT'}; % choice, conf, RT, multinom (choice+conf)

%%
profile off;
[X, err_final, fit, parsedFit, fitInterp] = dots3DMP_fitDDM(data,options,guess,fixed);
% fitInterp is in fact not obsolete, and needs fixing in ^^

%% plot it!
dots3DMP_plots_fit_byCoh(data,fitInterp,options.conftask,options.RTtask);

dots3DMP_plots_fit_byConf(data,parsedFit,options.conftask,options.RTtask);

%% in progress

% check for fit-then-predict (missing values for +/- delta, etc)

% if any(unique(fit.delta(isnan(fit.choice)))==0) %--> should be nonzeros only
%     keyboard
% end
% this doesn't help, need to do the predict step in the fitDDM func

