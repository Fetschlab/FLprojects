% dots3DMP_fitDDM_wrapper.m

% Generalized wrapper script for dots3DMP behavior fitting
% SJ 12-2022

%% try fitting simulated data to recover the generative parameters

cd /Users/stevenjerjian/Desktop/FetschLab/Analysis/data/dots3DMP_DDM
load tempsim_sepConfMaps_1000reps.mat

%% or real data

% load lucio_20220301-20221006_clean

%% set some vars

% set allNonHB == 0 to ignore unabsorbed probability, ie no need for weighted sum with max_dur
% in RT calculation, e.g. when sim excluded trials that failed to hit bound
if ~exist('allowNonHB','var'); allowNonHB=0; end

mods   = unique(data.modality);
cohs   = unique(data.coherence);
deltas = unique(data.delta);

options.RTtask   = 1; 
options.conftask = 2; % 1=continuous/rating, 2=PDW


%% select model + options

% ==== method ====
modelID=1; % unused for now
options.errfcn    = @dots3DMP_errfcn_DDM_2D_wConf_noMC; % 2D DDM aka anticorrelated race, for RT+conf [Kiani 14 / van den Berg 16 (uses Wolpert's images_dtb_2d (method of images, from Moreno-Bote 2010))]
options.fitMethod = 'global'; %'fms','global','multi','pattern','bads'
options.whichFit  = {'choice','RT'}; % choice, conf, RT, multinom (choice+conf)

% ==== implementation ====
options.dummyRun = 0; % dummyRun=1
options.ignoreUnabs = ~allowNonHB;
options.useVelAcc = 0; % use stimulus physical profiles in drift rates
% options.confModel = 'evidence+time'; % unused for now

% ==== diagnostics and output ====
options.plot      = 0;      % plot confidence maps
options.feedback  = 1;      % 0 - none, 1 - text feedback on LLs/err
options.runInterpFit = 0;   % model predictions for interpolated headings? for nice plotting


%% initialize parameters

guess = [origParams.kmult, origParams.B, origParams.theta, origParams.alpha, origParams.TndMean/1000];

guess = [50, 1.0, origParams.theta, origParams.alpha, 0.3, 0.3, 0.3];

fixed = zeros(1,length(guess));

% ************************************
% set all fixed to 1 for hand-tuning, or 0 for full fit
fixed(:)=1;
% ************************************

% or select some parameters to hold fixed
fixed = [0 0 1 1 1 1 0 0 0];


%% fit the model to data

f = parfeval(@dots3DMP_fitDDM,1,data,options,guess,fixed);

% X = dots3DMP_fitDDM(data,options,guess,fixed);

%% evaluate fitted parameters at actual headings (get overall fit error)

fixed = true(size(X)); % no actual fitting here

% use all stimulus conditions (i.e. including delta), evaluate fit on all outcomes
options.dummyRun = 0;
options.whichFit = {'choice','conf','RT'}; % choice, conf, RT, multinom (choice+conf)

options.feedback = 1;
[err_final, fit, parsedFit] = feval(options.errfcn,X,X,fixed,data,options);

%% plot the model predicted data points at this stage

dots3DMP_plots_fit_byCoh(data,fit,options.conftask,options.RTtask);
% dots3DMP_plots_fit_byConf(data,parsedFit,options.conftask,options.RTtask);

%% hold initial fit params fixed, refit for others
% e.g. we fit choice and RT, then come back for conf

fixed = [1 1 0 0 0 0 1 1 1]; % thetas and alpha are free params
options.whichFit = {'conf'}; % choice, conf, RT, multinom (choice+conf)

options.feedback = 2;
[err_final, fit, parsedFit] = feval(options.errfcn,X,X,fixed,data,options);

if options.runInterpFit

    fixed = true(size(X)); % again, no actual fitting

    % generate smooth fit curve with simulated linspaced headings
    hdgs = linspace(min(data.heading),max(data.heading),33); % odd so that there's a zero

    % Dfit is a dummy dataset with the same proportions of all trial types,
    [Dfit.heading, Dfit.modality, Dfit.coherence, Dfit.delta, ~] = dots3DMP_create_trial_list(hdgs,mods,cohs,deltas,1,0);

    % and the observables are just placeholders because we don't care about error
    % here, we just want the model predictions from fixed parameters
    Dfit.choice  = ones(size(Dfit.heading));
    Dfit.RT      = ones(size(Dfit.heading));
    Dfit.conf    = ones(size(Dfit.heading));
    Dfit.PDW     = ones(size(Dfit.heading));
    Dfit.correct = ones(size(Dfit.heading));

    % now calculate the interpolated fit
    [~, fitInterp] = feval(options.errfcn,X,X,fixed,Dfit,options);
else
    fitInterp = fit; % placeholder
end

%% plot the model predicted data points again now

dots3DMP_plots_fit_byCoh(data,fitInterp,options.conftask,options.RTtask);
% dots3DMP_plots_fit_byConf(data,parsedFit,options.conftask,options.RTtask);



