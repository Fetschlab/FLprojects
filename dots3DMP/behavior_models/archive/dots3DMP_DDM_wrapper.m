% dots3DMP DDM wrapper


% WIP

%% MODEL SPECIFICATIONS and CONDITIONS

modelID  = 1; % 1 will be 2Dacc model ('Candidate model' against which others are tested).
% modelVar = 1; % variation within model ID, e.g. velocity acceleration coding, confidence model, cue weighting in combined..


% task type
RTtask   = 1;
conftask = 2; % 1 - sacc endpoint, 2 - PDW

nreps = 300; % number of repetitions of each unique trial type % (ntrials depends on num unique trial types)

% stimulus conditions
mods  = [1 2 3]; % stimulus modalities: ves, vis, comb
cohs  = [0.4 0.8]; % visual coherence levels (these are really just labels, since k's are set manually)
hdgs  = [-12 -6 -3 -1.5 0 1.5 3 6 12];
% deltas = [-3 0 3]; % conflict angle; positive means vis to the right
deltas  = 0;


% time information
dT      = 1; % time step, ms
max_dur = 2100; % stimulus duration (ms)

% maybe these become obsolete with modelVar eventually
% confModel = 'evidence+time'; % 'evidence+time','evidence_only','time_only'
useVelAcc = 1; 

allowNonHB = 0; % allow non-hit-bound trials? if set to 0 and a trial lasts
% longer than max_dur, it is discarded. If set to 1, those trials are
% assigned RT = max_dur (affects comparison with mean RT in images_dtb, 
% which is calculated only for bound crossings)


%% PARAMS - these are things we will actually fit!

% drift rate and bound
kmult       = 150;              % drift rate multiplier
kvis        = kmult*cohs;       % assume drift proportional to coh, reduces nParams
kves        = mean(kvis);       % for now, assume 'straddling'
% knoise      = [0.07 0.07];    % additional variability added to drift rate; unused for now
B           = 1.8;              % assume a single bound, but different Tnds for each modality

% diffusion
sigma       = 1;                % unit variance (Moreno-Bote 2010), not a free param!         
sigmaVes    = sigma;            % assume same sigma for all modalities, for now
sigmaVis    = [sigma sigma];    % [at the very least, need to assume their average is 1]

% PDW
theta       = [1.2 0.9 1.05];   % threshold for high bet in logOdds [ves vis comb]
alpha       = 0.03;             % base rate of low bets (offset to PDW curve, as seen in data)
Tconf       = 0;                % (ms), delay between choice and conf report (not implemented yet)

 % Tnd = non-decision time (ms), to account for sensory/motor latencies
TndMean     = [300 500 400];    % must have different Tnds for [ves, vis, comb]
TndSD       = [0 0 0];          % 50-100 works well; set to 0 for fixed Tnd 
TndMin      = TndMean/2;        % need to truncate the Tnd dist
TndMax      = TndMean+TndMin;


%% store the original set parameters, to use e.g. for (pre)param recovery
origParams.kmult = kmult;
origParams.kvis  = kvis;
origParams.kves  = kves;
origParams.B     = B;
origParams.sigmaVes = sigmaVes;
origParams.sigmaVis = sigmaVis;
origParams.theta = theta;
origParams.alpha = alpha;
origParams.TndMean = TndMean;
origParams.TndSD   = TndSD;
origParams.TndMin  = TndMin;
origParams.TndMax  = TndMax;
origParams.Tconf   = Tconf;