% simulate "CCC" (cue-combination + confidence) experiment
% specifically in the case of the vis-ves heading task

% CF started it 2016
% heavily modified early 2019
% switched to 2D accumulator model 06/2020
% modified again 11/2022

%% build expt and hand-pick some model params

% at some point, build these settings out into a proper wrapper script

clear; close all

% all these folder settings are user specific and shouldn't be in the script

% datafolder = '/Users/chris/Documents/MATLAB';
% codefolder = '/Users/chris/Documents/MATLAB/Projects/offlineTools/dots3DMP/behav_models';

datafolder = '/Users/stevenjerjian/Desktop/FetschLab/Analysis/data/dots3DMP_DDM';
codefolder = '/Users/stevenjerjian/Desktop/FetschLab/Analysis/codes/FLprojects/dots3DMP/behav_models';

cohs = [0.4 0.8]; % visual coherence levels (these are really just labels, since k's are set manually)
% hdgs = [-12 -6 -3 -1.5 -eps eps 1.5 3 6 12]; % don't know if we realy need two zeroes
hdgs = [-12 -6 -3 -1.5 0 1.5 3 6 12];
% deltas = [-3 0 3]; % conflict angle; positive means vis to the right
deltas = 0; % conflict angle; positive means vis to the right
mods = [1 2 3]; % stimulus modalities: ves, vis, comb

%% MODEL SPECIFICATIONS and CONDITIONS

% === model type ==== unused for now
modelID = 1; % 1 will be 2Dacc model ('Candidate model' against which others are tested).
% modelVar = 1; % variation within model ID, e.g. velocity acceleration coding, confidence model, cue weighting in combined..

% task type
RTtask   = 1;
conftask = 2; % 1 - sacc endpoint, 2 - PDW

nreps = 50; % number of repetitions of each unique trial type % (ntrials depends on num unique trial types)

% stimulus conditions
mods  = [1 2 3];        % stimulus modalities: ves, vis, comb
cohs  = [0.4 0.8];      % visual coherence levels (these are really just labels, since k's are set manually)
hdgs  = [-12 -6 -3 -1.5 0 1.5 3 6 12];
% deltas = [-3 0 3];    % conflict angle; positive means vis to the right
deltas  = 0;

% time information
dT      = 1; % time step, ms
max_dur = 2100; % stimulus duration (ms)

% maybe these become supplanted by modelVar eventually
confModel  = 'evidence+time'; % 'evidence+time','evidence_only','time_only'
useVelAcc  = 0; 
allowNonHB = 0; % allow non-hit-bound trials? if set to 0 and a trial lasts
% longer than max_dur, it is discarded. If set to 1, those trials are
% assigned RT = max_dur (affects comparison with mean RT in images_dtb, 
% which is calculated only for bound crossings)

urgency = 1;
plotSimulatedData = 1; % plot average psychometric curves etc.


%% generative model parameters

% drift rate and bound
kmult       = 25;               % drift rate multiplier
kvis        = kmult*cohs;       % assume drift proportional to coh, reduces nParams
kves        = mean(kvis);       % for now, assume 'straddling'
% knoise      = [0.07 0.07];    % additional variability added to drift rate, unused for now
B           = 2;              % assume a single bound, but different Tnds for each modality

% diffusion
sigma       = 1;                % unit variance (Moreno-Bote 2010), not a free param!         
sigmaVes    = sigma;            % assume same sigma for all modalities, for now
sigmaVis    = [sigma sigma];    % [at the very least, need to assume their average is 1]

% PDW
theta       = [0.8 0.7 1.0];    % threshold for high bet in logOdds [ves vis comb]
alpha       = 0.03;             % base rate of low bets (offset to PDW curve, as seen in data)
% Tconf       = 0;              % (ms), delay between choice and conf report, unused for now

 % Tnd = non-decision time (ms), to account for sensory/motor latencies
TndMean     = [300 300 300];    % must have different Tnds for [ves, vis, comb]
TndSD       = [0 0 0];          % 50-100 works well; set to 0 for fixed Tnd 
TndMin      = TndMean/2;        % need to truncate the Tnd dist
TndMax      = TndMean+TndMin;

%% build trial list

[hdg, modality, coh, delta, ntrials] = dots3DMP_create_trial_list(hdgs,mods,cohs,deltas,nreps,0); % don't shuffle
dur = ones(ntrials,1) * max_dur;

%% create acceleration and velocity profiles

if useVelAcc
    % SJ 04/2020
    % Hou et al. 2019, peak vel = 0.37m/s, SD = 210ms
    % 07/2020 lab settings...16cm in 1.3s, sigma=0.14
    
    % our theoretical settings (before tf)
    ampl = 0.16; % movement in metres
    pos = normcdf(1:max_dur,max_dur/2,140)*ampl;
    vel = gradient(pos); % metres/s
    acc = gradient(vel); 

    % normalize (by max or by mean?) and take abs of acc 
    vel = vel/mean(vel);
    acc = abs(acc)/mean(abs(acc));
%     vel = vel/max(vel);
%     acc = acc/max(abs(acc));
%     %acc(acc<0) = 0;
%     acc = abs(acc);

% step functions, for testing param recovery and LL comparisons
%     acc = [zeros(max_dur/2,1); ones(max_dur/2,1)];
%     vel = [ones(max_dur/2,1); zeros(max_dur/2,1)];

    if useVelAcc==1 % ves follow acc, vis follows vel
        sves = acc; svis = vel;
    elseif useVelAcc==2 % both vel
        sves = vel; svis = vel;
    elseif useVelAcc==3 % both acc
        sves = acc; svis = acc; 
    end
else % or fixed, i.e. no vel/acc weighting
    sves = ones(1,max_dur);
    svis = sves;
end

%% store the generative parameters, to use e.g. for (pre)param recovery

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
% origParams.Tconf   = Tconf;


%% calculate log odds corr map using Wolpert Method of Images code (van den Berg 2016, after Moreno-Bote 2010)

% assume the mapping is based on an equal amount of experience with the 
% *three* levels of reliability (ves, vis-low, vis-high) hence k and sigma
% are their averages, for the purpose of expected logOddsCorr
    %... but only if all three mods are present!
if all(mods==1)
    k = kves;
else
    k = mean([kves kvis]);
end

R.t = (dT/1000:dT/1000:max_dur/1000)';
R.Bup = B;
R.lose_flag = 1;
R.plotflag  = 0; % 1 = plot, 2 = plot and export_fig

% new for MOIcollapse

if urgency
    R.grid   = linspace(-4*R.Bup,0,500); % this is now an argument in MOIcollapse
    R.k_urg  = 1;
    R.low_th = -R.Bup-R.Bup/4;
    
%     if useVelAcc
        images_func = @images_dtb_2d_urg_varDrift;
%     else
%         images_func = @images_dtb_2d_urg;
%     end
else
%     if useVelAcc
        images_func = @images_dtb_2d_varDrift;
%     else
%         images_func = @images_dtb_2d;
%     end
end

% R.drift = k * sind(hdgs(hdgs>=0)); % takes only unsigned drift rates
% P = images_dtb_2d(R);


% SJ 02-2023 separate logOdds map for each modality, variable drift rates

R.drift = sves .* kves .* sind(hdgs(hdgs>=0))';
P(1) = images_func(R); % ves

R.drift = svis .* mean(kvis) .* sind(hdgs(hdgs>=0))';
P(2)    = images_func(R); % vis

R.drift = sqrt(sves.*kves.^2 + svis.*mean(kvis).^2) .* sind(hdgs(hdgs>=0))';
P(3)    = images_func(R); % comb


%% simulate bounded evidence accumulation

% preallocate
dv_all = cell(ntrials,1); % shouldn't need to store every trial's DV, but if you want to, it's here
choice          = nan(ntrials,1); % choices (left = -1, right = 1);
RT              = nan(ntrials,1); % reaction time (or time-to-bound for fixed/variable duration task)
finalV          = nan(ntrials,1); % now this is the value of the losing accumulator
hitBound        = zeros(1,ntrials); % hit bound or not on that trial
logOddsCorr     = nan(ntrials,1); % log odds correct
expectedPctCorr = nan(ntrials,1); % expected probability correct (converted to confidence rating)
conf            = nan(ntrials,1); % confidence rating
pdw             = nan(ntrials,1); % post-decision wager

tic
for n = 1:ntrials
    
    % momentary evidence is a draw from bivariate normal distribution
    % with mean Mu (2 x time) and covariance matrix V (2x2).
    
    % Start with a desired correlation matrix:
    S = [1 -1/sqrt(2) ; -1/sqrt(2) 1];
    % -1/sqrt(2) is the correlation for our version of images_dtb;
    % can only be changed with an update to the 'flux' file

    % the next step depends on modality:
    switch modality(n)
        case 1
            % assume momentary evidence is proportional to sin(heading) (Drugowitsch14)
            mu = sves * kves * sind(hdg(n)) * dT/1000; % mean of momentary evidence, scaled by delta-t
            
            % convert correlation to covariance matrix:
            % s is the standard deviaton vector,
            s = [sigmaVes*sqrt(dT/1000) sigmaVes*sqrt(dT/1000)]; 
                     %^ variance is sigma^2*dt, so stdev is sigma*sqrt(dt),
                     % after converting from ms to seconds
        case 2
            mu = svis .* kvis(cohs==coh(n)) * sind(hdg(n)) * dT/1000;
            s = [sigmaVis(cohs==coh(n))*sqrt(dT/1000) sigmaVis(cohs==coh(n))*sqrt(dT/1000)];
        case 3
            % positive delta defined as ves to the left, vis to the right
            muVes = sves .* kves               * sind(hdg(n)-delta(n)/2) * dT/1000;
            muVis = svis .* kvis(cohs==coh(n)) * sind(hdg(n)+delta(n)/2) * dT/1000;
            
            % optimal weights (Drugo et al.)
            wVes = sqrt( kves^2 / (kvis(cohs==coh(n))^2 + kves^2) );
            wVis = sqrt( kvis(cohs==coh(n))^2 / (kvis(cohs==coh(n))^2 + kves^2) );
            
            % corrupt optimal weights with noise
%             wVes = wVes + randn*knoise(1);
%             wVis = wVis + randn*knoise(2);

            % or randomize (TMS idea, for Amir grant)
%             wVes = rand; wVis = 1 - wVes;       

            mu = wVes.*muVes + wVis.*muVis;
            
            % the DV is a sample from a dist with mean = weighted sum of
            % means. thus the variance is the weighted sum of variances
            % (error propagation formula):
            sigmaComb = sqrt(wVes.^2 .* sigmaVes^2 + wVis.^2 .* sigmaVis(cohs==coh(n))^2); % assume zero covariance
            s = [sigmaComb*sqrt(dT/1000) sigmaComb*sqrt(dT/1000)];
    end

    Mu = [mu; -mu]'; % mean vector for 2D DV

    % convert correlation to covariance matrix
    V = diag(s)*S*diag(s);
    dv = [0 0; cumsum(mvnrnd(Mu,V))]; % bivariate normrnd, where Mu is a

    dv_all{n} = dv; % uncomment if saving DV

    % decision outcome
    cRT1 = find(dv(1:dur(n),1)>=B, 1);
    cRT2 = find(dv(1:dur(n),2)>=B, 1);
    
    % the options are:
    % (1) only right accumulator hits bound,
    if ~isempty(cRT1) && isempty(cRT2)
        RT(n) = cRT1;
        finalV(n) = dv(cRT1,2); % only 1 hit, so 2 is the loser
        hitBound(n) = 1;
        choice(n) = 1;
    % (2) only left accumulator hits bound,
    elseif isempty(cRT1) && ~isempty(cRT2)
        RT(n) = cRT2;
        finalV(n) = dv(cRT2,1); % only 2 hit, so 1 is the loser
        hitBound(n) = 1;
        choice(n) = -1;
    % (3) neither hits bound,
    elseif isempty(cRT1) && isempty(cRT2)
        hitBound(n) = 0;
        if allowNonHB
            RT(n) = dur(n);

            % which DV matters for confidence if neither hits bound? 
            % SJ 07/2020 logOddsCorrMap is fixed, so just shift finalV up
            % so that 'winner' did hit bound,
            whichWon = dv(dur(n),:)==max(dv(dur(n),:));
            finalV(n) = dv(end,~whichWon) + B-dv(end,whichWon);
            % ^  shifting the losing dv up by whatever the
            % difference is between the bound and the winning dv
            a = [1 -1];
            choice(n) = a(whichWon);
        else
            RT(n) = NaN;
            finalV(n) = NaN;
            choice(n) = NaN;
        end
    % (4) or both do
    else
        RT(n) = min([cRT1 cRT2]);
        whichWon = [cRT1<=cRT2 cRT1>cRT2];
        finalV(n) = dv(min([cRT1 cRT2]),~whichWon); % the not-whichWon is the loser
        hitBound(n) = 1;
        a = [1 -1];
        choice(n) = a(whichWon);
    end
    
    % calculate confidence
    if hitBound(n)==0 && allowNonHB==0
        logOddsCorr(n) = NaN;
        conf(n) = NaN;
        pdw(n) = NaN;
    else    
        Pm = P(modality(n));
        diffV = abs((Pm.y+B)-finalV(n));
        diffT = abs(R.t*1000-RT(n));
        
        switch confModel
            case 'evidence+time'
                % use map to look up log-odds that the motion is rightward
                thisV = find(diffV==min(diffV));
                thisT = find(diffT==min(diffT));
                
                logOddsCorr(n) = Pm.logOddsCorrMap(thisV(1), thisT(1));
                expectedPctCorr(n) = logistic(logOddsCorr(n)); % convert to pct corr
                
                conf(n) = 2*expectedPctCorr(n) - 1; % convert to 0..1
                pdw(n) = logOddsCorr(n) > theta(modality(n)); % threshold
            case 'evidence_only'
                conf(n) = max(diffV) ./ range(Pm.y);
                pdw(n) = max(diffV) > theta(modality(n));
            case 'time_only'
                conf(n) = 1 - (RT(n)/1000) ./ range(Pm.t);
                pdw(n) = (RT(n)/1000) < theta(modality(n));
        end
        if isnan(conf(n)), conf(n)=0; end % if dvs are almost overlapping, force conf to zero as it can sometimes come out as NaN
    end
end
toc


%% post-processing

% quick sanity check that params are reasonable
correct = (choice==1 & hdg>0) | (choice==-1 & hdg<0) | ...
    (rand<0.5 & (hdg==0 | abs(hdg)<abs(delta)));
pCorrect_total = sum(correct) / ntrials

% remove NaNs (ie non-HBs, if not allowed)
remove = isnan(choice);
dv_all(remove) = [];
modality(remove)=[];
hdg(remove)=[];
coh(remove)=[];
delta(remove)=[];
choice(remove)=[];
RT(remove)=[];
conf(remove)=[];
pdw(remove)=[];
correct(remove)=[];
ntrials = length(choice);

% adjust wager probability for base rate of low bets, as seen in data
% ('compresses' the curve, not an offset, because P(high) varies with hdg)

% first, save original PDW for plotting choice/RT splits (bc relationship
% between conf and choice/RT is, by construction, independent of alpha)
pdw_preAlpha = pdw;
    % Puzzling why this is needed, actually.
    % Although we might wonder what those relationships would have looked
    % like if we had access to actual confidence independent of alpha, in
    % practice it is irrelevant because we cannot disentangle this in data.
    % Yet the pre-param recovery misses if we don't use this, in both the
    % 'data' and model calculations...

pdw(pdw==1 & rand(length(pdw),1)<alpha) = 0;

% add non-decision time (truncated normal dist)
Tnd = zeros(ntrials,1);
for n = 1:ntrials
    while Tnd(n)<=TndMin(modality(n)) || Tnd(n)>=TndMax(modality(n)) % simple trick for truncating, requires integers (ms)
        Tnd(n) = round(normrnd(TndMean(modality(n)),TndSD(modality(n))));
    end
end
DT = RT; % rename this 'decision time'
RT = DT+Tnd;

% should be >0.95 or else maxDur isn't long enough (or need urgency signal!)
pHitBound = sum(hitBound)/length(hitBound)
Z = abs(hdg)<0.0001;
pHitBound_zeroHdg = sum(hitBound(Z))/sum(Z)


%% format data as in experimental data files and generate output structs

choice(choice==0) = sign(randn)*eps; % should have no actual zeros, but if so, sign them randomly;
                               % this is just to assign a direction and correct/error
choice(choice==1) = 2; choice(choice==-1) = 1; % 1=left, 2=right

data.modality   = modality;
data.heading    = hdg;
data.coherence  = coh;
data.delta      = delta;
data.choice     = choice;
data.RT         = RT/1000; % change to seconds
data.conf       = conf;
data.PDW        = pdw;
data.PDW_preAlpha = pdw_preAlpha;
data.correct    = correct;
subject         = 'simul';

data.oneTargConf = false(size(data.heading));


%% save it
save(fullfile(datafolder,[savefilename '.mat']), 'data', 'origParams', 'allowNonHB', 'sves', 'svis');

%% plots, optional

if plotSimulatedData
    mods   = unique(data.modality);
    cohs   = unique(data.coherence);
    deltas = unique(data.delta);
    hdgs   = unique(data.heading);
    
    % means per condition, logistic fits
    parsedData = dots3DMP_parseData(data,mods,cohs,deltas,hdgs,conftask,RTtask);
    
    % plot it
    dots3DMP_plots(parsedData,mods,cohs,deltas,hdgs,conftask,RTtask)

    % gaussian fits
%     gfit = dots3DMP_fit_cgauss(data,mods,cohs,deltas,conftask,RTtask);
%     dots3DMP_plots_cgauss_byCoh(gfit,parsedData,mods,cohs,deltas,hdgs,conftask,RTtask)
end

