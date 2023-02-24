function [err,fit] = dots3DMP_fit_2Dacc_err_optk_sepbounds(param, guess, fixed, data, options)

global call_num

% set parameters for this run based on guess and fixed params flag
param = getParam(param, guess, fixed);

if isfield(data,'PDW')
    data.conf = data.PDW;
end

mods   = unique(data.modality)';
cohs   = unique(data.coherence)'; 
hdgs   = unique(data.heading)';
% deltas = unique(data.delta)';
deltas = 0;

% don't just use trials from data to run simulations, use a fixed, large amount!
nreps = options.nreps;
[hdg, modality, coh, delta, ntrials] = dots3DMP_create_trial_list(hdgs,mods,cohs,deltas,nreps,0); % don't shuffle

% modality = data.modality;
% coh      = data.coherence;
% hdg      = data.heading;
% delta    = data.delta;

duration = 2000; % stimulus duration (ms)

kves  = param(1);
kvis  = param([2 3]);
sigmaVes = param(4);
sigmaVis = param([5 6]);
BVes     = abs(param(7)); % don't accept negative bound heights
BVis     = abs(param(8)); % fixed across cohs
BComb    = abs(param(9));
muTndVes = param(10);
muTndVis = param(11); % fixed across cohs
muTndComb = param(12);

if numel(param)==13, theta = param(13); end % only relevant for PDW

if options.conftask==2
    timeToConf = 350; % additional processing time for confidence
else
    timeToConf = 0;
end
duration = duration + timeToConf;

dur = ones(ntrials,1) * duration; % fixed for all sims

sdTnd = 60; % fixed SD
% Tnds = muTnd + randn(ntrials,1).*sdTnd; % fixed for all sims of given trial

% assume the mapping is based on an equal amount of experience with the 
% *three* levels of reliability (ves, vis-low, vis-high) hence k and sigma
% are their averages
k = mean([kves kvis]);

RVes.t = 0.001:0.001:duration/1000;
RVes.Bup = BVes;
RVes.drift = k * sind(hdgs(hdgs>=0)); % takes only unsigned drift rates
RVes.lose_flag = 1;
RVes.plotflag = 0; % 1 = plot, 2 = plot and export_fig

PVes =  images_dtb_2d(RVes);

RVis.t = 0.001:0.001:duration/1000;
RVis.Bup = BVis;
RVis.drift = k * sind(hdgs(hdgs>=0)); % takes only unsigned drift rates
RVis.lose_flag = 1;
RVis.plotflag = 0; % 1 = plot, 2 = plot and export_fig

PVis =  images_dtb_2d(RVis);

RComb.t = 0.001:0.001:duration/1000;
RComb.Bup = BComb;
RComb.drift = k * sind(hdgs(hdgs>=0)); % takes only unsigned drift rates
RComb.lose_flag = 1;
RComb.plotflag = 0; % 1 = plot, 2 = plot and export_fig

PComb =  images_dtb_2d(RComb);

% create acceleration and velocity profiles (arbitrary for now)
% SJ 04/2020
% Hou et al. 2019, peak vel = 0.37m/s, SD = 210ms

ampl = 0.16; % movement in metres
sigma = 0.14; %
pos = normcdf(1:duration,duration/2,sigma*duration)*ampl;
vel = gradient(pos)*1000; % meters/s
acc = gradient(vel);

% vel = normpdf(1:duration,duration/2,210);
% vel = 0.37*vel./max(vel);
% acc = gradient(vel)*1000; % multiply by 1000 to get from m/s/ms to m/s/s

% normalize
vel = vel./max(vel);
acc = abs(acc./max(acc)); % (and abs)


%% bounded evidence accumulation

choice          = nan(ntrials,1);
RT              = nan(ntrials,1);
finalV          = nan(ntrials,1);
hitBound        = zeros(ntrials,1);
logOddsCorr     = nan(ntrials,1);
expectedPctCorr = nan(ntrials,1);
conf            = nan(ntrials,1);

% ME is now a draw from bivariate normal with mean vector Mu and covariance matrix V
% start with correlation matrix:
S = [1 -1/sqrt(2) ; -1/sqrt(2) 1];
% -1/sqrt(2) is the correlation for our version of images_dtb;
% can only be changed with an update to the flux file

for n = 1:ntrials
%     Tnd = Tnds(n) / 1000; % non-decision time for nth trial in seconds

    switch modality(n)
        case 1
            Tnd = (muTndVes + randn.*sdTnd)/1000;
            B = BVes; P = PVes; R = RVes;

            mu = acc .* kves * sind(hdg(n)) / 1000; % mean of momentary evidence
                % (I'm guessing drift rate in images_dtb is per second, hence div by 1000)
            s = [sigmaVes sigmaVes]; % standard deviaton vector (see below)
            
        case 2
            Tnd = (muTndVis + randn.*sdTnd)/1000;
            B = BVis; P = PVis; R = RVis;
            
            mu = vel .* kvis(cohs==coh(n)) * sind(hdg(n)) / 1000;
            s = [sigmaVis(cohs==coh(n)) sigmaVis(cohs==coh(n))];
        case 3
            Tnd = (muTndComb + randn.*sdTnd)/1000;
            B = BComb; P = PComb; R = RComb;
            
            % positive delta defined as ves to the left, vis to the right
            muVes = acc .* kves               * sind(hdg(n)-delta(n)/2) / 1000;
            muVis = vel .* kvis(cohs==coh(n)) * sind(hdg(n)+delta(n)/2) / 1000;

            % optimal weights (Drugo et al.) 
            wVes = sqrt( kves^2 / (kvis(cohs==coh(n))^2 + kves^2) );
            wVis = sqrt( kvis(cohs==coh(n))^2 / (kvis(cohs==coh(n))^2 + kves^2) );
            mu = wVes.*muVes + wVis.*muVis;
            
            % the DV is a sample from a dist with mean = weighted sum of
            % means. thus the variance is the weighted sum of variances
            % (error propagation formula):
            sigmaComb = sqrt(wVes.^2 .* sigmaVes^2 + wVis.^2 .* sigmaVis(cohs==coh(n))^2); % assume zero covariance
            s = [sigmaComb sigmaComb];
    end

    Mu = [mu; -mu]';
    V = diag(s)*S*diag(s);% convert correlation to covariance matrix
    dv = [0 0; cumsum(mvnrnd(Mu,V))];
    % because Mu is signed according to heading (positive=right),
    % dv(:,1) corresponds to evidence favoring rightward, not evidence
    % favoring the correct decision (as in Kiani eqn. 3 and images_dtb)
    
    % decision outcome
    cRT1 = find(dv(1:dur(n)-timeToConf,1)>=B, 1);
    cRT2 = find(dv(1:dur(n)-timeToConf,2)>=B, 1);
    % the options are:
    % (1) only right accumulator hits bound,
    if ~isempty(cRT1) && isempty(cRT2)
        RT(n) = cRT1/1000 + Tnd;
        finalV(n) = dv(cRT1+timeToConf,2); % only 1 hit, so 2 is the loser
        hitBound(n) = 1;
        choice(n) = 1;
    % (2) only left accumulator hits bound, 
    elseif isempty(cRT1) && ~isempty(cRT2)
        RT(n) = cRT2/1000 + Tnd;
        finalV(n) = dv(cRT2+timeToConf,1); % only 2 hit, so 1 is the loser
        hitBound(n) = 1;
        choice(n) = -1;
    % (3) neither hits bound,
    elseif isempty(cRT1) && isempty(cRT2)
        RT(n) = dur(n)/1000 + Tnd;
            % this is interesting: which DV matters for confidence
            % if neither hits bound? for now take the (abs) maximum,
            % ie whichever was closest to hitting bound. (alternative
            % would be their average?)
            % SJ 07/2020 finalV is technically distance of loser from bound (when winner
            % hits), so in this case, should also account for where winner is wrt bound 
        whichWon = dv(dur(n)-timeToConf,:)==max(dv(dur(n)-timeToConf,:));
        
%         finalV(n) = dv(end,~whichWon); % the not-whichWon is the loser
        % % finalV(n) = mean(dvEnds); 
        finalV(n) = dv(end,~whichWon) + B-dv(end,whichWon); % 

        hitBound(n) = 0;
        a = [1 -1];
        choice(n) = a(whichWon);
    % (4) or both do
    else
        RT(n) = min([cRT1 cRT2])/1000 + Tnd;
        whichWon = [cRT1<=cRT2 cRT1>cRT2];
        finalV(n) = dv(min([cRT1 cRT2]+timeToConf),~whichWon); % the not-whichWon is the loser
        hitBound(n) = 1;
        a = [1 -1];
        choice(n) = a(whichWon);
    end
        
    % use map to look up log-odds that the motion is rightward
    diffV = abs((P.y+B)-finalV(n));
    diffT = abs(R.t-RT(n));
        
    thisV = find(diffV==min(diffV));
    thisT = find(diffT==min(diffT));
    
    switch options.confModel
        case 'evidence+time'
            % use map to look up log-odds that the motion is rightward
            logOddsCorr(n) = P.logOddsCorrMap(thisV(1), thisT(1));
            
            if ~exist('theta','var') % sacc endpoint
                expectedPctCorr(n) = logistic(logOddsCorr(n)); % convert to pct corr
                conf(n) = 2*expectedPctCorr(n) - 1; % convert to 0..1
            else % PDW
                conf(n) = logOddsCorr(n) > theta;
            end
        case 'evidence_only'
            if ~exist('theta','var')
                conf(n) = max(diffV) ./ range(P.y);
            else
                conf(n) = max(diffV) > theta;
            end
        case 'time_only'
            if ~exist('theta','var')
                conf(n) = 1 - RT(n) ./ range(P.t);
            else
                conf(n) = RT(n) < theta;
            end
    end
end

choice(choice==0) = sign(randn); % not needed under usual circumstances
choice(choice==1) = 2; choice(choice==-1) = 1; % 1=left, 2=right

% output var 'fit' gets the same values as data for the conditions, but the
% simulated trial outcomes for the observables!
fit.heading   = hdg;
fit.coherence = coh;
fit.modality  = modality;
fit.delta     = delta;
fit.choice    = choice;
fit.RT        = RT;
fit.conf      = conf;

fit.correct = ((choice==2 & hdg>0) | (choice==1 & hdg<0)) | ((hdg==0 | abs(hdg)<abs(delta)) & rand(size(hdg))<0.5);
data.correct = (data.choice==2 & data.heading>0) | (data.choice==1 & data.heading<0) | ...
    ((data.heading==0 | abs(data.heading)<abs(data.delta)) & rand(size(data.heading))<0.5);

%% calculate error


RTdata = data.RT .* 1000;
RTfit  = mean(fit.RT,2)  .* 1000;

if options.conftask==1
    conffit = mean(fit.conf,2) .* 100;
    confdata = data.conf .* 100;
else 
    conffit = fit.conf;
    confdata = data.conf;
end

% convert to 0/1
choiceD  = logical(data.choice-1);
choiceM  = logical(fit.choice-1);

n = nan(length(mods),length(cohs),length(deltas),length(hdgs));

meanRT_fit = n;     meanRT_data = n;    sigmaRT = n;
meanConf_fit = n;   meanConf_data = n;  sigmaConf = n;
nCor = n;
nmodel = n;

% nRight_data   = n;
% nHighBet_data = n;
pRight_model = nan(length(data.heading),1);
pHigh_model  = nan(length(data.heading),1);

for m = 1:length(mods)
for c = 1:length(cohs)
for d = 1:length(deltas)

    for h = 1:length(hdgs)
        Jdata  = data.modality==mods(m) & data.coherence==cohs(c) & data.heading==hdgs(h) & data.delta==deltas(d);
        Jmodel = modality==mods(m) & coh==cohs(c) & hdg==hdgs(h) & delta==deltas(d);
        n(m,c,d,h) = sum(Jdata);
        nmodel(m,c,d,h) = sum(Jmodel);

        % pRight_model stores the predictions made by the model for the
        % probability of a rightward choice on each trial condition in the
        % data
        % and then assigned to the actual data trials of each condition in pRight_model for
        % log likelihood calculation
        
        if m==1, Ptemp = PVes; muTnd = muTndVes;
        elseif m==2, Ptemp = PVis; muTnd = muTndVis;
        elseif m==3, Ptemp = PComb; muTnd = muTndComb;
        end
        
%         if hdgs(h)<0
%             % if hdg is leftward, then pRight is p(incorrect response)
%             pRight_model(Jdata) = Ptemp.lo.p(abs(hdgs(h))==(hdgs(hdgs>=0)));  
%         else
%             pRight_model(Jdata) = Ptemp.up.p(abs(hdgs(h))==(hdgs(hdgs>=0)));  
%         end
        pRight_model(Jdata) = sum(Jmodel & fit.choice==2)/nmodel(m,c,d,h); % 2 is rightward
        
        if options.conftask==2
            pHigh_model(Jdata) = sum(Jmodel & conffit==1)/nmodel(m,c,d,h); % 1 is high
        end
        
        % RT and SEP fits
        % use only correct trials (or 0 hdg)
        usetrs_data  = data.correct | data.heading==0;
        usetrs_model = fit.correct | hdg==0;
        
        nCor(m,c,d,h) = sum(Jdata & usetrs_data); % use only correct trials for RT and conf (sacc EP) fits
        
        if options.RTtask
            meanRT_fit(m,c,d,h) = mean(RTfit(Jmodel & usetrs_model));
            
%             if hdgs(h)<0
%             % if hdg is leftward, then pRight is p(incorrect response)
%                 meanRT_fit(m,c,d,h) = Ptemp.lo.mean_t(abs(hdgs(h))==(hdgs(hdgs>=0)))+muTnd/1000;  
%             else
%                 meanRT_fit(m,c,d,h) = Ptemp.up.mean_t(abs(hdgs(h))==(hdgs(hdgs>=0)))+muTnd/1000; 
%             end
            
            meanRT_data(m,c,d,h) = mean(RTdata(Jdata & usetrs_data));
            sigmaRT(m,c,d,h) = std(RTdata(Jdata & usetrs_data)) / sqrt(nCor(m,c,d,h));
        end
        
        if options.conftask == 1 % sacc endpoint
            meanConf_fit(m,c,d,h) = mean(conffit(Jmodel & usetrs_model));
            meanConf_data(m,c,d,h) = mean(confdata(Jdata & usetrs_data));
            sigmaConf(m,c,d,h) = std(confdata(Jdata & usetrs_data)) / sqrt(nCor(m,c,d,h));
        end
    end
    
end
end
end


%% calculate log likelihoods

% % likelihood of rightward choice on each trial, under binomial assumptions
% based on Palmer et al. 2005 coh version, but this ignores mod/coh in 3DMPtask
% Pr_model = 1 ./ (1 + exp(-2*k*B*sind(hdg)));

% kluge to avoid log(0) issues
pRight_model(pRight_model==0) = eps; 
pRight_model(pRight_model==1) = 1-eps;

% log likelihood of choice is summed log likelihood across all trials (log
% probability of observing choice given model)
LL_choice = (nansum(log(pRight_model(choiceD))) + nansum(log(1-pRight_model(~choiceD))));

% log likelihood of mean RTs, under Gaussian approximation
LL_RT = 0;
if options.RTtask
    L_RT = 1./(sigmaRT.*sqrt(2*pi)) .* exp(-(meanRT_fit - meanRT_data).^2 ./ (2*sigmaRT.^2));
    L_RT(L_RT==0) = eps; %min(L_RT(L_RT~=0));
    LL_RT = nansum(log(L_RT(:))); % sum over all conditions
end

% likelihood of mean confidence (for SEP), or log probabilities for PDW
LL_conf = 0; % initialize, if no conf task will stay as 0
if options.conftask == 1 % SEP
    L_conf = 1./(sigmaConf*sqrt(2*pi)) .* exp(-(meanConf_fit - meanConf_data).^2 ./ (2*sigmaConf.^2));
    L_conf(L_conf==0) = eps; %min(L_conf(L_conf~=0));
    LL_conf = nansum(log(L_conf(:))); % sum over all conditions
    
elseif options.conftask == 2 % PDW
    pHigh_model(pHigh_model==0) = eps; 
    pHigh_model(pHigh_model==1) = 1-eps;
    PDW = logical(data.conf);
    LL_conf = (nansum(log(pHigh_model(PDW))) + nansum(log(1-pHigh_model(~PDW))));   
end


err = -(LL_choice + LL_conf + LL_RT); % negate for minimization
% err = -(LL_choice + LL_conf);
% err = -LL_choice;

%% print progress report!
fprintf('\n\n\n****************************************\n');
fprintf('run %d\n', call_num);
fprintf('param values\n')
for i=1:length(param)
    fprintf('\t%g\t',param(i))
end
fprintf('\n');

fprintf('err: %f\n', err);
fprintf('split errs: Choice: %.2f, Conf: %.2f, RT: %.2f\n',LL_choice,LL_conf,LL_RT)
if options.ploterr && strcmp(options.fitMethod,'fms')
    figure(options.fh); hold on
    plot(call_num, err, '.','MarkerSize',12,'color','k');
    drawnow;
end
call_num = call_num + 1;


end


% retrieve the full parameter set given the adjustable and fixed parameters 
function param2 = getParam ( param1 , guess , fixed )
  
  if all(fixed), param2 = param1; return; end

  param2(fixed==0) = param1;            %get adjustable parameters from param1
  param2(fixed==1) = guess(fixed==1);   %get fixed parameters from guess (initial point)
end
