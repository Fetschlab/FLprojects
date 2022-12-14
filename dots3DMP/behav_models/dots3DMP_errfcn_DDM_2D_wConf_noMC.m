function [err,fit,parsedFit,logOddsMap] = dots3DMP_errfcn_DDM_2D_wConf_noMC(param, guess, fixed, data, options)

% updated model calculations to Steven's method: 
% SJ 10-11-2021 no Monte Carlo simulation for fitting, it's redundant!
% just use model predictions directly

% uses UNSIGNED cohs for MOI calculation
% ONLY FIT zero delta

% SJ spawned from dots version, for dots3DMP

tic

global call_num

% retrieve the full parameter set given the adjustable and fixed parameters 
param = getParam(param, guess, fixed);
if isfield(data,'dur')
    max_dur = max(data.dur)/1000; % max stimulus duration (s)
else
    max_dur = 2.3;
end

kmult = param(1)*100; % scale up
B = abs(param(2)); % don't accept negative bound heights
theta = abs(param(3:5)); % or negative thetas, one theta per mod
alpha = param(6); % base rate of low-conf choices
Tnds = param(7:9); % separate Tnd for each modality
% perhaps this will become obsolete too with acc/vel

cohs = unique(data.coherence);
mods = unique(data.modality);
hdgs = unique(data.heading);

if options.dummyRun == 0
    deltas = 0;                 % fit only zero delta
else
    deltas = unique(data.delta); % predict all deltas
end

% separate kvis and kves
kvis  = kmult*cohs'; % assume drift proportional to coh, reduces nParams
kves  = mean(kvis); % for now, assume 'straddling'

% ...but k for MOI for conf is constant
if all(mods==1), k = kves;
else, k = mean([kves kvis]);
end


if options.useVelAcc
    % SJ 04/2020
    % Hou et al. 2019, peak vel = 0.37m/s, SD = 210ms
    % 07/2020 lab settings...16cm in 1.3s, sigma=0.14
    
    % our theoretical settings (before tf)
    ampl = 0.16; % movement in metres
    pos = normcdf(1:max_dur*1000,max_dur*1000/2,0.14*1000)*ampl;
    vel = gradient(pos)*1000; % metres/s
    acc = gradient(vel); 

    % normalize (by max or by mean?) and take abs of acc 
%     vel = vel/mean(vel);
%     acc = abs(acc)/mean(abs(acc));
    vel = vel/max(vel);
    acc = abs(acc)/max(abs(acc));

    if options.useVelAcc==1 
        sves = acc; svis = vel;
    elseif options.useVelAcc==2 % both vel
        sves = vel; svis = vel;
    elseif options.useVelAcc==3 % both acc
        sves = acc; svis = acc; 
    end
else
    sves = ones(1,max_dur*1000);
    svis = sves;
end

% generate or load log odds correct map

tConf = 0; % placeholder, should be a fit param eventually

% use method of images to calculate PDFs of DV and mapping to log odds corr
R.t = 0.001:0.001:max_dur+tConf;
R.Bup = B;
R.drift = k * sind(hdgs(hdgs>=0)); % takes only unsigned drift rates

R.lose_flag = 1; % we always need the losing densities
R.plotflag = options.plot; % 1 = plot, 2 = plot nicer and export_fig (eg for talk)
% R.plotflag = 1; % manual override
Pconf = images_dtb_2d(R); % /WolpertMOI

% dt = R.t(2)-R.t(1);
% skipT = floor(tConf/dt);
% Pconf.t(1:skipT) = [];
% Pconf.logOddsCorrMap(:,1:skipT) = [];

% use previously calculated logOddsMap if passed in
if options.dummyRun==0
    logOddsMap = Pconf.logOddsCorrMap;
else
    logOddsMap = options.logOddsMap;
end


%% now compute a separate P for model choices and RTs (where drift rate is modality-specific)

% ves
RVes.t = 0.001:0.001:max_dur;
RVes.Bup = B;
RVes.drift = sves .* kves .* sind(hdgs(hdgs>=0));
RVes.driftSigned = kves * sind(hdgs); % sves is not needed here

RVes.lose_flag = 1;
RVes.plotflag = 0;
PVes =  images_dtb_2d_varDrift(RVes);

% vis and comb, separate drift for each coh
for c = 1:length(cohs)
    RVis.t = 0.001:0.001:max_dur;
    RVis.Bup = B;
    RVis.drift = svis .* kvis(c) .* sind(hdgs(hdgs>=0));
    RVis.driftSigned = kvis(c) * sind(hdgs);

    RVis.lose_flag = 1;
    RVis.plotflag = 0;
    PVis(c) =  images_dtb_2d_varDrift(RVis);

    RComb.t = 0.001:0.001:max_dur;
    RComb.Bup = B;
    RComb.lose_flag = 1;
    RComb.plotflag = 0;

    if options.dummyRun == 0
      
        kcomb(c) = sqrt(kves.^2 + kvis(c).^2); % optimal per Drugo
        RComb.driftSigned = kcomb(c) * sind(hdgs);

        kcombT(c,:) = sqrt(sves.*kves.^2 + svis.*kvis(c).^2); % optimal per Drugo
        RComb.drift = kcombT(c,:) .* sind(hdgs(hdgs>=0));
        PComb(c,1) =  images_dtb_2d_varDrift(RComb);

    else     % compute w and mu to capture biases under cue conflict
        
        wVes = sqrt( kves^2 / (kves^2 + kvis(c)^2) );
        wVis = sqrt( kvis(c)^2 / (kves^2 + kvis(c)^2) );
        
        for d = 1:length(deltas)
            % positive delta defined as ves to the left, vis to the right
            muVes = kves    * sind(hdgs-deltas(d)/2);
            muVis = kvis(c) * sind(hdgs+deltas(d)/2);

            RComb.driftSigned = wVes.*muVes + wVis.*muVis;
            RComb.drift = abs(sves.*wVes.*muVes + svis.*wVis.*muVis);

%             RComb.drift = abs(RComb.driftSigned);
            PComb(c,d) =  images_dtb_2d_varDrift(RComb);
        end
    end
end


%% calculate likelihood of the observed data given the model parameters

% usetrs_data = data.correct | data.coherence<1e-6; % use only correct (OR ~0 hdg) trials for RT fits
usetrs_data = true(size(data.correct)); % try all trials

% (for parsedFit, and some likelihood calculations)
n = nan(length(mods),length(cohs),length(deltas),length(hdgs));
pRight_model        = n;
pHigh_model         = n;
pRightHigh_model    = n;
pRightLow_model     = n;
pLeftHigh_model     = n;
pLeftLow_model      = n;
pHighCorr_model     = n;
pHighErr_model      = n;
meanConf_model      = n;
meanRT_model        = n;
meanRThigh_model    = n;
meanRTlow_model     = n;

meanRT_data         = n;
meanRThigh_data     = n;
meanRTlow_data      = n;
sigmaRT_data        = n;
sigmaRThigh_data    = n;
sigmaRTlow_data     = n;
meanConf_data       = n;
sigmaConf_data      = n;
nCor                = n;

nTrials             = n;

% trialwise vars
m = nan(length(data.heading),1);
pRight_model_trialwise      = m;
pHigh_model_trialwise       = m;
pRightHigh_model_trialwise  = m;
pRightLow_model_trialwise   = m;
pLeftHigh_model_trialwise   = m;
pLeftLow_model_trialwise    = m;
conf_model_trialwise        = m;
RT_model_trialwise          = m;
RThigh_model_trialwise      = m;
RTlow_model_trialwise       = m;
sigmaRT_data_trialwise      = m;
sigmaRThigh_data_trialwise  = m;
sigmaRTlow_data_trialwise   = m;

for m = 1:length(mods)

    if mods(m)==1, P = PVes; end
    Tnd = Tnds(m);

for c = 1:length(cohs) 


for d = 1:length(deltas)
    
    if deltas(d)~=0 && mods(m)~=3, continue, end

    if mods(m)==2, P = PVis(c);
    elseif mods(m)==3, P = PComb(c,d);
    end


for h = 1:length(hdgs)

%     if deltas(d)==-3 && mods(m)==3, keyboard, end
    % index for images_dtb (unsigned hdg/drift corresponding to this signed hdg)                    
    
    Jdata  = data.modality==mods(m) & data.coherence==cohs(c) & data.heading==hdgs(h) & data.delta==deltas(d);
    nTrials(m,c,d,h) = sum(Jdata);

    % empirical probabilities from data (actually just take counts
    % 'successes', for pdf functions
    nRight_data(m,c,d,h) = sum(Jdata & data.choice==1); % choice is logical 0..1!!
    
    if options.conftask==2
        nHigh_data(m,c,d,h)  = sum(Jdata & data.PDW==1); 

        nRightHigh_data(m,c,d,h) = sum(Jdata & data.choice==1 & data.PDW==1);
        nRightLow_data(m,c,d,h)  = sum(Jdata & data.choice==1 & data.PDW==0);
        nLeftHigh_data(m,c,d,h)  = sum(Jdata & data.choice==0 & data.PDW==1);
        nLeftLow_data(m,c,d,h)   = sum(Jdata & data.choice==0 & data.PDW==0);
    end

    if options.dummyRun==1 && mods(m)==3
        uh = h;
    else
        uh = abs(hdgs(h))==(hdgs(hdgs>=0));  % *** marking differences in signed vs. unsigned ver ***
    end

    % grab the distributions from images_dtb:
    Pxt1 = squeeze(P.up.distr_loser(uh,:,:))'; % losing race pdf given correct
    Pxt2 = squeeze(P.lo.distr_loser(uh,:,:))'; % losing race pdf given incorrect
   
    % CHOICE & PDW
    % yes, we can fit non-zero bias, just need to use the bias-corrected signed
    % drift rate, rather than the experimenter-defined headings!

    if P.driftSigned(h)<0 % leftward motion
        % if hdg is negative, then pRight is p(incorrect), aka P.lo
        pRight = P.lo.p(uh)/(P.up.p(uh)+P.lo.p(uh)); % normalized by total bound-crossing probability 
        pLeft = 1-pRight; % (evidently, the unabsorbed probability can be ignored when it comes to pRight)

        % calculate probabilities of the four outcomes: right/left x high/low
        % these are the intersections, e.g. P(R n H)
        pRightHigh = sum(sum(Pxt2.*(logOddsMap>=theta(m))));
        pRightLow = sum(sum(Pxt2.*(logOddsMap<theta(m))));
        pLeftHigh = sum(sum(Pxt1.*(logOddsMap>=theta(m))));
        pLeftLow =  sum(sum(Pxt1.*(logOddsMap<theta(m))));                                                             
    
    else % rightward motion (and zero motion is symmetrical so gets taken care of here too)

        % if hdg is positive, then pRight is p(correct), aka P.up
        pRight = P.up.p(uh)/(P.up.p(uh)+P.lo.p(uh)); % normalized by total bound-crossing probability
        pLeft = 1-pRight;

        % calculate probabilities of the four outcomes: right/left x high/low
        % these are the intersections, e.g. P(R n H)
        pRightHigh = sum(sum(Pxt1.*(logOddsMap>=theta(m))));
        pRightLow = sum(sum(Pxt1.*(logOddsMap<theta(m))));
        pLeftHigh = sum(sum(Pxt2.*(logOddsMap>=theta(m))));
        pLeftLow =  sum(sum(Pxt2.*(logOddsMap<theta(m))));                                                     
    end

    %ensure total prob sums to one (when it doesn't, it's because of unabsorbed prob)
    % THIS STEP WASN'T THOUGHT TO BE NECESSARY FOR PARAM RECOV, SO BE
    % CAREFUL AND REMOVE IT IF IT MUCKS IT UP
    Ptot = pRightHigh + pRightLow + pLeftHigh + pLeftLow;
    pRightHigh = pRightHigh/Ptot;
    pRightLow  = pRightLow/Ptot;
    pLeftHigh  = pLeftHigh/Ptot;
    pLeftLow   = pLeftLow/Ptot;

    pRightHigh_model(m,c,d,h) = pRightHigh;
    pRightLow_model(m,c,d,h) = pRightLow;
    pLeftHigh_model(m,c,d,h) = pLeftHigh;
    pLeftLow_model(m,c,d,h) = pLeftLow;

    % copy intersections to trials, for multinomial likelihood
    pRightHigh_model_trialwise(Jdata) = pRightHigh;
    pRightLow_model_trialwise(Jdata) = pRightLow;
    pLeftHigh_model_trialwise(Jdata) = pLeftHigh;
    pLeftLow_model_trialwise(Jdata) = pLeftLow;
    
    % by definition of conditional probability: P(A|B) = P(A n B) / P(B)
    pHigh_Right = pRightHigh / pRight;
    pLow_Right = pRightLow / pRight; 
    pHigh_Left = pLeftHigh / pLeft;
    pLow_Left = pLeftLow / pLeft;  
    
    % by law of total probability
    pHigh = pHigh_Right*pRight + pHigh_Left*pLeft;
    pLow = pLow_Right*pRight + pLow_Left*pLeft; % why isn't this 1-pHigh? probability is so weird.

    % by Bayes' rule:
    pRight_High = pHigh_Right * pRight / pHigh;
    pRight_Low = pLow_Right * pRight / pLow;
    pLeft_High = pHigh_Left * pLeft / pHigh;
    pLeft_Low = pLow_Left * pLeft / pLow; %#ok<NASGU>

    % adjust the probabilities for the base rate of low-conf bets (alpha)
    pHigh_wAlpha = pHigh - alpha*pHigh;

    % Lastly, the PDW split; use Bayes to adjust P(H|R) and P(H|L)
    pHigh_Right_wAlpha = pRight_High * pHigh_wAlpha / pRight;
    pHigh_Left_wAlpha = pLeft_High * pHigh_wAlpha / pLeft;
    if hdgs(h)<0 % leftward motion
        pHigh_Corr = pHigh_Left_wAlpha;
        pHigh_Err = pHigh_Right_wAlpha;
    else % rightward motion
        pHigh_Corr = pHigh_Right_wAlpha;
        pHigh_Err = pHigh_Left_wAlpha;    
    end
    
    % copy to vectors for parsedFit struct
    pRight_model(m,c,d,h) = pRight;
    pHigh_model(m,c,d,h) = pHigh_wAlpha; % includes the effect of alpha
%     pRightHigh_model(m,c,d,h) = pRight_High; % does NOT include the effect of alpha
%     pRightLow_model(m,c,d,h) = pRight_Low; % does NOT include the effect of alpha
    
%     pLeftHigh_model(m,c,d,h) = pLeftHigh; % MAYBE TEMP: SEE LIKELIHOOD
%     pLeftLow_model(m,c,d,h) = pLeftLow;  % MAYBE TEMP: SEE LIKELIHOOD
    
    pHighCorr_model(m,c,d,h) = pHigh_Corr; % includes the effect of alpha
    pHighErr_model(m,c,d,h) = pHigh_Err; % includes the effect of alpha

    % copy to trials for fit struct, and likelihood
    pRight_model_trialwise(Jdata) = pRight_model(m,c,d,h);
    pHigh_model_trialwise(Jdata) = pHigh_model(m,c,d,h);
    pRightHigh_model_trialwise(Jdata) = pRightHigh_model(m,c,d,h);
    pRightLow_model_trialwise(Jdata) = pRightLow_model(m,c,d,h);
% % %     pHighCorr_model_trialwise(Jdata) = pHighCorr_model(m,c,d,h);
% % %     pHighErr_model_trialwise(Jdata) = pHighErr_model(m,c,d,h);
        
    % RT
%     nCor(c) = sum(Jdata & (data.correct | data.coherence<1e-6));
    if options.RTtask
        
        % different Tnds for different choices, NOT for different cohs,
        % so we simply take a weighted average based on P(right)
%         Tnd = TndR*pRight + TndL*(1-pRight);
        
        % for param-recov, model RT needs adjusted to account for unabsorbed
        % probability (trials that don't hit the bound, usually only at low
        % cohs, when bound is high, etc). Turns out the calculation of
        % P.up.mean_t only applies to absorbed prob (bound-crossings), so
        % we simply adjust it by averaging mean_t with maxdur, weighted by
        % one minus the probability of bound crossing
            % (this was a cool exercise when it worked, but seems easier to
            % exclude those trials and pretend there is no unabsorbed prob;
            % also more consistent with real behavior where there's always a
            % choice before max_dur is reached  --thanks to MVL for this)
        if options.ignoreUnabs
            pHB = 1;
        else
            pHB = P.up.p(uh)+P.lo.p(uh); % prob of crossing either bound
        end
        
        meanRT_model(m,c,d,h) = pHB*P.up.mean_t(uh) + (1-pHB)*max_dur + Tnd; % weighted average
 
        % for RT conditioned on wager, it seems we don't need to do any
        % scaling by Ptb; the correct distribution is captured by the
        % conditional Pxts, we just sum them, dot product w t and normalize
        PxtAboveTheta = sum(Pxt1.*(logOddsMap>=theta(m))); % shouldn't matter if Pxt1 or 2, it's symmetric
        PxtBelowTheta = sum(Pxt1.*(logOddsMap<theta(m)));
        meanRThigh = PxtAboveTheta * R.t' / sum(PxtAboveTheta);
        meanRTlow = PxtBelowTheta * R.t' / sum(PxtBelowTheta);
        
        % BUT these also need to factor in unabsorbed probability
        % (ie weighted average with max_dur). Intuitively, we know that we    
        % must scale pHB up for high bets, down for low bets. But I was
        % unable to calculate the adjusted pHBs and gave up, for now (see
        % old notes)
        % Instead, just set both to pHB, which is close enough when pHB>.98
        pHBhigh = pHB;
        pHBlow = pHB;

        meanRThigh_model(m,c,d,h) = pHBhigh*meanRThigh + (1-pHBhigh)*max_dur + Tnd;
        meanRTlow_model(m,c,d,h) = pHBlow*meanRTlow + (1-pHBlow)*max_dur + Tnd;
        
        % copy to trials
        RT_model_trialwise(Jdata) = meanRT_model(m,c,d,h);
        RThigh_model_trialwise(Jdata) = meanRThigh_model(m,c,d,h);
        RTlow_model_trialwise(Jdata) = meanRTlow_model(m,c,d,h);
        
        % grab the mean and sd(se?) for corresponding data
        I = Jdata & usetrs_data;
        meanRT_data(m,c,d,h) = mean(data.RT(I));
        sigmaRT_data(m,c,d,h) = std(data.RT(I))/sqrt(sum(I)); % SD or SEM?
        sigmaRT_data_trialwise(Jdata) = sigmaRT_data(c); % copy to trials, for alternate LL calculation below
                                                         % (should this be Jdata or I?)            
        % repeat for high/low
%         I = Jdata & usetrs_data & data.PDW_preAlpha==1; % preAlpha or not???
        I = Jdata & usetrs_data & data.PDW==1;
        meanRThigh_data(m,c,d,h) = mean(data.RT(I));
        sigmaRThigh_data(m,c,d,h) = std(data.RT(I))/sqrt(sum(I)); % SD or SEM?
        sigmaRThigh_data_trialwise(Jdata) = sigmaRThigh_data(h); % copy to trials, for alternate LL calculation below
%         I = Jdata & usetrs_data & data.PDW_preAlpha==0;
        I = Jdata & usetrs_data & data.PDW==0;
        meanRTlow_data(m,c,d,h) = mean(data.RT(I));
        sigmaRTlow_data(m,c,d,h) = std(data.RT(I))/sqrt(sum(I)); % SD or SEM?
        sigmaRTlow_data_trialwise(Jdata) = sigmaRTlow_data(h); % copy to trials, for alternate LL calculation below

    end

end

end
end
end

% copy trial params to data struct 'fit', then replace data with model vals
fit = data;
fit = rmfield(fit,'choice');
try fit = rmfield(fit,'correct'); end %#ok<TRYNC>
fit.pRight = pRight_model_trialwise;
% fit.pRightHigh = pRightHigh_model_trialwise;
% fit.pRightLow = pRightHigh_model_trialwise;

if options.conftask==1 % SEP
    fit.conf = conf_model_trialwise;
    fit.PDW = nan(size(fit.conf));
elseif options.conftask==2 % PDW
    fit.conf = pHigh_model_trialwise;
    fit.PDW = pHigh_model_trialwise; % legacy
end
if options.RTtask            
    fit.RT = RT_model_trialwise;
%     fit.RThigh = RThigh_model_trialwise;
%     fit.RTlow = RTlow_model_trialwise;
end

% also stored the 'parsed' values, for later plotting
parsedFit = struct();
parsedFit.pRight = pRight_model;
if options.conftask==1 % SEP
    parsedFit.confMean = meanConf_model;
elseif options.conftask==2 % PDW
    parsedFit.pHigh = pHigh_model;
    parsedFit.pRightHigh = pRightHigh_model;
    parsedFit.pRightLow = pRightLow_model;
    parsedFit.pHighCorr = pHighCorr_model;
    parsedFit.pHighErr = pHighErr_model;
end
if options.RTtask
    parsedFit.RTmean = meanRT_model;
    if options.conftask==2 % PDW
        parsedFit.RTmeanHigh = meanRThigh_model;
        parsedFit.RTmeanLow = meanRTlow_model;
    end
end


%% Next, calculate error (negative log-likelihood) under binomial/gaussian assumptions

% SJ 12-8-2022 this was missing before! Need choice to be 0..1
if max(data.choice)==2
    data.choice = data.choice-1;
end

% convert data vars to logicals
choice = logical(data.choice);
PDW    = logical(data.PDW);

% to avoid log(0) issues:
% minP = eps;
minP = 1e-300;
pRight_model_trialwise(pRight_model_trialwise==0) = minP; 
% % % pRight_model_trialwise(pRight_model_trialwise==1) = 1-minP;
pHigh_model_trialwise(pHigh_model_trialwise<=0) = minP; 
% % % pHigh_model_trialwise(pHigh_model_trialwise>=1) = 1-minP;


% CHOICE
% log likelihood of rightward choice on each trial, under binomial assumptions:
% LL_choice = sum(log(pRight_model_trialwise(choice))) + sum(log(1-pRight_model_trialwise(~choice)));

L_choice = binopdf(nRight_data,nTrials,pRight_model); 
L_choice(L_choice==0) = minP;
LL_choice = nansum(log(L_choice(:)));

% RT
LL_RT = 0;
if options.RTtask     
    
    % Option 1a:
    % likelihood of mean RTs for each coherence, under Gaussian
    % approximation, NOT separated by high/low bet
    L_RT = 1./(sigmaRT_data*sqrt(2*pi)) .* exp(-(meanRT_model-meanRT_data).^2 ./ (2*sigmaRT_data.^2));
%     L_RT = exp(-(meanRT_model-meanRT_data).^2 ./ (2*sigmaRT_data.^2));

    % Option 1b:
    % assign means to trials and sum over those, to keep same order of magnitude as other LLs
%     I = usetrs_data;
%     L_RT = 1./(sigmaRT_data_trialwise(I)*sqrt(2*pi)) .* exp(-(RT_model_trialwise(I) - data.RT(I)).^2 ./ (2*sigmaRT_data_trialwise(I).^2));
%     L_RT(L_RT==0) = min(L_RT(L_RT~=0)); % this seems weird; there shouldn't be any 0 vals anyway...

    L_RT(L_RT==0) = minP;
    LL_RT = nansum(log(L_RT(:))); % sum over all conditions (or trials)
    
    % Option 2a:
    % separate by high/low bet, fit mean RT for each coherence
%     L_RT_high = 1./(sigmaRThigh_data*sqrt(2*pi)) .* exp(-(meanRThigh_model-meanRThigh_data).^2 ./ (2*sigmaRThigh_data.^2));
%     L_RT_low = 1./(sigmaRTlow_data*sqrt(2*pi)) .* exp(-(meanRTlow_model-meanRTlow_data).^2 ./ (2*sigmaRTlow_data.^2));
%     LL_RT = log(max(minP,L_RT_high)) + log(max(minP,L_RT_low)); % this doesn't seem right; likelihoods are greater than 1

    % Option 2b:
    % separate by high/low bet and assign to trials, to keep order of mag consistent
%     % I = usetrs_data & data.PDW_preAlpha==1;
%     I = usetrs_data & data.PDW==1;
%     L_RT_high = 1./(sigmaRThigh_data_trialwise(I)*sqrt(2*pi)) .* exp(-(RThigh_model_trialwise(I) - data.RT(I)).^2 ./ (2*sigmaRThigh_data_trialwise(I).^2));    
%     L_RT_high(L_RT_high==0) = minP;
%     % I = usetrs_data & data.PDW_preAlpha==0;    
%     I = usetrs_data & data.PDW==0;    
%     L_RT_low = 1./(sigmaRTlow_data_trialwise(I)*sqrt(2*pi)) .* exp(-(RTlow_model_trialwise(I) - data.RT(I)).^2 ./ (2*sigmaRTlow_data_trialwise(I).^2));
%     L_RT_low(L_RT_low==0) = minP;
% 
%     LL_RT = nansum(log(L_RT_high(:))) + nansum(log(L_RT_low(:)));
%     LL_RT = nansum(log(L_RT_high(:))) / 385; % TEMP: ignore RTlow for now, and scale to same order of mag

    
    % Option 3: fit full RT distributions!
    
        
    
end

% CONF
LL_conf = 0;
if options.conftask==1 % SEP
    
    keyboard % this is unfinished
    % likelihood of mean conf ratings for each condition, under Gaussian approximation
    L_conf = 1./(sigmaConf_data*sqrt(2*pi)) .* exp(-(meanConf_model - meanConf_data).^2 ./ (2*sigmaConf_data.^2));
    L_conf(L_conf==0) = min(L_conf(L_conf~=0));
    LL_conf = nansum(log(L_conf(:))); % sum over all conditions

elseif options.conftask==2 % PDW

    % log likelihood of high bet on each trial, under binomial assumptions:
%     LL_conf = sum(log(pHigh_model_trialwise(PDW))) + sum(log(1-pHigh_model_trialwise(~PDW)));

    L_conf = binopdf(nHigh_data,nTrials,pHigh_model);
    L_conf(L_conf==0) = minP;
    LL_conf = nansum(log(L_conf(:)));

   
    % OR multinomial likelihood with n=1 and k=4 (categorial distribution)
%     RH = data.choice==1 & data.PDW==1;
%     RL = data.choice==1 & data.PDW==0;
%     LH = data.choice==0 & data.PDW==1;
%     LL = data.choice==0 & data.PDW==0;
%     LL_multinom = sum(log(max(pRightHigh_model_trialwise(RH),minP))) + sum(log(max(pRightLow_model_trialwise(RL),minP))) + ...
%                   sum(log(max(pLeftHigh_model_trialwise(LH),minP)))  + sum(log(max(pLeftLow_model_trialwise(LL),minP)));

    % reorganize for mnpdf input (needs k cols, one per outcome)
    X = cat(5,nRightHigh_data,nRightLow_data,nLeftLow_data,nLeftHigh_data);
    P = cat(5,pRightHigh_model,pRightLow_model,pLeftLow_model,pLeftHigh_model);

    X(1,2,:,:,:) = NaN; % vestibular 'high' coh
    P(1,2,:,:,:) = NaN; % not strictly necessary but for consistency

    X = reshape(X,[],4);
    P = reshape(P,[],4);

    L_multinom = mnpdf(X,P);
    L_multinom(L_multinom==0) = minP;
    LL_multinom = nansum(log(L_multinom(:)));
end

% total -LL
% err = -(LL_choice + LL_conf + LL_RT);
% OR
% err = -(LL_choice + LL_conf);
% OR
%err = -(LL_multinom + LL_RT);
% OR
err = -LL_multinom;


% OR: fit choice and RT first, then hold those fixed to find theta
% (Miguel's method)
% err = -(LL_choice + LL_RT);


%% print progress report!
if options.feedback
    % TODO, optionally list only non-fixed params
    fprintf('\n\n\n****************************************\n');
    fprintf('run %d\n', call_num);
%     fprintf('\tk= %g\n\tB= %g\n\tthetaVes= %g\n\tthetaVis= %g\n\tthetaComb= %g\n\talpha= %g\n\tTndVes= %g\n\tTndVis= %g\n\tTndComb= %g\n', kmult, B, theta, alpha, Tnds(1),Tnds(2),Tnds(3));
    fprintf('\tk= %g\n\tB= %g\n', kmult, B);
    fprintf('err: %f\n', err);
end
if options.feedback==2 && strcmp(options.fitMethod,'fms')
    figure(400); hold on;
    plot(call_num, err, '.','color','k','markersize',14);
    xlim([0 call_num+1])
    drawnow;
end
call_num = call_num + 1;
toc


%************
% temp
if options.dummyRun==0
% fprintf('\nerrChoice= %g\nerrRT= %g\nerrRTlo= %g\nerrRThi= %g\nerrConf= %g\nerrMult= %g\n', -LL_choice, -LL_RT, -nansum(log(L_RT_low(:))), -nansum(log(L_RT_high(:))), -LL_conf, -LL_multinom);
end
% LOW has smaller error than HIGH even though it looks way off by eye!
%************


end


% retrieve the full parameter set given the adjustable and fixed parameters 
function param2 = getParam ( param1 , guess , fixed )
  param2(fixed==0) = param1(fixed==0);  %get adjustable parameters from param1
  param2(fixed==1) = guess(fixed==1);   %get fixed parameters from guess
end
