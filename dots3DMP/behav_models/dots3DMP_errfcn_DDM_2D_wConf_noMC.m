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

kmult = param(1);
B = abs(param(2)); % don't accept negative bound heights
theta = abs(param(3:5)); % or negative thetas, one theta per mod
alpha = param(6); % base rate of low-conf choices
% % Tnd = param(7); % fixed Tnd, can't see any other way in this model
Tnds = param(7:9);
% if length(param)==7
%     TndR = param(7); 
%     TndL = param(7);
% else
%     TndR = param(7); % optionally, separate Tnds for L and R
%     TndL = param(8);
% end

cohs = unique(data.coherence);
mods = unique(data.modality);
hdgs = unique(data.heading);

deltas = 0; % fit only zero delta with unsigned MOI

% separate kvis and kves
kvis  = kmult*cohs'; % assume drift proportional to coh, reduces nParams
kves  = mean(kvis); % for now, assume 'straddling'

% ...but k for MOI for conf is constant
if all(mods==1), k = kves;
else, k = mean([kves kvis]);
end

timeToConf = 0; % placeholder
    % use method of images to calculate PDFs of DV and mapping to log odds corr
    R.t = 0.001:0.001:max_dur+timeToConf;
    R.Bup = B;
    R.drift = k * sind(hdgs(hdgs>=0)); % takes only unsigned drift rates
    R.lose_flag = 1; % we always need the losing densities
    R.plotflag = options.plot; % 1 = plot, 2 = plot nicer and export_fig (eg for talk)
    R.plotflag = 1; % manual override
    Pconf = images_dtb_2d(R); % /WolpertMOI

    % SJ need to check this works, and is in the right direction
    dt = R.t(2)-R.t(1);
    skipT = floor(timeToConf/dt);
    Pconf.t = circshift(Pconf.t,skipT);
    Pconf.logOddsCorrMap = circshift(Pconf.logOddsCorrMap,skipT,2);

if options.dummyRun==0
    logOddsMap = Pconf.logOddsCorrMap;
else
    logOddsMap = options.logOddsMap;
end




% now compute a separate P for model choices and RTs (modality-specific)
RVes.t = 0.001:0.001:max_dur;
RVes.Bup = B;
% *** marking differences in signed vs. unsigned ver ***
RVes.drift = kves * sind(hdgs(hdgs>=0));
RVes.lose_flag = 1;
RVes.plotflag = 0;
PVes =  images_dtb_2d(RVes);

for c = 1:length(cohs) % vis and comb
    RVis.t = 0.001:0.001:max_dur;
    RVis.Bup = B;
% *** marking differences in signed vs. unsigned ver ***
    RVis.drift = kvis(c) * sind(hdgs(hdgs>=0));
    RVis.lose_flag = 1;
    RVis.plotflag = 0;
    PVis(c) =  images_dtb_2d(RVis);
    
    kcomb(c) = sqrt(kves.^2 + kvis(c).^2); % optimal per Drugo
    RComb.t = 0.001:0.001:max_dur;
    RComb.Bup = B;
% *** marking differences in signed vs. unsigned ver ***
    RComb.drift = kcomb(c) * sind(hdgs(hdgs>=0));
    RComb.lose_flag = 1;
    RComb.plotflag = 0;
    PComb(c) =  images_dtb_2d(RComb);

    % compute w and mu to capture biases under cue conflict
    % not needed/possible here without signed heading, right?
    % so with unsigned hdg we can only fit delta=0
    % then we can use fit params to predict conflict data
%     wVes = sqrt( kves^2 / (kves^2 + kvis(c)^2) );
%     wVis = sqrt( kvis(c)^2 / (kves^2 + kvis(c)^2) );
%     for d = 1:length(deltas)
%         clear Rcomb
%         % positive delta defined as ves to the left, vis to the right
%         muVes = kves    * sind(hdgs-deltas(d)/2);
%         muVis = kvis(c) * sind(hdgs+deltas(d)/2);
%         
%         RComb.t = 0.001:0.001:duration;
%         RComb.Bup = B;
%         % *** marking differences in signed vs. unsigned ver ***
%         RComb.drift = wVes.*muVes + wVis.*muVis;
%         RComb.lose_flag = 1;
%         RComb.plotflag = 0;
%         PComb{c}{d} =  images_dtb_2d(RComb);
%     end
end


%% calculate likelihood of the observed data given the model parameters

% usetrs_data = data.correct | data.coherence<1e-6; % use only correct (OR ~0 hdg) trials for RT fits
usetrs_data = true(size(data.correct)); % try all trials

% (for parsedFit, and some likelihood calculations)
n = nan(length(mods),length(cohs),length(hdgs));
pRight_model        = n;
pHigh_model         = n;
pRightHigh_model    = n;
pRightLow_model     = n;
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

    if mods(m)==2, P = PVis(c);
    elseif mods(m)==3, P = PComb(c);
    end

for d = 1:length(deltas)
for h = 1:length(hdgs)

    % index for images_dtb (unsigned hdg/drift corresponding to this signed hdg)                    
    uh = abs(hdgs(h))==(hdgs(hdgs>=0));  % *** marking differences in signed vs. unsigned ver ***
      
    % grab the distributions from images_dtb:
    Pxt1 = squeeze(P.up.distr_loser(uh,:,:))'; % losing race pdf given correct
    Pxt2 = squeeze(P.lo.distr_loser(uh,:,:))'; % losing race pdf given incorrect

    % CHOICE & PDW
    if hdgs(h)<0 % leftward motion
        % if coh is negative, then pRight is p(incorrect), aka P.lo
        pRight = P.lo.p(uh)/(P.up.p(uh)+P.lo.p(uh)); % normalized by total bound-crossing probability 
        pLeft = 1-pRight; % (evidently, the unabsorbed probability can be ignored when it comes to pRight)

        % calculate probabilities of the four outcomes: right/left x high/low
        % these are the intersections, e.g. P(R n H)
        pRightHigh = sum(sum(Pxt2.*(logOddsMap>=theta(m))));
        pRightLow = sum(sum(Pxt2.*(logOddsMap<theta(m))));
        pLeftHigh = sum(sum(Pxt1.*(logOddsMap>=theta(m))));
        pLeftLow =  sum(sum(Pxt1.*(logOddsMap<theta(m))));                                                             
    
    else % rightward motion
        % if coh is positive, then pRight is p(correct), aka P.up
        pRight = P.up.p(uh)/(P.up.p(uh)+P.lo.p(uh)); % normalized by total bound-crossing probability
        pLeft = 1-pRight;

        % calculate probabilities of the four outcomes: right/left x high/low
        % these are the intersections, e.g. P(R n H)
        pRightHigh = sum(sum(Pxt1.*(logOddsMap>=theta(m))));
        pRightLow = sum(sum(Pxt1.*(logOddsMap<theta(m))));
        pLeftHigh = sum(sum(Pxt2.*(logOddsMap>=theta(m))));
        pLeftLow =  sum(sum(Pxt2.*(logOddsMap<theta(m))));                                                     
    end
    
%     if m==1 && hdgs(h)==0 && options.dummyRun, keyboard, end

    %ensure total prob sums to one (when it doesn't, it's because of unabsorbed prob)
    % THIS STEP WASN'T THOUGHT TO BE NECESSARY FOR PARAM RECOV, SO BE
    % CAREFUL AND REMOVE IT IF IT MUCKS IT UP
    Ptot = pRightHigh + pRightLow + pLeftHigh + pLeftLow;
    pRightHigh = pRightHigh/Ptot;
    pRightLow  = pRightLow/Ptot;
    pLeftHigh  = pLeftHigh/Ptot;
    pLeftLow   = pLeftLow/Ptot;

    % copy intersections to trials, for multinomial likelihood
    Jdata  = data.modality==mods(m) & data.coherence==cohs(c) & data.heading==hdgs(h) & data.delta==deltas(d);
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
    pRightHigh_model(m,c,d,h) = pRight_High; % does NOT include the effect of alpha
    pRightLow_model(m,c,d,h) = pRight_Low; % does NOT include the effect of alpha
    
% % %     pLeftHigh_model(m,c,d,h) = pLeftHigh; % MAYBE TEMP: SEE LIKELIHOOD
% % %     pLeftLow_model(m,c,d,h) = pLeftLow;  % MAYBE TEMP: SEE LIKELIHOOD
    
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
if options.conftask==1 % SEP
    fit.conf = conf_model_trialwise;
    fit.PDW = nan(size(fit.conf));
elseif options.conftask==2 % PDW
    fit.conf = pHigh_model_trialwise;
    fit.PDW = pHigh_model_trialwise; % legacy
end
if options.RTtask            
    fit.RT = RT_model_trialwise;
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

% convert data vars to logicals
choice = logical(data.choice);
PDW = logical(data.PDW);

% to avoid log(0) issues:
% minP = eps;
minP = 1e-300;
pRight_model_trialwise(pRight_model_trialwise==0) = minP; 
% % % pRight_model_trialwise(pRight_model_trialwise==1) = 1-minP;
pHigh_model_trialwise(pHigh_model_trialwise<=0) = minP; 
% % % pHigh_model_trialwise(pHigh_model_trialwise>=1) = 1-minP;


% CHOICE
% log likelihood of rightward choice on each trial, under binomial assumptions:
LL_choice = sum(log(pRight_model_trialwise(choice))) + sum(log(1-pRight_model_trialwise(~choice)));

% RT
LL_RT = 0;
if options.RTtask     
    
    % Option 1a:
    % likelihood of mean RTs for each coherence, under Gaussian
    % approximation, NOT separated by high/low bet
    L_RT = 1./(sigmaRT_data*sqrt(2*pi)) .* exp(-(meanRT_model-meanRT_data).^2 ./ (2*sigmaRT_data.^2));
    
    % Option 1b:
    % assign means to trials and sum over those, to keep same order of magnitude as other LLs
    I = usetrs_data;
    L_RT = 1./(sigmaRT_data_trialwise(I)*sqrt(2*pi)) .* exp(-(RT_model_trialwise(I) - data.RT(I)).^2 ./ (2*sigmaRT_data_trialwise(I).^2));
    L_RT(L_RT==0) = min(L_RT(L_RT~=0)); % this seems weird; there shouldn't be any 0 vals anyway...
    LL_RT = nansum(log(L_RT(:))); % sum over all conditions (or trials)
    
    % Option 2a:
    % separate by high/low bet, fit mean RT for each coherence
    L_RT_high = 1./(sigmaRThigh_data*sqrt(2*pi)) .* exp(-(meanRThigh_model-meanRThigh_data).^2 ./ (2*sigmaRThigh_data.^2));
    L_RT_low = 1./(sigmaRTlow_data*sqrt(2*pi)) .* exp(-(meanRTlow_model-meanRTlow_data).^2 ./ (2*sigmaRTlow_data.^2));
    LL_RT = log(max(minP,L_RT_high)) + log(max(minP,L_RT_low)); % this doesn't seem right; likelihoods are greater than 1

    % Option 2b:
    % separate by high/low bet and assign to trials, to keep order of mag consistent
%     I = usetrs_data & data.PDW_preAlpha==1;
    I = usetrs_data & data.PDW==1;
    L_RT_high = 1./(sigmaRThigh_data_trialwise(I)*sqrt(2*pi)) .* exp(-(RThigh_model_trialwise(I) - data.RT(I)).^2 ./ (2*sigmaRThigh_data_trialwise(I).^2));    
    L_RT_high(L_RT_high==0) = minP;
%     I = usetrs_data & data.PDW_preAlpha==0;    
    I = usetrs_data & data.PDW==0;    
    L_RT_low = 1./(sigmaRTlow_data_trialwise(I)*sqrt(2*pi)) .* exp(-(RTlow_model_trialwise(I) - data.RT(I)).^2 ./ (2*sigmaRTlow_data_trialwise(I).^2));
    L_RT_low(L_RT_low==0) = minP;

    
% % % %     LL_RT = nansum(log(L_RT_high(:))) + nansum(log(L_RT_low(:)));
% % %     LL_RT = nansum(log(L_RT_high(:))) / 385; % TEMP: ignore RTlow for now, and scale to same order of mag

    
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
    LL_conf = sum(log(pHigh_model_trialwise(PDW))) + sum(log(1-pHigh_model_trialwise(~PDW)));
    
    % OR multinomial likelihood with n=1 and k=4 (categorial distribution)
    RH = data.choice==1 & data.PDW==1;
    RL = data.choice==1 & data.PDW==0;
    LH = data.choice==0 & data.PDW==1;
    LL = data.choice==0 & data.PDW==0;
    LL_multinom = sum(log(max(pRightHigh_model_trialwise(RH),minP))) + sum(log(max(pRightLow_model_trialwise(RL),minP))) + ...
                  sum(log(max(pLeftHigh_model_trialwise(LH),minP)))  + sum(log(max(pLeftLow_model_trialwise(LL),minP)));
    
%     % sanity check: go coh by coh, should be the same
%     LL_choice_pdw2 = 0;
%     for c = 1:length(cohs)
%         I = data.scoh==cohs(c);
%         RH = I & data.choice==1 & data.PDW==1;
%         RL = I & data.choice==1 & data.PDW==0;
%         LH = I & data.choice==0 & data.PDW==1;
%         LL = I & data.choice==0 & data.PDW==0;
%         LL_c = sum(log(max(pRightHigh_model_trialwise(RH),minP))) + sum(log(max(pRightLow_model_trialwise(RL),minP))) + ...
%                sum(log(max(pLeftHigh_model_trialwise(LH),minP)))  + sum(log(max(pLeftLow_model_trialwise(LL),minP)));
%         LL_choice_pdw2 = LL_choice_pdw2 + LL_c;
%         % yep, close enough.
%            
% %         % is this the same as summing over k, using the total number of
% %         % occurrences of each outcome (as in "bag of words" blog post)?
% %         % [cannot include the factorial terms, xi is too large, turns inf]
% %         LL_choice_pdw3 = ...
% %         sum(RH)*log(max(pRightHigh_model(c),minP)) + ...
% %         sum(RL)*log(max(pRightLow_model(c),minP)) + ...
% %         sum(LH)*log(max(pLeftHigh_model(c),minP)) + ...
% %         sum(LL)*log(max(pLeftLow_model(c),minP));
% %           % definitely not. oh well.
%     
%     end

end

% total -LL
% err = -(LL_choice + LL_conf + LL_RT);
% OR
% err = -(LL_choice + LL_conf);
% OR
% err = -(LL_multinom + LL_RT);
% OR
err = -LL_multinom;


% OR: fit choice and RT first, then hold those fixed to find theta
% (Miguel's method)
% err = -(LL_choice + LL_RT);


%% print progress report!
if options.feedback
    fprintf('\n\n\n****************************************\n');
    fprintf('run %d\n', call_num);
    fprintf('\tk= %g\n\tB= %g\n\tthetaVes= %g\n\tthetaVis= %g\n\tthetaComb= %g\n\talpha= %g\n\tTndVes= %g\n\tTndVis= %g\n\tTndComb= %g\n', kmult, B, theta, alpha, Tnds(1),Tnds(2),Tnds(3));
    fprintf('err: %f\n', err);
end
if options.feedback==2 && strcmp(options.fitMethod,'fms')
    figure(400);
    plot(call_num, err, '.','MarkerSize',14);
    drawnow;
end
call_num = call_num + 1;
toc


%************
% temp
fprintf('\nerrChoice= %g\nerrRT= %g\nerrRTlo= %g\nerrRThi= %g\nerrConf= %g\nerrMult= %g\n', -LL_choice, -LL_RT, -nansum(log(L_RT_low(:))), -nansum(log(L_RT_high(:))), -LL_conf, -LL_multinom);
% LOW has smaller error than HIGH even though it looks way off by eye!
%************


end


% retrieve the full parameter set given the adjustable and fixed parameters 
function param2 = getParam ( param1 , guess , fixed )
  param2(fixed==0) = param1(fixed==0);  %get adjustable parameters from param1
  param2(fixed==1) = guess(fixed==1);   %get fixed parameters from guess
end
