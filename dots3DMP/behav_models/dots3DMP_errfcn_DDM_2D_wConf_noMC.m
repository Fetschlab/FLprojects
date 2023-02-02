function [err,fit,parsedFit] = dots3DMP_errfcn_DDM_2D_wConf_noMC(param, guess, fixed, data, options)

% updated model calculations to Steven's method:
% SJ 10-11-2021 no Monte Carlo simulation for fitting, it's redundant!
% just use model predictions directly

% uses UNSIGNED cohs for MOI calculation
% ONLY FIT zero delta

% SJ spawned 12-2022 from dots version, for dots3DMP

% SJ 02-2023 
% added R.lose_flag for all PDW related calculations - we don't need to
% compute losing distributions if not fitting conf

errFuncStart = tic;

% SJ 02-2023 change to persistent
persistent call_num
if isempty(call_num)
    call_num = 0;
end
call_num = call_num + 1;

persistent err_iter
if isempty(err_iter)
    err_iter = 0;
    
end



%% parse parameter inputs & relevant conditions, and set up for images_dtb

% retrieve the full parameter set given the adjustable and fixed parameters
param = getParam(param, guess, fixed);


paramNames = {'kmult','B',sprintf('%s_Ves',char(952)),sprintf('%s_Vis',char(952)),sprintf('%s_Comb',char(952)),'alpha','TndVes','TndVis','TndComb'};
kmult = param(1);
B     = abs(param(2)); % don't accept negative bound heights
theta = abs(param(3:5)); % or negative thetas, one theta per mod
alpha = param(6); % base rate of low-conf choices, same for all mods
Tnds  = param(7:9); % separate Tnd for each modality to account for RT offsets
% but perhaps this will become obsolete too with acc/vel?

% SJ 12-8-2022 LL calcs assume choice is 0..1
if max(data.choice)==2
    data.choice = data.choice-1;
end

cohs = unique(data.coherence);
mods = unique(data.modality);
hdgs = unique(data.heading);

% separate kvis and kves
kvis  = kmult*cohs'; % assume drift proportional to coh, reduces nParams
kves  = mean(kvis); % for now, assume 'straddling'

if options.dummyRun == 0
    deltas = 0;                  % fit only zero delta
else
    deltas = unique(data.delta); % predict all deltas
end


if isfield(data,'dur'), max_dur = max(data.dur)/1000; % max stimulus duration (s)
else,                   max_dur = 2.3;
end

% going to calculate separate confidence mapping for each modality since we
% are using temporal weighting
% can preset these R structs to be the same for all conditions
% SJ 01-2023 

R.t = 0.001:0.001:max_dur;
R.Bup = B;
R.lose_flag = any(contains(options.whichFit,{'conf','multinom'})); % only get losing distributions if we are fitting confidence
R.plotflag  = 0;

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


% SJ borrowed from Dots version 01-2023
% for fitting RT dists, need to  clip any data.RTs that may be past the limit in the mapping

dataRT_forPDF = data.RT;
if sum(dataRT_forPDF>max_dur)/length(dataRT_forPDF) > 0.05
    warning('Lots of long RTs, check the data or extend the mapping')
end
dataRT_forPDF(dataRT_forPDF>max_dur) = max_dur;


%% calculate likelihood of the observed data given the model parameters

% usetrs_data = data.correct | data.coherence<1e-6; % use only correct (OR ~0 hdg) trials for RT fits
usetrs_data = true(size(data.correct)); % try all trials

n = nan(length(mods),length(cohs),length(deltas),length(hdgs)); % (for parsedFit, and some likelihood calculations)
m = nan(length(data.heading),1); % trialwise vars

pRight_model  = n;
meanRT_model  = n;
meanRT_data   = n;

if R.lose_flag
    pHigh_model         = n;
   
    pHigh_Corr_model    = n;
    pHigh_Err_model     = n;

    meanConf_model      = n;
    meanRThigh_model    = n;
    meanRTlow_model     = n;

    pRight_High_model   = n;
    pRight_Low_model    = n;

    meanRThigh_data     = n;
    meanRTlow_data      = n;
    % sigmaRT_data        = n;
    % sigmaRThigh_data    = n;
    % sigmaRTlow_data     = n;

    sigmaRT_model     = n;
    sigmaRThigh_model = n;
    sigmaRTlow_model  = n;

    meanConf_data       = n;
    sigmaConf_data      = n;
end

% nCor  = n;
nTrials = n;


% marginals
pRight_model_trialwise   = m;
RT_model_trialwise       = m;
sigmaRT_model_trialwise  = m;
L_RT_fromPDF             = m;

if R.lose_flag
    pHigh_model_trialwise  = m;

    % intersections
    pRightHigh_model_trialwise = m;
    pRightLow_model_trialwise  = m;
    pLeftHigh_model_trialwise  = m;
    pLeftLow_model_trialwise   = m;

    % conditionals
    pRight_High_model_trialwise = m;
    pRight_Low_model_trialwise  = m;
    pHigh_Corr_model_trialwise  = m;
    pHigh_Err_model_trialwise   = m;

    % other
    conf_model_trialwise   = m;

    RThigh_model_trialwise = m;
    RTlow_model_trialwise  = m;
    sigmaRThigh_model_trialwise = m;
    sigmaRTlow_model_trialwise  = m;

    L_RT_PDW_fromPDF  = m;
end


for m = 1:length(mods)

    Tnd = Tnds(m);
    TndSD = 0;

    % for convolution with RT PDF later
    TndDist = normpdf(R.t,Tnd,2e-4)/sum(normpdf(R.t,Tnd,2e-4));

    for c = 1:length(cohs)

        for d = 1:length(deltas)

            if deltas(d)~=0 && mods(m)~=3, continue, end

            if mods(m)==1
                R.drift = sves .* kves .* sind(hdgs(hdgs>=0));
                R.driftSigned = kves * sind(hdgs); % sves is not needed here, since we only need this variable to determine overall sign of motion 
                                                    % (for predicting cue conflict effects, and appropriate use of p.up and p.lo)
            
            elseif mods(m)==2
                R.drift = svis .* kvis(c) .* sind(hdgs(hdgs>=0));
                R.driftSigned = kvis(c) * sind(hdgs);
            
            elseif mods(m)==3

                if options.dummyRun == 0
                    % 'kcomb'
                    R.driftSigned = sqrt(kves.^2 + kvis(c).^2) * sind(hdgs);
                    R.drift     = sqrt(sves.*kves.^2 + svis.*kvis(c).^2) .* sind(hdgs(hdgs>=0));

                else     % compute w and mu, and use signed headings, to capture biases under cue conflict

                    wVes = sqrt( kves^2 / (kves^2 + kvis(c)^2) );
                    wVis = sqrt( kvis(c)^2 / (kves^2 + kvis(c)^2) );

                    % positive delta defined as ves to the left, vis to the right
                    muVes = kves    * sind(hdgs-deltas(d)/2);
                    muVis = kvis(c) * sind(hdgs+deltas(d)/2);

                    R.driftSigned = wVes.*muVes + wVis.*muVis;
                    R.drift       = abs(sves.*wVes.*muVes + svis.*wVis.*muVis);

                end

            end
            
            P = images_dtb_2d_varDrift(R);

%             if options.dummyRun==0
                logOddsMap = P.logOddsCorrMap;
%             else
%                 % if not a dummyRun, we don't want to recalculate the
%                 maps with cue-conflict headings, we want ot use the same
%                 maps as before (either pass them in, or fix vars)
%                 for some reason, also an issue with interpolated
%                 headings, conf map is different even for same range
%             end


            for h = 1:length(hdgs)

                Jdata  = data.modality==mods(m) & data.coherence==cohs(c) & data.heading==hdgs(h) & data.delta==deltas(d);
                nTrials(m,c,d,h) = sum(Jdata);

                % SJ 01-2023 removing unnecessary matrices
                % if we do use them, make sure they are initialized
                % empirical probabilities from data (actually just take counts
                % 'successes', for pdf functions
                %     nRight_data(m,c,d,h) = sum(Jdata & data.choice==1); % choice is logical 0..1!!
                %
                %     if options.conftask==2
                %         nHigh_data(m,c,d,h)  = sum(Jdata & data.PDW==1);
                %
                %         nRightHigh_data(m,c,d,h) = sum(Jdata & data.choice==1 & data.PDW==1);
                %         nRightLow_data(m,c,d,h)  = sum(Jdata & data.choice==1 & data.PDW==0);
                %         nLeftHigh_data(m,c,d,h)  = sum(Jdata & data.choice==0 & data.PDW==1);
                %         nLeftLow_data(m,c,d,h)   = sum(Jdata & data.choice==0 & data.PDW==0);
                %     end

                % index for images_dtb (unsigned hdg/drift corresponding to this signed hdg)
                if options.dummyRun==1 && mods(m)==3
                    uh = h; % for cue conflict, drift rates are not symmetrical around 0, use the signed heading to get appropriate drift
                else
                    uh = abs(hdgs(h))==(hdgs(hdgs>=0));
                end

                if R.lose_flag
                    % grab the distributions from images_dtb:
                    Pxt1 = squeeze(P.up.distr_loser(uh,:,:))'; % losing race pdf given correct
                    Pxt2 = squeeze(P.lo.distr_loser(uh,:,:))'; % losing race pdf given incorrect
                end

                % CHOICE & PDW
                
                % SJ: yes, we can fit non-zero bias, just need to use the bias-corrected, signed
                % drift rate, rather than the experimenter-defined headings!
                if P.driftSigned(h)<0 % leftward motion

                    % if hdg is negative, then pRight is p(incorrect), aka P.lo
                    pRight = P.lo.p(uh)/(P.up.p(uh)+P.lo.p(uh)); % normalized by total bound-crossing probability
                    pLeft = 1-pRight; % (evidently, the unabsorbed probability can be ignored when it comes to pRight)

                    if R.lose_flag
                        % calculate probabilities of the four outcomes: right/left x high/low
                        % these are the intersections, e.g. P(R n H)
                        pRightHigh = sum(sum(Pxt2.*(P.logOddsCorrMap >= theta(m))));
                        pRightLow  = sum(sum(Pxt2.*(P.logOddsCorrMap < theta(m))));
                        pLeftHigh  = sum(sum(Pxt1.*(P.logOddsCorrMap >= theta(m))));
                        pLeftLow   = sum(sum(Pxt1.*(P.logOddsCorrMap < theta(m))));
                    end

                else % rightward motion (and P is symmetrical at 0 so it's fine that that gets taken care of here too)

                    % if hdg is positive, then pRight is p(correct), aka P.up
                    pRight = P.up.p(uh)/(P.up.p(uh)+P.lo.p(uh)); % normalized by total bound-crossing probability
                    pLeft = 1-pRight;

                    if R.lose_flag
                        % calculate probabilities of the four outcomes: right/left x high/low
                        % these are the intersections, e.g. P(R n H)
                        pRightHigh = sum(sum(Pxt1.*(logOddsMap>=theta(m))));
                        pRightLow  = sum(sum(Pxt1.*(logOddsMap<theta(m))));
                        pLeftHigh  = sum(sum(Pxt2.*(logOddsMap>=theta(m))));
                        pLeftLow   =  sum(sum(Pxt2.*(logOddsMap<theta(m))));
                    end
                end


                if R.lose_flag % && options.conftask==2 | 'fitConditionals'

                    % SJ 02-2023 reducing lines here and number of scalar vars, all calcs can be run on array
                    % actually  doing it this way, we could do a
                    % similar multinomial for SEP task (with >2 columns for
                    % binned confidence reports)

                    pChoiceAndWager = [pRightHigh, pRightLow; pLeftHigh, pLeftLow];
                    pChoice = [pRight;pLeft];

                    [pWager_Choice,pChoice_Wager,pWager,pChoiceAndWager] = int_to_margconds(pChoiceAndWager,pChoice);
                    
                    % adjust marginal for the base rate of low-conf bets (alpha)
                    pWager_wAlpha = pWager + [-1 1]*alpha*pWager(1);
                    
                    % and adjust conditionals and marginals accordingly
                    [pWager_Choice_wAlpha, pChoiceAndWager_wAlpha] = margconds_to_int(pChoice_Wager,pChoice,pWager_wAlpha);

                    %{
                    %ensure total prob sums to one (when it doesn't, it's because of unabsorbed prob)
                    % THIS STEP WASN'T THOUGHT TO BE NECESSARY FOR PARAM RECOV, SO BE
                    % CAREFUL AND REMOVE IT IF IT MUCKS IT UP
                    Ptot = pRightHigh + pRightLow + pLeftHigh + pLeftLow;
                    %if abs(Ptot-1)>1e-6; warning('Ptot'); keyboard, fprintf('m=%d, c=%g, h=%g\tPtot=%.2f\n',mods(m),cohs(c),hdgs(h),Ptot); end
                    pRightHigh = pRightHigh/Ptot;
                    pRightLow  = pRightLow/Ptot;
                    pLeftHigh  = pLeftHigh/Ptot;
                    pLeftLow   = pLeftLow/Ptot;

                    %     % these are intersections, for model fitting
                    %     pRightHigh_model(m,c,d,h) = pRightHigh;
                    %     pRightLow_model(m,c,d,h) = pRightLow;
                    %     pLeftHigh_model(m,c,d,h) = pLeftHigh;
                    %     pLeftLow_model(m,c,d,h) = pLeftLow;
                    %
                    %     % copy intersections to trials, for multinomial likelihood
                    %     pRightHigh_model_trialwise(Jdata) = pRightHigh;
                    %     pRightLow_model_trialwise(Jdata) = pRightLow;
                    %     pLeftHigh_model_trialwise(Jdata) = pLeftHigh;
                    %     pLeftLow_model_trialwise(Jdata) = pLeftLow;
                    % this step now occurs below, after incorporating alpha

                    % by definition of conditional probability: P(A|B) = P(A n B) / P(B)
                    pHigh_Right = pRightHigh / pRight;
                    pLow_Right  = pRightLow / pRight;
                    pHigh_Left  = pLeftHigh / pLeft;
                    pLow_Left   = pLeftLow / pLeft;

                    % by law of total probability
                    pHigh = pHigh_Right*pRight + pHigh_Left*pLeft;
                    pLow  = pLow_Right*pRight + pLow_Left*pLeft; % why isn't this 1-pHigh? probability is so weird.

                    % by Bayes' rule: choive GIVEN wager
                    pRight_High = pHigh_Right * pRight / pHigh;
                    pRight_Low  = pLow_Right * pRight / pLow;
                    pLeft_High  = pHigh_Left * pLeft / pHigh;
                    pLeft_Low   = pLow_Left * pLeft / pLow;

                    % adjust the probabilities for the base rate of low-conf bets (alpha)
                    pHigh_wAlpha = pHigh - alpha*pHigh;

                    % Lastly, the PDW split; use Bayes to adjust P(H|R) and P(H|L)
                    pHigh_Right_wAlpha = pRight_High * pHigh_wAlpha / pRight;
                    pHigh_Left_wAlpha  = pLeft_High * pHigh_wAlpha / pLeft;
                    if hdgs(h)<0 % leftward motion
                        pHigh_Corr = pHigh_Left_wAlpha;
                        pHigh_Err  = pHigh_Right_wAlpha;
                    else % rightward motion
                        pHigh_Corr = pHigh_Right_wAlpha;
                        pHigh_Err  = pHigh_Left_wAlpha;
                    end

                    % and use (rearranged) definition of conditional probability to adjust
                    % the intersections
                    pRightHigh_wAlpha = pHigh_wAlpha * pRight_High;
                    pRightLow_wAlpha  = (1-pHigh_wAlpha) * pRight_Low;
                    pLeftHigh_wAlpha  = pHigh_wAlpha * pLeft_High;
                    pLeftLow_wAlpha   = (1-pHigh_wAlpha) * pLeft_Low;

                    % still should sum to 1
                    Ptot = pRightHigh_wAlpha + pRightLow_wAlpha + pLeftHigh_wAlpha + pLeftLow_wAlpha;
                    if abs(Ptot-1)>1e-6; warning('Ptot'); fprintf('m=%d, c=%g, h=%g\tPtot=%.2f\n',mods(m),cohs(c),hdgs(h),Ptot); end
                    pRightHigh_wAlpha = pRightHigh_wAlpha/Ptot;
                    pRightLow_wAlpha  = pRightLow_wAlpha/Ptot;
                    pLeftHigh_wAlpha  = pLeftHigh_wAlpha/Ptot;
                    pLeftLow_wAlpha   = pLeftLow_wAlpha/Ptot;

                    %}
                                      
                    % pWager_Choice_wAlpha: [Hi_R,Lo_R;Hi_L,Lo_L]
                    % pWager_Acc_wAlpha:    [Hi_Corr,Lo_Corr;Hi_Err,Lo_Err]
                    
                    if hdgs(h)<0 % leftward motion, flipud Choice matrix
                        pWager_Acc_wAlpha = flipud(pWager_Choice_wAlpha);
                    else % rightward motion (or zero)
                        pWager_Acc_wAlpha = pWager_Choice_wAlpha;
                    end                    

                    % % could generalize below assignments to multi-split
                    % confidence, but leaving for now
                    % e.g. pChoiceWager_model_trialwise(Jdata,:) = pChoiceAndWager_wAlpha;

                    % copy intersections to trials, for multinomial likelihood
                    % *NOW INCLUDES effect of alpha
                    % rows of pChoiceAndWager are R/L, columns are Hi/Lo
                    
                    pRightHigh_model_trialwise(Jdata) = pChoiceAndWager_wAlpha(1,1);
                    pRightLow_model_trialwise(Jdata)  = pChoiceAndWager_wAlpha(1,2);
                    pLeftHigh_model_trialwise(Jdata)  = pChoiceAndWager_wAlpha(2,1);
                    pLeftLow_model_trialwise(Jdata)   = pChoiceAndWager_wAlpha(2,2);

                    % copy to arrays for parsedFit struct
                    pHigh_model(m,c,d,h)       = pWager_wAlpha(1);   % includes the effect of alpha
                    pRight_High_model(m,c,d,h) = pChoice_Wager(1,1); % does NOT include the effect of alpha, since bet high/low is already the given
                    pRight_Low_model(m,c,d,h)  = pChoice_Wager(1,2); % does NOT include the effect of alpha

                    pHigh_Corr_model(m,c,d,h)  = pWager_Acc_wAlpha(1,1);   % includes the effect of alpha
                    pHigh_Err_model(m,c,d,h)   = pWager_Acc_wAlpha(2,1);   % includes the effect of alpha

                    % copy to trials for fit struct
                    pHigh_model_trialwise(Jdata)       = pHigh_model(m,c,d,h);
                    pRight_High_model_trialwise(Jdata) = pRight_High_model(m,c,d,h);
                    pRight_Low_model_trialwise(Jdata)  = pRight_Low_model(m,c,d,h);
                    pHigh_Corr_model_trialwise(Jdata)  = pHigh_Corr_model(m,c,d,h);
                    pHigh_Err_model_trialwise(Jdata)   = pHigh_Err_model(m,c,d,h);
                end

                pRight_model(m,c,d,h)         = pRight;
                pRight_model_trialwise(Jdata) = pRight_model(m,c,d,h);

                % ************
                % RT
                % ************
                %     nCor(c) = sum(Jdata & (data.correct | abs(data.heading)<1e-6));
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


                    %         meanRT_model(m,c,d,h) = P.up.p(uh)*P.up.mean_t(uh) + P.lo.p(uh)*P.lo.mean_t(uh) + (1-pHB)*max_dur + Tnd; % weighted average
                    meanRT_model(m,c,d,h) = pHB*P.up.mean_t(uh) + (1-pHB)*max_dur + Tnd; % weighted average
                    %         varRT_model(m,c,d,h)  = pHB*P.up.var_t(uh)  + (1-pHB)*0 + TndSD.^2; % TndSD?

                    RT_model_trialwise(Jdata) = meanRT_model(m,c,d,h);

                    if R.lose_flag %&& options.conftask==2
                        % for RT conditioned on wager, it seems we don't need to do any
                        % scaling by Ptb; the correct distribution is captured by the
                        % conditional Pxts, we just sum them, dot product w t and normalize
                        PxtAboveTheta = sum(Pxt1.*(logOddsMap>=theta(m))); % shouldn't matter if Pxt1 or 2, it's symmetric
                        PxtBelowTheta = sum(Pxt1.*(logOddsMap<theta(m)));

%                         PxtAboveTheta2 = sum(Pxt2.*(logOddsMap>=theta(m))); % shouldn't matter if Pxt1 or 2, it's symmetric
%                         PxtBelowTheta2 = sum(Pxt2.*(logOddsMap<theta(m)));

                        meanRThigh = PxtAboveTheta * R.t' / sum(PxtAboveTheta);
                        meanRTlow  = PxtBelowTheta * R.t' / sum(PxtBelowTheta);

%                         meanRThigh = ( (PxtAboveTheta * P.up.p(uh) / sum(PxtAboveTheta)) + (PxtAboveTheta2 * P.lo.p(uh) / sum(PxtAboveTheta2))) * R.t';
%                         meanRTlow = ( (PxtBelowTheta * P.up.p(uh) / sum(PxtBelowTheta)) + (PxtBelowTheta2 * P.lo.p(uh) / sum(PxtBelowTheta2))) * R.t';

                        % BUT these also need to factor in unabsorbed probability
                        % (ie weighted average with max_dur). Intuitively, we know that we
                        % must scale pHB up for high bets, down for low bets. But I was
                        % unable to calculate the adjusted pHBs and gave up, for now (see
                        % old notes)
                        % Instead, just set both to pHB, which is close enough when pHB>.98
                        pHBhigh = pHB;
                        pHBlow  = pHB;

                        meanRThigh_model(m,c,d,h) = pHBhigh*meanRThigh + (1-pHBhigh)*max_dur + Tnd;
                        meanRTlow_model(m,c,d,h)  = pHBlow*meanRTlow + (1-pHBlow)*max_dur + Tnd;

                        % copy to trials
                        RThigh_model_trialwise(Jdata) = meanRThigh_model(m,c,d,h);
                        RTlow_model_trialwise(Jdata)  = meanRTlow_model(m,c,d,h);

                    end


                    %{
                    % TEMP(?): compute the variances of these PDFs for the
                    % fitting-by-mean RT likelihood calculations below. This may become
                    % obsolete if we can use the PDFs themselves

                    % variance = sum((x - mean).^2 .* y)
                    % where x is the domain of the PDF, y is the PDF values over the
                    % domain x, and mean is the mean of the PDF, which you can
                    % calculate using the formula you provided: mean = sum(y.*x). -- thanks ChatGPT!

                    % % %         % SJ 12-2022 Is diag is what we want here? Var(X) = E((X-E(X))^2)
                    % % %         P.up.var_t = diag(P.up.pdf_t * ((R.t - P.up.mean_t).^2)' ./ P.up.p);
                    % % %         P.lo.var_t = diag(P.lo.pdf_t * ((R.t - P.lo.mean_t).^2)' ./ P.lo.p);

                    % ___note the PDFs are normalized___
                    sigmaRT_model(m,c,d,h) = sqrt(sum((R.t - P.up.mean_t(uh)).^2 .* P.up.pdf_t(uh,:)/sum(P.up.pdf_t(uh,:))));
                    
                    if R.lose_flag
                        sigmaRThigh_model(m,c,d,h) = sqrt(sum((R.t - meanRThigh).^2 .* PxtAboveTheta/sum(PxtAboveTheta)));
                        sigmaRTlow_model(m,c,d,h) = sqrt(sum((R.t - meanRTlow).^2 .* PxtBelowTheta/sum(PxtBelowTheta)));                  
                    end

                    % (note these do not include any variance added from convolution w Tnd dist)
                    sigmaRT_model_trialwise(Jdata) = sigmaRT_model(m,c,d,h);

                    if R.lose_flag
                        sigmaRThigh_model_trialwise(Jdata) = sigmaRThigh_model(m,c,d,h);
                        sigmaRTlow_model_trialwise(Jdata) = sigmaRTlow_model(m,c,d,h);
                    end

                    % grab the mean and sd(se?) for corresponding data
                    %I = Jdata & usetrs_data;
                    %n_RT_data(m,c,d,h) = sum(I);
                    %meanRT_data(m,c,d,h) = mean(data.RT(I));
                    %sigmaRT_data(m,c,d,h) = std(data.RT(I))/sqrt(sum(I)); % SD or SEM? --*obsolete even if not using full dists, use Sigma from the dist itself
                    %sigmaRT_data_trialwise(Jdata) = sigmaRT_data(c); % copy to trials, for alternate LL calculation below
                    % (should this be Jdata or I?)
                    %}

                    % but we can just get the likelihood directly from the PDF!
                    PDF = P.up.pdf_t(uh,:)/sum(P.up.pdf_t(uh,:)); % renormalize % WHY ONLY P.UP HERE??
%                     PDF2 = P.lo.pdf_t(uh,:)/sum(P.lo.pdf_t(uh,:)); % renormalize


                    PDFwTnd = conv(PDF,TndDist); % convolve with ^Tnd dist (currently zero variance, approx something tiny)
                    PDFwTnd = PDFwTnd(1:length(P.t))/sum(PDFwTnd(1:length(P.t))); % normalize again and trim the extra length from conv
                    % figure;plot(PDF,'b');hold on; plot(PDFwTnd,'r'); % temp, just to see that the conv makes sense
                    L_RT_fromPDF(I) = PDFwTnd(round(dataRT_forPDF(I)*1000))';

                    % repeat for high/low bet separately
                    if R.lose_flag %&& options.conftask==2
                        %         I = Jdata & usetrs_data & data.PDW_preAlpha==1;
                        %I = Jdata & usetrs_data & data.PDW==1;
                        %n_RThigh_data(m,c,d,h) = sum(I);
                        %meanRThigh_data(m,c,d,h) = mean(data.RT(I));
                        %         sigmaRThigh_data(m,c,d,h) = std(data.RT(I))/sqrt(sum(I)); % SD or SEM?
                        %         sigmaRThigh_data_trialwise(Jdata) = sigmaRThigh_data(h); % copy to trials, for alternate LL calculation below

                        % but we can just get the likelihood directly from the PDF!
                        PDF = PxtAboveTheta/sum(PxtAboveTheta);
                        PDFwTnd = conv(PDF,TndDist);
                        PDFwTnd = PDFwTnd(1:length(P.t))/sum(PDFwTnd(1:length(P.t)));
                        L_RT_PDW_fromPDF(I) = PDFwTnd(round(dataRT_forPDF(I)*1000))';
                        if any(isnan(L_RT_PDW_fromPDF(I))); keyboard; end % TEMP

                        %         I = Jdata & usetrs_data & data.PDW_preAlpha==0;
                        %I = Jdata & usetrs_data & data.PDW==0;
                        %n_RTlow_data((m,c,d,h)) = sum(I);
                        %meanRTlow_data(m,c,d,h) = mean(data.RT(I));
                        %         sigmaRTlow_data(m,c,d,h) = std(data.RT(I))/sqrt(sum(I)); % SD or SEM?
                        %         sigmaRTlow_data_trialwise(Jdata) = sigmaRTlow_data(h); % copy to trials, for alternate LL calculation below

                        % but we can just get the likelihood directly from the PDF!
                        PDF = PxtBelowTheta/sum(PxtBelowTheta);
                        PDFwTnd = conv(PDF,TndDist);
                        PDFwTnd = PDFwTnd(1:length(P.t))/sum(PDFwTnd(1:length(P.t)));
                        L_RT_PDW_fromPDF(I) = PDFwTnd(round(dataRT_forPDF(I)*1000))';
                        if any(isnan(L_RT_PDW_fromPDF(I))); keyboard; end % TEMP

                    end
                end
            end
        end
    end
end

% copy trial params to data struct 'fit',
% replacing observed data with model predictions

fit = data;
fit = rmfield(fit,'choice');
try fit = rmfield(fit,'correct'); end %#ok<TRYNC>
fit.pRight = pRight_model_trialwise;

if R.lose_flag
    fit.pRightHigh = pRightHigh_model_trialwise;
    fit.pRightLow = pRightHigh_model_trialwise;

    if options.conftask==1 % SEP
        fit.conf = conf_model_trialwise;
        fit.PDW  = nan(size(fit.conf));
    elseif options.conftask==2 % PDW
        fit.conf = pHigh_model_trialwise;
        fit.PDW  = pHigh_model_trialwise; % legacy
    end
    if options.RTtask
        fit.RT = RT_model_trialwise;
        fit.RThigh = RThigh_model_trialwise;
        fit.RTlow = RTlow_model_trialwise;
    end
end

% also stored the 'parsed' values, for later plotting
parsedFit = struct();
parsedFit.pRight = pRight_model;

if R.lose_flag
    if options.conftask==1 % SEP
        parsedFit.confMean = meanConf_model;
    elseif options.conftask==2 % PDW
        parsedFit.pHigh = pHigh_model;
        parsedFit.pRightHigh = pRight_High_model;
        parsedFit.pRightLow = pRight_Low_model;
        parsedFit.pHighCorr = pHigh_Corr_model;
        parsedFit.pHighErr = pHigh_Err_model;
    end
end
if options.RTtask
    parsedFit.RTmean = meanRT_model;
    if options.conftask==2 && R.lose_flag % PDW
        parsedFit.RTmeanHigh = meanRThigh_model;
        parsedFit.RTmeanLow = meanRTlow_model;
    end
end


%% Calculate error (negative log-likelihood) 
% choice and PDW indep binomials, or one multinomial
% RT originally fit means under gaussian assumptions. but now try fit PDF

% convert data vars to logicals
choice = logical(data.choice);
PDW    = logical(data.PDW);

% to avoid log(0) issues:
minP = 1e-300;
pRight_model_trialwise(pRight_model_trialwise==0) = minP;

% CHOICE
% log likelihood of rightward choice on each trial, under binomial assumptions:
LL_choice = sum(log(pRight_model_trialwise(choice))) + sum(log(1-pRight_model_trialwise(~choice)));

% binomial over conditions, rather than bernoulli over trials
% SJ 12-2022, I guess this is equivalent to above bernoulli
% L_choice = binopdf(nRight_data,nTrials,pRight_model);
% L_choice(L_choice==0) = minP; % avoid log(0)
% LL_choice = nansum(log(L_choice(:)));

% CHOICE, from conditionals, currently unused
%{
pRight_Low_model_trialwise(pRight_Low_model_trialwise==0) = minP;
pp = pRight_Low_model_trialwise(data.PDW==0); % conditionalize first
cc = choice(data.PDW==0);
LL_choice_low = sum(log(pp(cc))) + sum(log(1-pp(~cc)));

pRight_High_model_trialwise(pRight_High_model_trialwise==0) = minP;
pp = pRight_High_model_trialwise(data.PDW==1); % conditionalize first
cc = choice(data.PDW==1);
LL_choice_high = sum(log(pp(cc))) + sum(log(1-pp(~cc)));
LL_choice_cond = LL_choice_low + LL_choice_high; % currently unused
%}

% RT
LL_RT = 0;

% RT
if options.RTtask

    %     sigmaRT_data = sigmaRT_data.^2 % ***this seems more comparable to the theoretical SEM from Palmer, keep it as an option.

    %     L_RToption = 'means_nosep';
    %     L_RToption = 'trials_nosep';
    %     L_RToption = 'means_sep';
    %     L_RToption = 'trials_sep';
    L_RToption = 'full_dist';

    switch L_RToption
        case 'means_nosep'

            % Option 1a: means, not sep
            % (fit mean RTs for each coh under Gaussian approximation, NOT separated by high/low bet)
            L_RT = 1./(sigmaRT_model*sqrt(2*pi)) .* exp(-(meanRT_model-meanRT_data).^2 ./ (2*sigmaRT_model.^2)) / 1000; % div by 1000 to account for units (sec/ms)
            % *** NOT SURE ABOUT THE FACTOR OF 1000 ***
            L_RT(L_RT==0) = minP; % kluge to remove zeros, which are matlab's response to  exp(-[too large a number]), otherwise the log will give you infinity
            LL_RT = sum(log(L_RT(:))); % sum over all conditions (or trials)
            %             LL_RT_high = NaN; LL_RT_low = NaN;

        case 'trials_nosep'

            % Option 1b: trials, not sep
            % (assign means to trials and sum over those, to keep same order of mag as other LLs)
            I = usetrs_data;
            L_RT = 1./(sigmaRT_model_trialwise(I)*sqrt(2*pi)) .* exp(-(RT_model_trialwise(I) - data.RT(I)).^2 ./ (2*sigmaRT_model_trialwise(I).^2)) / 1000;
            L_RT(L_RT==0) = minP;
            LL_RT = sum(log(L_RT(:)));
            %             LL_RT_high = NaN; LL_RT_low = NaN;

        case 'means_sep'

            % Option 2a: means, sep hi/lo
            % (fit mean RT for each coherence, separately for high/low bet)
            L_RT_high = 1./(sigmaRThigh_model*sqrt(2*pi)) .* exp(-(meanRThigh_model-meanRThigh_data).^2 ./ (2*sigmaRThigh_model.^2)) / 1000;
            L_RT_high(L_RT_high==0) = minP;
            % WEIGHT LL by nTrials, e.g. since there are fewer low bets
            % at high coh, a miss there should penalize less
            weightedLL = n_RThigh_data/sum(n_RThigh_data) .* log(L_RT_high) * length(cohs);
            LL_RT_high = sum(weightedLL);

            L_RT_low = 1./(sigmaRTlow_model*sqrt(2*pi)) .* exp(-(meanRTlow_model-meanRTlow_data).^2 ./ (2*sigmaRTlow_model.^2)) / 1000;
            L_RT_low(L_RT_low==0) = minP;
            weightedLL = n_RTlow_data/sum(n_RTlow_data) .* log(L_RT_low) * length(cohs);
            LL_RT_low = sum(weightedLL);

            LL_RT = LL_RT_high + LL_RT_low;

        case 'trials_sep'

            % Option 2b: trials, sep hi/lo
            % (separate by high/low bet and assign to trials, to keep order of mag)
            %     I = usetrs_data & data.PDW_preAlpha==1; %??
            I = usetrs_data & data.PDW==1;
            L_RT_high = 1./(sigmaRThigh_model_trialwise(I)*sqrt(2*pi)) .* exp(-(RThigh_model_trialwise(I) - data.RT(I)).^2 ./ (2*sigmaRThigh_model_trialwise(I).^2)) / 1000;
            L_RT_high(L_RT_high==0) = minP;
            LL_RT_high = sum(log(L_RT_high));

            %     I = usetrs_data & data.PDW_preAlpha==0; %??
            I = usetrs_data & data.PDW==0;
            L_RT_low = 1./(sigmaRTlow_model_trialwise(I)*sqrt(2*pi)) .* exp(-(RTlow_model_trialwise(I) - data.RT(I)).^2 ./ (2*sigmaRTlow_model_trialwise(I).^2)) / 1000;
            L_RT_low(L_RT_low==0) = minP;
            LL_RT_low = sum(log(L_RT_low));

            LL_RT = LL_RT_high + LL_RT_low;

        case 'full_dist'

            % Option 3: fit full RT distributions
            L_RT_fromPDF(L_RT_fromPDF==0) = minP;
            LL_RT = sum(log(L_RT_fromPDF)); % no need to use uncond

            if R.lose_flag
%             L_RT_PDW_fromPDF(L_RT_PDW_fromPDF==0) = minP;
%             LL_RT = nansum(log(L_RT_PDW_fromPDF));
            end

    end

end

% CONF
if R.lose_flag

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

        % CONF, from conditionals
        %     pHigh_Corr_model_trialwise(pHigh_Corr_model_trialwise==0) = minP;
        %     pp = pHigh_Corr_model_trialwise(data.correct==1); % conditionalize first
        %     cc = PDW(data.correct==1);
        %     LL_conf_corr = sum(log(pp(cc))) + sum(log(1-pp(~cc)));
        %
        %     pHigh_Err_model_trialwise(pHigh_Err_model_trialwise==0) = minP;
        %     pp = pHigh_Err_model_trialwise(data.correct==0); % conditionalize first
        %     cc = PDW(data.correct==0);
        %     LL_conf_err = sum(log(pp(cc))) + sum(log(1-pp(~cc)));
        %
        %     LL_conf_cond = LL_conf_corr + LL_conf_err; % currently unused


        % CONF, binomial method
        %     L_conf = binopdf(nHigh_data,nTrials,pHigh_model);
        %     L_conf(L_conf==0) = minP;
        %     LL_conf = nansum(log(L_conf(:)));


        % OR multinomial likelihood with n=1 and k=4 (categorial distribution)
        RH = data.choice==1 & data.PDW==1;
        RL = data.choice==1 & data.PDW==0;
        LH = data.choice==0 & data.PDW==1;
        LL = data.choice==0 & data.PDW==0;
        LL_multinom = sum(log(max(pRightHigh_model_trialwise(RH),minP))) + sum(log(max(pRightLow_model_trialwise(RL),minP))) + ...
            sum(log(max(pLeftHigh_model_trialwise(LH),minP)))  + sum(log(max(pLeftLow_model_trialwise(LL),minP)));

        % CONF multinomial, using mnpdf
        % reorganize for mnpdf input (needs k cols, one per outcome)
        %     X = cat(5,nRightHigh_data,nRightLow_data,nLeftLow_data,nLeftHigh_data);
        %     P = cat(5,pRightHigh_model,pRightLow_model,pLeftLow_model,pLeftHigh_model);
        %
        %     X(1,2,:,:,:) = NaN; % vestibular 'high' coh
        %     P(1,2,:,:,:) = NaN; % not strictly necessary but for consistency
        %
        %     X = reshape(X,[],4);
        %     P = reshape(P,[],4);
        %
        %     L_multinom = mnpdf(X,P);
        %     L_multinom(L_multinom==0) = minP;
        %     LL_multinom = nansum(log(L_multinom(:)));

    end
end

% total -LL
for f = 1:length(options.whichFit)
    err = err - eval(['LL_',options.whichFit{1}]);
end
err_iter(end+1) = err;

%% print progress report!
if options.feedback
    fprintf('\n\n\n****************************************\n');
    fprintf('run %d\n', call_num);
    %     fprintf('\tkmult= %g\n\tB= %g\n\tthetaVes= %g\n\tthetaVis= %g\n\tthetaComb= %g\n\talpha= %g\n\tTndVes= %g\n\tTndVis= %g\n\tTndComb= %g\n', kmult, B, theta, alpha, Tnds(1),Tnds(2),Tnds(3));

    for p = 1:length(param)
        if ~fixed(p) % only print fitted parameter values
            fprintf('\t%s\t= %g\n',paramNames{p},param(p));
        end
    end

    fprintf('err: %f\n', err);

end
if options.feedback==2 && strcmp(options.fitMethod,'fms')
    if call_num == 1
        figure(400); h = plot(call_num,err_iter,'.','color','k','markersize',14);
        h.Parent.XLabel.String = 'Iteration';
        h.Parent.YLabel.String = 'Total Error';
        set(h,'YDataSource','err_all');
        set(h,'XDataSource','call_num');
    else
        set(h,'XData',call_num,'YData',err_iter); drawnow;
    end
%     figure(400); hold on;
%     plot(call_num, err, '.','color','k','markersize',14);
%     xlim([0 call_num])
%     drawnow;
end
errFuncElapsed = toc(errFuncStart)


%************
% temp, diagnosing LL issues
% fprintf('\nLL_choice= %g\nLL_choice_low= %g\nLL_choice_high= %g\nLL_conf= %g\nLL_conf_corr= %g\nLL_conf_err= %g\nLL_multinom= %g\nLL_RT= %g\n', ...
%     LL_choice,     LL_choice_low,     LL_choice_high,     LL_conf,     LL_conf_corr,     LL_conf_err,     LL_multinom,     LL_RT   );
%************

end

%% some helper functions

% retrieve the full parameter set, given the adjustable and fixed parameters
function param2 = getParam ( param1 , guess , fixed )
    param2(fixed==0) = param1(fixed==0);  %get adjustable parameters from param1
    param2(fixed==1) = guess(fixed==1);   %get fixed parameters from guess
end


% go from choice-wager intersections to conditionals, and marginal on wager
function [pB_A, pA_B, pB, pAB] = int_to_margconds(pAB,pA)

% pAB:  NxM matrix with intersections for all possibilities of A and B
%       e.g. [pRightHigh, pRightLow; 
%               pLeftHigh, pLeftLow]
%       rows of pAB should contain all unique A, columns all unique B
% pA :  Nx1 column vector with probabilities of A (e.g. pRight, pLeft)

% Returns -

% pB_A : NxM matrix with probabilities of B given A i.e. P(B|A)
% pA_B : NxM matrix with probabilities of A given B i.e. P(A|B)
% pAB  : normalized pAB from input (after division by Ptot)
% pB   : calculated marginal probability of B

% force A and B to be column and row vectors respectively
% this is critical
if ~iscolumn(pA), pA = pA'; end

Ptot = sum(pAB(:));
pAB  = pAB ./ Ptot;

% by definition of conditional probability: P(B|A) = P(A n B) / P(A)
pB_A = pAB ./ pA;

% by law of total probability

% pB(1) = pB_A(1,1)*pA(1) + pB_A(2,1)*pA(2) etc.
pB = pA' * pB_A; % do it with matrix multiplication

% by Bayes' rule: A GIVEN B
pA_B = pB_A .* pA ./ pB;

end

% inverse of int_to_margconds
% go from conditional and marginals to intersection
function [pB_A, pAB] = margconds_to_int(pA_B,pA,pB)

% pA_B:  NxM matrix with probabilities of A given B i.e. P(A|B)
%       e.g. [pRight|High, pRight|Low; 
%               pLeft|High, pLeft|Low]
%       rows of pA_B should contain all unique A, columns all unique B
% pB :  1xM row vector with probability of B (e.g. pHigh, pLow)


pAB  = pA_B .* pB;
pB_A = pAB  ./ pA;

% assign Ptot because we maybe want a warning if Ptot~=1
Ptot = sum(pAB(:));
pAB  = pAB ./ Ptot;

end

