function [err,fit] = dots3DMP_fit_2Dacc_err(param, guess, fixed, data, options)

global call_num

% set parameters for this run based on guess and fixed params flag
param = getParam(param, guess, fixed);

if isfield(data,'PDW')
    data.conf = data.PDW;
end

ntrials = length(data.heading);
            
mods = unique(data.modality);
cohs = unique(data.coherence); % visual coherence levels
hdgs = unique(data.heading);
deltas = unique(data.delta);

modality = data.modality;
coh      = data.coherence;
hdg      = data.heading;
delta    = data.delta;

duration = 2000; % stimulus duration (ms)
dur = ones(ntrials,1) * duration;

ks    = param(1);
% sigma = param(2);
% B     = abs(param(3)); % don't accept negative bound heights
% muTnd = param(4);
% if numel(param)==5, theta = param(5); end % only relevant for PDW

% sigmaVes = sigma; % std of momentary evidence
% sigmaVis = [sigma sigma]; % allow for separate sigmas for condition, coherence

sigmaVes = param(2);
sigmaVis = param([3 4]);
B     = abs(param(5)); % don't accept negative bound heights
muTnd = param(6);

if numel(param)==7, theta = param(7); end % only relevant for PDW

kves = ks;
kvis = ks * [2 4.5]/3; % should 2/4.5 be another set of parameters? - relationship between ves and vis is fixed otherwise

sdTnd = 60; % fixed SD
% sdTnd = 0;
Tnds = muTnd + randn(ntrials,1).*sdTnd;


% % assume the mapping is based on an equal amount of experience with the 
% % *three* levels of reliability (ves, vis-low, vis-high) hence k and sigma
% % are their averages
% k = mean([kves kvis]);

% assume the mapping is based on an equal amount of experience with the 
% *three* levels of reliability (ves, vis-low, vis-high) hence k and sigma
% are their averages
k = mean([kves kvis]);

            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            
% TEMP: simulate a few things
kves = 0.33*ks; % vestibular deficit [try this before vs after mapping!]
% B = B*1.5; % bound height change

            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP



R.t = 0.001:0.001:duration/1000;
R.Bup = B;
R.drift = k * sind(hdgs(hdgs>=0)); % takes only unsigned drift rates
R.lose_flag = 1; % need pdf of the losing race(s)
R.plotflag = 0; % 1 = plot, 2 = plot and export_fig

P = images_dtb_2d(R);

% create acceleration and velocity profiles (arbitrary for now)
% SJ 04/2020
% Hou et al. 2019, peak vel = 0.37m/s, SD = 210ms

ampl = 0.16; % movement in metres
pos = normcdf(1:duration,duration/2,0.14*duration)*ampl;
vel = gradient(pos)*1000; % metres/s
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
hitBound        = zeros(1,ntrials);
logOddsCorr     = nan(ntrials,1);
expectedPctCorr = nan(ntrials,1);
conf            = nan(ntrials,1);


% ME is now a draw from bivariate normal with mean vector Mu and covariance matrix V
% start with correlation matrix:
S = [1 -1/sqrt(2) ; -1/sqrt(2) 1];
% -1/sqrt(2) is the correlation for our version of images_dtb;
% can only be changed with an update to the flux file

for n = 1:ntrials
    Tnd = Tnds(n) / 1000; % non-decision time for nth trial in seconds

    switch modality(n)
        case 1
            mu = acc .* kves * sind(hdg(n)) / 1000; % mean of momentary evidence
                % (I'm guessing drift rate in images_dtb is per second, hence div by 1000)
            s = [sigmaVes sigmaVes]; % standard deviaton vector (see below)
        case 2
            mu = vel .* kvis(cohs==coh(n)) * sind(hdg(n)) / 1000;
            s = [sigmaVis(cohs==coh(n)) sigmaVis(cohs==coh(n))];
        case 3
            % positive delta defined as ves to the left, vis to the right
            muVes = acc .* kves               * sind(hdg(n)-delta(n)/2) / 1000;
            muVis = vel .* kvis(cohs==coh(n)) * sind(hdg(n)+delta(n)/2) / 1000;

            % optimal weights (Drugo et al.) 
            wVes = sqrt( kves^2 / (kvis(cohs==coh(n))^2 + kves^2) );
            wVis = sqrt( kvis(cohs==coh(n))^2 / (kvis(cohs==coh(n))^2 + kves^2) );
            
            
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            
%            wVes = 0.1; wVis = 0.9; % vis overweighting, no coh shift
            % OR
            if coh(n)==cohs(end) % use high coh as placeholder for TMS/Pt

            %    wVes = 0.5; wVis = 0.5;
                wVes = rand; wVis = 1-wVes; % TMS

            else
                wVes = 0.1; wVis = 0.9;
                

                % OR
% %                 kvesPt = 0.1*ks; % vestibular deficit
% %                 wVes = sqrt( kvesPt^2 / (kvis(cohs==coh(n))^2 + kvesPt^2) );
% %                 wVis = sqrt( kvis(cohs==coh(n))^2 / (kvis(cohs==coh(n))^2 + kvesPt^2) );
                % nope, this doesn't capture it because sensitivity isn't right
            end
            
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            
            mu = wVes.*muVes + wVis.*muVis;

            % the DV is a sample from a dist with mean = weighted sum of
            % means. thus the variance is the weighted sum of variances
            % (error propagation formula):
            sigmaComb = sqrt(wVes.^2 .* sigmaVes^2 + wVis.^2 .* sigmaVis(cohs==coh(n))^2); % assume zero covariance
            s = [sigmaComb sigmaComb];
            
            
            
%             % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
%             % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
%             % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
%             % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
%             % okay can do this here, just repeat above
%             if coh(n)==cohs(end) % use high coh as placeholder for TMS/Pt
%                 kves = 0.1*ks; % vestibular deficit
%                 muVes = acc .* kves               * sind(hdg(n)-delta(n)/2) / 1000;
%                 muVis = vel .* kvis(cohs==coh(n)) * sind(hdg(n)+delta(n)/2) / 1000;
%                 wVes = sqrt( kves^2 / (kvis(cohs==coh(n))^2 + kves^2) );
%                 wVis = sqrt( kvis(cohs==coh(n))^2 / (kvis(cohs==coh(n))^2 + kves^2) );
%                 mu = wVes.*muVes + wVis.*muVis;
%                 sigmaComb = sqrt(wVes.^2 .* sigmaVes^2 + wVis.^2 .* sigmaVis(cohs==coh(n))^2); % assume zero covariance
%                 s = [sigmaComb sigmaComb];
%                 
%                 % restore kves for next trial
%                 kves = ks;
%             end
%             % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
%             % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
%             % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
%             % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            
    end

%     Mu = [mu,-mu]; % mean vector for 2D DV
    Mu = [mu; -mu]';
    
    % convert correlation to covariance matrix
    V = diag(s)*S*diag(s);
%     dv = [0 0; cumsum(mvnrnd(Mu,V,dur(n)-1))]; % bivariate normrnd
    dv = [0 0; cumsum(mvnrnd(Mu,V))];
    % because Mu is signed according to heading (positive=right),
    % dv(:,1) corresponds to evidence favoring rightward, not evidence
    % favoring the correct decision (as in Kiani eqn. 3 and images_dtb)
    
    % decision outcome
    cRT1 = find(dv(:,1)>=B, 1);
    cRT2 = find(dv(:,2)>=B, 1);
    % the options are:
    % (1) only right accumulator hits bound,
    if ~isempty(cRT1) && isempty(cRT2)
        RT(n) = cRT1/1000 + Tnd;
        finalV(n) = dv(cRT1,2); % only 1 hit, so 2 is the loser
        hitBound(n) = 1;
        choice(n) = 1;
    % (2) only left accumulator hits bound, 
    elseif isempty(cRT1) && ~isempty(cRT2)
        RT(n) = cRT2/1000 + Tnd;
        finalV(n) = dv(cRT2,1); % only 2 hit, so 1 is the loser
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
        whichWon = dv(end,:)==max(dv(end,:));
        
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
        finalV(n) = dv(min([cRT1 cRT2]),~whichWon); % the not-whichWon is the loser
        hitBound(n) = 1;
        a = [1 -1];
        choice(n) = a(whichWon);
    end
        
    % use map to look up log-odds that the motion is rightward
    diffV = abs((P.y+B)-finalV(n));
    diffT = abs(R.t-RT(n));
        
    thisV = find(diffV==min(diffV));
    thisT = find(diffT==min(diffT));
    logOddsCorr(n) = P.logOddsCorrMap(thisV(1), thisT(1));
    
    if ~exist('theta','var')
        expectedPctCorr(n) = logistic(logOddsCorr(n)); % convert to pct corr
        conf(n) = 2*expectedPctCorr(n) - 1; % convert to 0..1
    else
        conf(n) = logOddsCorr(n) > theta;
    end
end


            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            
% TEMP SUPER KLUGE (need to add alpha param)
conf = conf-0.08;

            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
            % TEMP TEMP TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP

choice(choice==0) = sign(randn); % not needed under usual circumstances
choice(choice==1) = 2; choice(choice==-1) = 1; % 1=left, 2=right

% output var 'fit' gets the same values as data for the conditions, but the
% simulated trial outcomes for the observables!
fit.heading = data.heading;
fit.coherence = data.coherence;
fit.modality = data.modality;
fit.delta = data.delta;
fit.choice = choice;
fit.RT = RT;
fit.conf = conf;

fit.correct = (choice==2 & data.heading>0) | (choice==1 & data.heading<0);
data.correct = (data.choice==2 & data.heading>0) | (data.choice==1 & data.heading<0);

%% calculate error

% Do ibs_basic?

% % likelihood of rightward choice on each trial, under binomial assumptions

% based on Palmer et al. 2005 coh version, but this ignores mod/coh in 3DMP
% task
% Pr_model = 1 ./ (1 + exp(-2*k*B*sind(hdg)));

% RTdata = (data.RT - min(data.RT)) ./ (max(data.RT)-min(data.RT));
% RTfit = (fit.RT - min(fit.RT)) ./ (max(fit.RT)-min(fit.RT));
% 
% confdata = (data.conf - min(data.conf)) ./ (max(data.conf)-min(data.conf));
% conffit = (fit.conf - min(fit.conf)) ./ (max(fit.conf)-min(fit.conf));

RTdata = data.RT .* 1000;
RTfit  = fit.RT  .* 1000;

confdata = data.conf .* 100;
conffit  = fit.conf .* 100;

choiceD  = logical(data.choice-1);

pRight_trialwise = nan(ntrials,1);
pHighBet_trialwise  = nan(ntrials,1);

n = nan(length(mods),length(cohs),length(deltas),length(hdgs));

pRight_model = n; 
pHighBet_model = n;
RTmean_fit = n;     RTmean_data = n;    sigmaRT = n;
confMean_fit = n;   confMean_data = n;  sigmaConf = n;
nCor = n;

nRight_data   = n;
nHighBet_data = n;

for m = 1:length(mods)
for c = 1:length(cohs)
for d = 1:length(deltas)

    for h = 1:length(hdgs)
        J = data.modality==mods(m) & data.coherence==cohs(c) & data.heading==hdgs(h) & data.delta==deltas(d);
        
        n(m,c,d,h) = sum(J);
        
        if sum(J)>0
            nRight_data(m,c,d,h) = sum(J & choiceD);
        end
        
        usetrs = data.correct | data.heading==0;
        nCor(m,c,d,h) = sum(J & usetrs); % use only correct trials for RT and conf (sacc EP) fits
        
        pRight_model(m,c,d,h) = sum(J & fit.choice==2) / n(m,c,d,h); % 2 is rightward
        pRight_trialwise(J) = pRight_model(m,c,d,h);
        
        if options.RTtask
            RTmean_fit(m,c,d,h) = mean(RTfit(J & usetrs));
            RTmean_data(m,c,d,h) = mean(RTdata(J & usetrs));
            sigmaRT(m,c,d,h) = std(RTdata(J & usetrs))/sqrt(nCor(m,c,d,h));
        end
        
        if options.conftask == 1 % sacc endpoint
            confMean_fit(m,c,d,h) = mean(conffit(J & usetrs));
            confMean_data(m,c,d,h) = mean(confdata(J & usetrs));
            sigmaConf(m,c,d,h) = std(confdata(J & usetrs))/sqrt(nCor(m,c,d,h));
        elseif options.conftask == 2 % PDW
            pHighBet_model(m,c,d,h) = sum(J & fit.conf==1) / n(m,c,d,h); % 1 is high
            pHighBet_trialwise(J) = pHighBet_model(m,c,d,h);
            
            if sum(J)>0
                nHighBet_data(m,c,d,h) = sum(J & data.conf==1);
            end
        end 
    end
    
end
end
end

% keyboard
%%
% kluge to avoid log(0) issues
% pRight_model(pRight_model==0) = min(pRight_model(pRight_model~=0)); 
% pRight_model(pRight_model==1) = max(pRight_model(pRight_model~=1));
pRight_trialwise(pRight_trialwise==0) = eps; 
pRight_trialwise(pRight_trialwise==1) = 1-eps;

% log likelihood of choice is summed log likelihood across all trials (log
% probability of observing choice given model)
LL_choice = (sum(log(pRight_trialwise(choiceD))) + sum(log(1-pRight_trialwise(~choiceD))));

% or use binomial directly
% X = nRight_data(:); N = n(:); P = pRight_model(:);
% P(P==0) = eps; P(P==1) = 1-eps;
% 
% L_choice = binopdf(X(~isnan(X)),N(~isnan(X)),P(~isnan(X)));
% L_choice(L_choice==0) = eps;
% L_choice(L_choice==1) = 1-eps;
% LL_choice = nansum(log(L_choice));

% log likelihood of mean RTs, under Gaussian approximation
LL_RT = 0;
if options.RTtask
    L_RT = 1./(sigmaRT(:).*sqrt(2*pi)) .* exp(-(RTmean_fit(:) - RTmean_data(:)).^2 ./ (2*sigmaRT(:).^2));
%     L_RT(L_RT==0) = min(L_RT(L_RT~=0));
    L_RT(L_RT==0) = eps;
    LL_RT = nansum(log(L_RT(:))); % sum over all conditions
end

% likelihood of mean confidence (for SEP), or log probabilities for PDW
LL_conf = 0;
if options.conftask == 1 % SEP
    L_conf = 1./(sigmaConf*sqrt(2*pi)) .* exp(-(confMean_fit - confMean_data).^2 ./ (2*sigmaConf.^2));
%     L_conf(L_conf==0) = min(L_conf(L_conf~=0));
    L_conf(L_conf==0) = eps;
    LL_conf = nansum(log(L_conf(:))); % sum over all conditions
    
elseif options.conftask ==2 % PDW
%     PHighBet_trialwise(PHighBet_trialwise==0) = min(PHighBet_trialwise(PHighBet_trialwise~=0)); 
%     PHighBet_trialwise(PHighBet_trialwise==1) = max(PHighBet_trialwise(PHighBet_trialwise~=1));
    pHighBet_trialwise(pHighBet_trialwise==0) = eps; 
    pHighBet_trialwise(pHighBet_trialwise==1) = 1-eps;
    PDW = logical(data.conf);
    LL_conf = (sum(log(pHighBet_trialwise(PDW))) + sum(log(1-pHighBet_trialwise(~PDW)))); 
    
%     X = nHighBet_data(:); N = n(:); P = pHighBet_model(:);
%     P(P==0) = eps; P(P==1) = 1-eps;
%     
%     L_conf = binopdf(X(~isnan(X)),N(~isnan(X)),P(~isnan(X)));
%     L_conf(L_conf==0) = eps;
%     L_conf(L_conf==1) = 1-eps;
%     LL_conf = nansum(log(L_conf));
end

% err = -LL_choice;
err = -(LL_choice + LL_conf + LL_RT);
% err = -(LL_choice + LL_conf);

%% print progress report!
fprintf('\n\n\n****************************************\n');
fprintf('run %d\n', call_num);
if exist('theta','var')
%     fprintf('\tks= %g\n\tsigma= %g\n\tB= %g\n\tTnd= %g\n\ttheta= %g\n', param(1), param(2), param(3), param(4), param(5));
    fprintf('\tks= %g\n\tsigmas= %g, %g, %g\n\tB= %g\n\tTnd= %g\n\ttheta= %g\n', param(1), param(2), param(3), param(4), param(5), param(6), param(7));

else
%     fprintf('\tks= %g\n\tsigma= %g\n\tB= %g\n\tTnd= %g\n', param(1), param(2), param(3), param(4));
fprintf('\tks= %g\n\tsigmas= %g, %g, %g\n\tB= %g\n\tTnd= %g\n', param(1), param(2), param(3), param(4), param(5),param(6));

end
fprintf('err: %f\n', err);
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
