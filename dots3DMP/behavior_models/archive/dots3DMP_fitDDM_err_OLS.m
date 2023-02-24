function [err,fit] = dots3DMP_fitDDM_err_OLS(param, data, options)

global call_num

ntrials = length(data.heading);

cohs = unique(data.coherence); % visual coherence levels
hdgs = [-10 -5 -2.5 -1.25 -eps eps 1.25 2.5 5 10]; % heading angles
                            % (map fn seems to require even number of diff 
                            % levels (or just no zero), so we use +/- eps)
duration = 2000; % stimulus duration (ms)
dur = ones(ntrials,1) * duration;

kves = param(1);
kvis = param(2)*cohs;
B = abs(param(3)); % don't accept negative bound heights

% what about sigma? assume 1 for now (later may need not just separate 
% sigmas but full Drugowitsch paramaterization):
sigmaVes = 1; % std of momentary evidence
sigmaVis = [1 1]; % allow for separate sigmas for condition, coherence

Tnd = 300; % non-decision time (or make this the mean of a dist?)

maxdur = duration;
% assume the mapping is based on an equal amount of experience with the 
% *three* levels of reliability (ves, vis-low, vis-high) hence k and sigma
% are their averages
k = mean([kves kvis']);
sigma = mean([sigmaVes sigmaVis]);
[~, ~, logOddsCorrMap, tAxis, vAxis] = makeLogOddsCorrMap_3DMP(hdgs,k,B,sigma,maxdur,0);
% uses Fokker-Planck equation to propagate the probability density of the DV,
% as in Kiani & Shadlen 2009. Required for readout of confidence, although
% a simpler heuristic could be used (conf proportional to accum evidence)

% create acceleration and velocity profiles (arbitrary for now)
% SJ 04/2020
% Hou et al. 2019, peak vel = 0.37m/s, SD = 210ms
vel = normpdf(1:duration,duration/2,210);
vel = 0.37*vel./max(vel);
acc = gradient(vel)*1000; % multiply by 1000 to get from m/s/ms to m/s/s

% normalize
vel = vel./max(vel);
acc = abs(acc./max(acc)); % (and abs)


%% bounded evidence accumulation

choice = nan(ntrials,1);
RT = nan(ntrials,1);
finalV = nan(ntrials,1);
hitBound = zeros(1,ntrials);
logOddsCorr = nan(ntrials,1);
expectedPctCorr = nan(ntrials,1);
conf = nan(ntrials,1);

hdg = data.heading;
coh = data.coherence;
delta = data.delta;
modality = data.modality;


tic
for n = 1:ntrials

    switch modality(n)
        case 1
            mu = acc .* kves * sind(hdg(n)); % mean of momentary evidence
%             dv = [0, cumsum(normrnd(mu,sigmaVes,1,dur(n)-1))];
            dv = [0, cumsum(normrnd(mu,sigmaVes))];
        case 2
            mu = vel .* kvis(cohs==coh(n)) * sind(hdg(n));
%             dv = [0, cumsum(normrnd(mu,sigmaVis(cohs==coh(n)),1,dur(n)-1))];
            dv = [0, cumsum(normrnd(mu,sigmaVis(cohs==coh(n))))];
        case 3
            % positive delta defined as ves to the left, vis to the right
%             muVes = kves                      * sind(hdg(n) - delta(n)/2);
            muVes = acc .* kves               * sind(hdg(n) - delta(n)/2);
            muVis = vel .* kvis(cohs==coh(n)) * sind(hdg(n) + delta(n)/2);
            
            wVes = sqrt( kves^2 / (kvis(cohs==coh(n))^2 + kves^2) );
            wVis = sqrt( kvis(cohs==coh(n))^2 / (kvis(cohs==coh(n))^2 + kves^2) );

            mu = wVes.*muVes + wVis.*muVis;
            % here the DV is a sample from a dist with mean = weighted sum
            % of means. thus the variance is the weighted sum of variances
            % (error propagation formula):
            sigmaComb = sqrt(wVes.^2 .* sigmaVes^2 + wVis.^2 .* sigmaVis(cohs==coh(n))^2); % assume zero covariance
%             dv = [0, cumsum(normrnd(mu,sigmaComb,1,dur(n)-1))];
            dv = [0, cumsum(normrnd(mu,sigmaComb))]; % if mu is a vector
    end

    cRT = find(abs(dv)>=B, 1);
    if isempty(cRT) % did not hit bound
        RT(n) = dur(n) + Tnd;
        finalV(n) = dv(dur(n));
        hitBound(n) = 0;
    else % hit bound
        RT(n) = cRT + Tnd;
        finalV(n) = B*sign(dv(cRT));
        hitBound(n) = 1;
    end    
    choice(n) = sign(finalV(n));
    
    % use map to look up log-odds that the motion is rightward
    diffV = abs(vAxis-finalV(n));
    diffT = abs(tAxis-RT(n));
        
    thisV = find(diffV==min(diffV));
    thisT = find(diffT==min(diffT));
    logOddsCorr(n) = logOddsCorrMap(thisV(1), thisT(1));
    expectedPctCorr(n) = logistic(logOddsCorr(n)); % convert to pct corr
    conf(n) = 2*expectedPctCorr(n) - 1; % convert to 0..1

end
toc

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


%% calculate error

% What we want is to use maximum-likelihood estimation (MLE), so we need a
% way to calculate the likelihood.
% Here's a refresher (see also likelihood_tutorial.m):

% The likelihood is often written P(data|params), which looks like a 
% conditional density where data varies and params is fixed, but it's
% really kind of the opposite: it's a function of the params for fixed data.

% Confusingly, we don't actually need an equation for the full likelihood
% *function* -- that's the Nparams-dimensional surface which is mapped out
% over many iterations of our fitting routine. The optimization just needs
% to sample the surface sufficiently well to find the maximum.

% What each iteration needs is a functional form for the conditional
% density p(data|params), where data is a set of N trials x 3 observations,
% choice, RT, and confidence. In other words we need a concise way to
% describe the stochastic process that gives rise to the observations,
% given a set of params.

% In likelihood_tutorial.m, the stochastic process is Poisson, the
% observation is a spike count, and params is one-dimensional: some
% stimulus value theta. The probability p(obs|theta) for a given theta
% is either calculated from the equation for the Poisson probability
% mass function, or estimated using a 'histogram method' by generating many
% samples from this distribution.

% This approach seems like it might be amenable to what we're doing here,
% since we're generating the model fits/predictions through monte carlo
% simulation anyway. But there's a key difference: In likelihood_tutorial,
% there was only a scalar observation (the spike count) and an IID set of
% trials (poissrnd(theta)). Here, we have a vector of observations (choice,
% RT, conf) and a non-iid set of trials, each defined by the 4-vector:
% [heading, coh, delta, modality].

% In the Kiani DDM fitting code, it was possible to accumulate (log)
% likelihood trial by trial, because for any combination of coh and
% duration, the model specified a probability of rightward, leftward, and
% sure bet *for that trial*. To replicate this in our monte carlo sim-style
% fitting, we'd have to simulate every combination of hdg, mod, coh, delta
% many times to estimate these probabilities. This is probably too slow to
% be viable, but I may try it anyway.



% For now we'll use least squares just to get it working at all, but need a
% better long term fix... (least squares is equivalent to maximum
% likelihood only in the case of Gaussian error dist, which maybe is close
% enough?)


% normalize before calculating error:
% shift choice back to 0..1, normalize RT, and conf already is 0..1
choiceM = choice - 1;
choiceD = data.choice - 1;

RTm = RT/max(RT);
RTd = data.RT/max(data.RT);

err_choice = sum((choiceM-choiceD).^2);
err_RT = sum((RTm-RTd).^2);
err_conf = sum((conf-data.conf).^2);

err = err_choice + err_RT + err_conf;





%% print progress report!
fprintf('\n\n\n****************************************\n');
fprintf('run %d\n', call_num);
fprintf('\tkves= %g\n\tkvisMult= %g\n\tB= %g\n', param(1), param(2), param(3));
fprintf('err: %f\n', err);
if options.ploterr && strcmp(options.fitMethod,'fms')
    figure(options.fh);
    plot(call_num, err, '.','MarkerSize',12);
    drawnow;
end
call_num = call_num + 1;



end


