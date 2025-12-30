                                  
% Miguel's (Hanzo) data

% KEEPING THIS ALIVE TO DIAGNOSE FIT ISSUE: why by eye a set of params
% looks better but gives worse likelihoods

clear all; close all
load Hanzo_Confidence

validTrials = ~isnan(signCoherence) & ~isnan(RT) & ~isnan(choice) & ...
              ~isnan(correct) & ~isnan(confidenceChoice);
          
% rename some vars:
D.strength = signCoherence(validTrials);
D.dur = round(RT(validTrials)*1000);
D.choice = choice(validTrials); 
D.correct = correct(validTrials);
D.conf = confidenceChoice(validTrials);


options.fitMethod = 'fms';
% options.fitMethod = 'global';
% options.fitMethod = 'multi';
% options.fitMethod = 'pattern';


% params: 
%        1 2 3
%       [k B theta alpha]
fixed = [0 0 0 0];


% % initial guess (or hand-tuned params)
% k = 0.5; % sensitivity parameter
% B = 30; % bound height
% theta = 1.6; % criterion (in log odds correct) for betting high
% alpha = 0.1; % base rate of low-bet choices


% best-fit values, but visually worse! get to the bottom of this!
k = 0.62; % sensitivity parameter
B = 14.65; % bound height
theta = 1.21; % criterion (in log odds correct) for betting high
alpha = 0.12; % base rate of low-bet choices


guess = [k B theta alpha];

% ************************************
% set all fixed to 1 for hand-tuning:
fixed(:)=1;
% ************************************


options.feedback = 1; % 1 = text output to cmd window, 2 = that and plot LL across runs
options.plot = 0; % plot the marginal PDFs, logOddsCorr map, and high/low bet regions (only makes sense for fixed(:)=1)


[X, LL_final, data, fit] = fitDDM_wConfidence_simple(D,options,guess,fixed);


%% plot data and compare to fits

ucoh = unique(data.strength);
pRight = nan(length(ucoh),1);
pRight_se = nan(length(ucoh),1);
pHigh = nan(length(ucoh),1); % high bet (conf)
pHigh_se = nan(length(ucoh),1);
N = nan(length(ucoh),1);
    % SEs of the expected proportions from the model are probably not
    % meaningful. What matters are the SEs of the parameters themselves
for c = 1:length(ucoh)
    I = data.strength==ucoh(c);
    pRight(c) = sum(data.choice(I)==1) / sum(I); % 1 is rightward
    pRight_se(c) = sqrt(pRight(c)*(1-pRight(c)) / sum(I)); % formula for standard error of a proportion
    pHigh(c) = sum(data.conf(I)==1) / sum(I); % 1 is high bet
    pHigh_se(c) = sqrt(pHigh(c)*(1-pHigh(c)) / sum(I)); 
    N(c) = sum(I);
end

% now the model
ucoh_fit = unique(fit.strength);
pRight_model = nan(length(ucoh_fit),1);
pHigh_model = nan(length(ucoh_fit),1); % high bet (conf)
for c = 1:length(ucoh_fit)
    I = fit.strength==ucoh_fit(c);
    pRight_model(c) = mean(fit.expectedPright(I));
    pHigh_model(c) = mean(fit.expectedPhigh(I));
        % should separate this out by duration for a more fine-grained look.
        % Then there's also confidence conditioned on correct/error, and
        % accuracy conditioned on high/low bet...
end

figure; set(gcf,'Position',[86 925 1070 420]);
subplot(1,2,1); errorbar(ucoh,pRight,pRight_se,'bo-');
hold on; plot(ucoh_fit,pRight_model,'c--');
xlabel('motion strength (%coh)');
ylabel('proportion rightward choices');

subplot(1,2,2); errorbar(ucoh,pHigh,pHigh_se,'ro-');
hold on; plot(ucoh_fit,pHigh_model,'m--'); ylim([0 1]);
xlabel('motion strength (%coh)');
ylabel('proportion high bets');

