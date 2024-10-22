% Miguel Vivar-Lazo
% 10/08/2024

% Wrapper for fitting choice, RT, and wager data to a "serial" model, i.e.
% two sequential 1D drift-diffusion processes, the first for choice and the
% second for binary confidence or 'wager'.
%
% Uses FP4.c by Roozbeh Kiani to compute a numerical solution to the 1D
% Fokker-Planck equation, for calculating RT distributions and log odds
% correct (confidence) in a 1D bounded accumulator (aka DDM)
% 
% The algorithm, based on the method of Chang & Cooper, 1970, was used in:
% Kiani & Shadlen (2009). Representation of Confidence Associated with a
% Decision by Neurons in the Parietal Cortex, Science 324, 759
% and
% Fetsch, Kiani, Newsome, & Shadlen (2014). Effects of Cortical
% Microstimulation on Confidence in a Perceptual Decision.
% Neuron 83, 797–804 (2014).
  

% Code provided for peer review of:
% M. Vivar-Lazo & C.R. Fetsch, Neural basis of concurrent deliberation
% toward a choice and degree of confidence (2024)

clear all
plotflag = 0; % plot the data before fiting?

%% load example dataset (subset of trials from one monkey)
cd ../data
load exampleData_Genji
cd ../code
%*** struct 'data' with fields: ***
% direction (of motion on that trial: 0=right, 180=left)
% coherence (as a probability, unsigned)
% choice (0=left, 1=right)
% RT (in ms)
% correct (1=correct)
% wager (0=low, 1=high)
% cohs (signed coherence, negative=leftward)


%% Plot the data
if plotflag
    conftask = 2;
    RTtask = 1;
    RTCorrOnly = 0;
    parsedData = parseData(data,conftask,RTtask,RTCorrOnly);
    simplePlots(parsedData,unique(data.cohs),conftask,RTtask,0);
    drawnow;
    pause;
end

%% rename for backward compatibility
monkeyData = data;
clear data;


%% Flags for which algorithms to look at 
% Hyperparameters
cohs = unique(monkeyData.cohs);
tEnd = 3;

%% Fit Model! Fun

useUrgency = 1; % For now only for confidence
asymmetrical = 0; % Only for monkey H
wagerOffset = 0; % Only for monkey H

numParams = 5 + useUrgency + asymmetrical + wagerOffset;
extraParams = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set TolX and TolFun to large values (~0.1) and/or MaxIter and MaxFunEvals
% to small values (30-100) to get a quick-and-dirty fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% original settings:
% options = optimset('Display','iter','TolX', 1e-6, 'TolFun', 1e-6, 'MaxIter', 10000, 'MaxFunEvals', 10000);
% 'quick and dirty':
options = optimset('Display','iter','TolX', 1e-1, 'TolFun', 1e-1, 'MaxIter', 100, 'MaxFunEvals', 100);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Making Grid Search
iterationsSearch = 30;
gridDesign = lhsdesign(iterationsSearch, numParams);
gridDesign(:,1) = gridDesign(:,1) .* 30; % Drift Rate 
gridDesign(:,2) = gridDesign(:,2) .* 2; % Symmetrical Bounds for Choice
gridDesign(:,3) = gridDesign(:,3) .* 2; % Asymmetrical Bounds for High Confidence
gridDesign(:,4) = gridDesign(:,4) .* -2; % Asymmetrical Bounds for Low Confidence
gridDesign(:,5) = gridDesign(:,5) .* .400; % mean non-DT time

% Fminbound boundaries
LB = [1  .01 .01 -2   eps];
UB = [40 2   2   -.01 .500];

% Do we need extra params?
extraParamsName = [];
if useUrgency
    extraParams = extraParams+1;
    gridDesign(:,5+extraParams) = gridDesign(:,4+extraParams) .* 2; % linear urgency signal
    LB = [LB 0];
    UB = [UB 2]; 
    extraParamsName{extraParams} = 'urgencySignal';
end
if asymmetrical % For rt distribution
    extraParams = extraParams+1;
    gridDesign(:,5+extraParams) = gridDesign(:,3) .* .400; % mean non-DT time
    LB = [LB 0];
    UB = [UB .450]; 
    extraParamsName{extraParams} = 'asymmetrical';
end
if wagerOffset
    extraParams = extraParams+1;
    gridDesign(:,5+extraParams) = gridDesign(:,4+extraParams) .* .3; % wager offset
    LB = [LB 0];
    UB = [UB .3]; 
    extraParamsName{extraParams} = 'wagerOffset';
end

% Variables for memory storage 
allErrorLL = nan(1, iterationsSearch);
outcomeParamsMat = nan(iterationsSearch, numParams); % the 8 can change

for i = 1:iterationsSearch
    i
    howToCalculateError = 2;
    tic

    params0 = gridDesign(i,:);
    
    % Map with as high as Urgency signal of 5 still works well 
    [tempOutcomeParams, errorLL] = fminsearchbnd(@(params0) Serial_Fitting(monkeyData, cohs, tEnd, howToCalculateError, params0, extraParamsName), params0, LB, UB, options);
    
    % Save some parameters (Error and outcome parameters for that error)
    allErrorLL(i) = errorLL; 
    outcomeParamsMat(i,:) = tempOutcomeParams;
    cd ../results
    save(['genji_Serial_' char(datetime("today")) '.mat'], 'outcomeParamsMat', 'allErrorLL', 'gridDesign');
    cd ../code
    
    % Display time information and iteration value
    display(['Iteration Number =' num2str(i)])
    datetime 
end

[outcomeError_wTheta, b] = min(allErrorLL);
outcomeParams = outcomeParamsMat(b, :);

%% After finding the best outcome generate the Choice, Reaction Time, and Confidence model responses 

data = 0;
modelCoherence = unique(monkeyData.cohs); % -.5:.01:.5; % unique(monkeyData.cohs);
howToCalculateError = -1;
params0 = outcomeParams;
[~, avgChoice, avgRT, avgConf, avgRightHigh, avgRightLow, avgRTHigh, avgRTLow, avgHigh_Corr, avgHigh_Inco] ...
    = Serial_Fitting(data, modelCoherence, tEnd, howToCalculateError, params0, extraParamsName);


%% Plot figures

[avgChoiceData, ~, stdChoiceData]   = behavioralAverages(monkeyData.choice, monkeyData.cohs); % choice
[avgHighChoiceData, ~, stdHighChoiceData]   = behavioralAverages(monkeyData.choice(monkeyData.wager==1), monkeyData.cohs(monkeyData.wager==1)); % choice
[avgLowChoiceData, ~, stdLowChoiceData]   = behavioralAverages(monkeyData.choice(monkeyData.wager==0), monkeyData.cohs(monkeyData.wager==0)); % choice

[avgRTData, ~, stdRTData]           = behavioralAverages(monkeyData.RT, monkeyData.cohs); % RT
[avgRTHighData, ~, stdRTHighData]    = behavioralAverages(monkeyData.RT(monkeyData.wager==1), monkeyData.cohs(monkeyData.wager==1)); % RT High
[avgRTLowData, ~, stdRTLowData]      = behavioralAverages(monkeyData.RT(monkeyData.wager==0), monkeyData.cohs(monkeyData.wager==0)); % RT Low

[avgWagerData, ~, stdWagerData]             = behavioralAverages(monkeyData.wager, monkeyData.cohs); % Wager
[avgWagerCorrData, ~, stdWagerCorrData]     = behavioralAverages(monkeyData.wager(monkeyData.correct==1), monkeyData.cohs(monkeyData.correct==1)); % wager | Correct
[avgWagerIncorrData, ~, stdWagerIncorrData] = behavioralAverages(monkeyData.wager(monkeyData.correct==0), monkeyData.cohs(monkeyData.correct==0)); % wager | incorrect

uniqMonkeyCoh = unique(monkeyData.cohs);

f = figure;
f.Position = [0 500 1400 800];

subplot(2,3,1);
errorbar(uniqMonkeyCoh, avgChoiceData, stdChoiceData, 'k.', 'Markersize', 25); hold on;
plot(modelCoherence, avgChoice, 'k-', 'Color', [.5 .5 .5], 'Linewidth', 2.5); hold off; %Simulation calculations
xlim([-.512 .512]);
xticks([-.512 -.256 0 .256 .512]); 
xticklabels({'', '', '', '', ''})
% ylabel('Proportion rightward choices')

% prettifyGraphs(20, 'normal', 1)

subplot(2,3,2);
errorbar(uniqMonkeyCoh, avgRTData./1000, stdRTData./1000, 'k.', 'Markersize', 25); hold on; %Simulation calculations
% plot(modelCoherence(1:51), avgRT(1:51), 'k-', 'Color', [.5 .5 .5], 'Linewidth', 2.5);
% plot(modelCoherence(52:end), avgRT(52:end), 'k-', 'Color', [.5 .5 .5], 'Linewidth', 2.5); hold off; %Simulation calculations
plot(modelCoherence, avgRT, 'k-', 'Color', [.5 .5 .5], 'Linewidth', 2.5);
xlim([-.512 .512]);
% xticks([uniqueSignCoh([1,2, 6, 10, 11])]); 
xticks([-.512 -.256 0 .256 .512]); 
xticklabels({'', '', '', '', ''})
% ylabel('Reaction time (s)')

% prettifyGraphs(20, 'normal', 1)


subplot(2,3,3);
errorbar(uniqMonkeyCoh, avgWagerData, stdWagerData, 'k.', 'Markersize', 25); hold on; %Simulation calculations
plot(modelCoherence, avgConf, 'k-', 'Color', [.5 .5 .5], 'Linewidth', 2.5); hold off; %Simulation calculations
xlim([-.512 .512]);
% ylabel('Proportion high bets')
ylim([.5 1])
xticks([-.512 -.256 0 .256 .512]); 
xticklabels({'', '', '', '', ''})
% xlabel('Coherence (%)'); 

% prettifyGraphs(20, 'normal', 1)

%% Conditional Plots

% f = figure;
% f.Position = [0 500 500 800];

subplot(2,3,4)
errorbar(uniqMonkeyCoh, avgHighChoiceData, stdHighChoiceData, 'r.', 'Markersize', 25); hold on;
errorbar(uniqMonkeyCoh, avgLowChoiceData, stdLowChoiceData, 'b.', 'Markersize', 25);
plot(modelCoherence, avgRightHigh, '-', 'Color', [255 127 127]./255, 'Linewidth', 2.5);
plot(modelCoherence, avgRightLow, '-', 'Color', [127 127 255]./255, 'Linewidth', 2.5); hold off
ylabel('Proportion rightward choices'); 
legend('high bet', 'low bet');
legend boxoff
xticks([-.512 -.256 0 .256 .512]); 
xticklabels({'-51.2', '-25.6', '0', '25.6', '51.2'})
xlim([-.512 .512])

% prettifyGraphs(20, 'normal', 1)


subplot(2,3,5)
errorbar(unique(monkeyData.cohs), avgRTHighData./1000, stdRTHighData./1000, 'r.', 'Markersize', 25); hold on;
errorbar(unique(monkeyData.cohs), avgRTLowData./1000, stdRTLowData./1000, 'b.', 'Markersize', 25); 
% plot(modelCoherence(1:51), avgRTHigh(1:51), '-', 'Color', [255 127 127]./255, 'Linewidth', 2.5); 
% plot(modelCoherence(52:end), avgRTHigh(52:end), '-', 'Color', [255 127 127]./255, 'Linewidth', 2.5); 
plot(modelCoherence, avgRTHigh, '-', 'Color', [255 127 127]./255, 'Linewidth', 2.5); 
% plot(modelCoherence(1:51), avgRTLow(1:51), '-', 'Color', [127 127 255]./255, 'Linewidth', 2.5); 
% plot(modelCoherence(52:end), avgRTLow(52:end), '-', 'Color', [127 127 255]./255, 'Linewidth', 2.5); hold off
plot(modelCoherence, avgRTLow, '-', 'Color', [127 127 255]./255, 'Linewidth', 2.5); hold off
xlabel('Motion strength (% coh)'); 
ylabel('Reaction time (s)'); 
% legend('High Bet', 'Low Bet');
xticks([-.512 -.256 0 .256 .512]); 
xticklabels({'-51.2', '-25.6', '0', '25.6', '51.2'})
xlim([-.512 .512])

% prettifyGraphs(20, 'normal', 1)


subplot(2,3,6)
errorbar(uniqMonkeyCoh, avgWagerCorrData, stdWagerCorrData, '.', 'Color', [0.4940 0.1840 0.5560], 'Markersize', 25); hold on;
errorbar(uniqMonkeyCoh(2:end-1), avgWagerIncorrData(2:end-1), stdWagerIncorrData(2:end-1), '.', 'Color', [0.4660 0.6740 0.1880], 'Markersize', 25);
plot(modelCoherence, avgHigh_Corr, '-', 'Color', [0.4940 0.1840 0.5560], 'Linewidth', 2.5); hold on;
plot(modelCoherence, avgHigh_Inco, '-', 'Color', [0.4660 0.6740 0.1880], 'Linewidth', 2.5); 
ylabel('Proportion high bets'); 
legend('correct', 'incorrect');
legend boxoff
xticks([-.512 -.256 0 .256 .512]); 
xticklabels({'-51.2', '-25.6', '0', '25.6', '51.2'})
xlim([-.512 .512])
ylim([0 1])

% prettifyGraphs(20, 'normal', 1)


%% Run a BIC and AIC test
if 1
    indexModel = 0;
    aicVector = nan(1,10);
    bicVector = nan(1,10);
    nameOfModel = cell(1,10);
end

% AIC & BIC
[aic, bic] = aicbic(-outcomeError_wTheta , length(outcomeParams), length(monkeyData.choice)); 

indexModel = indexModel+1;
aicVector(indexModel) = aic;
bicVector(indexModel) = bic;
nameOfModel{indexModel} = 'Genji Serial';

%% Calculate LogLL for Wager|Accuracy
% avgHigh_Corr

errorWagerHigh = nan(1,length(modelCoherence));
errorWagerLow = nan(1,length(modelCoherence));

for i = 1:length(modelCoherence)
    % Find the empirical data
    empData_High  = sum(monkeyData.wager == 1 & monkeyData.correct == 1 & monkeyData.cohs == modelCoherence(i));
    empData_Low  = sum(monkeyData.wager == 0 & monkeyData.correct == 1 & monkeyData.cohs == modelCoherence(i));
    errorWagerHigh(i) = empData_High*log(avgHigh_Corr(i)) + empData_Low * log(1-avgHigh_Corr(i));
    % Low
    empData_High  = sum(monkeyData.wager == 1 & monkeyData.correct == 0 & monkeyData.cohs == modelCoherence(i));
    empData_Low  = sum(monkeyData.wager == 0 & monkeyData.correct == 0 & monkeyData.cohs == modelCoherence(i));
    errorWagerLow(i) = empData_High*log(avgHigh_Inco(i)) + empData_Low * log(1-avgHigh_Inco(i));

end

% Final values
errorWagerHigh = -sum(errorWagerHigh);
errorWagerLow = -sum(errorWagerLow);


%% Save the figure and workspace just in case
cd ../results
saveas(gcf,'data_wFit_Serial','fig')
clear f
save workspace_Serial.mat
cd ../code

