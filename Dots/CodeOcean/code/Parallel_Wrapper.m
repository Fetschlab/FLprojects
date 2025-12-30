% Miguel Vivar-Lazo
% 10/08/2024
% 
% Wrapper for fitting choice, RT, and wager data with a "parallel" model
% (2D accumulator, aka anticorrelated race) 

% Uses FP_Drugowitsch.c, which is our implementation of the method in:
% Shan, Moreno-Bote, and Drugowitsch (2019) Family of closed-form solutions
% for two-dimensional correlated diffusion processes. Phys Rev E100, 032132

% Model structure inspired by and very similar to:
% Kiani, Corthell & Shadlen (2014). Choice Certainty Is Informed by Both 
% Evidence and Decision Time. Neuron 84, 1329–1342.
  
% Code provided for peer review of:
% M. Vivar-Lazo & C.R. Fetsch, Neural basis of concurrent deliberation
% toward a choice and degree of confidence (2024)


clear all
plotflag = 0; % plot the data before fiting?
dbstop if error

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
allResponseVariables = 0;
gridSearch = 1;

% Hyperparameters
k = 4; % -cos(pi/k), number of images (grid search: k=[3, 4, 5, 7, 10] ~= rho=[-.5, -.707, -.809, -.901, -.951])
bound = -.025; % Closer to 0. Remember that starting points are near lowest negative numbers
totalTime = 3; %in seconds
deltaT = .01;
ti = 0:.01:totalTime;

%% Try fitting the Choice and RT even Wager (Using Global Search)
% Fit Choice, Reaction Time, and Wager
% This Runs just as fast as CDF so just run this in the the parameters
% placed here

useUrgency = 1; % Both monekys
asymmetrical = 1; % Only monkey H
wagerOffset = 1; % Only monkey H
thetaTimeCutoff = 0; 
useUrgencyHalfMax = 0; % Only monkey G

numParams = 4 + useUrgency + asymmetrical + wagerOffset + thetaTimeCutoff + useUrgencyHalfMax;
extraParams = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set TolX and TolFun to large values (~0.1) and/or MaxIter and MaxFunEvals
% to small values (~10) to get a quick-and-dirty fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% original settings:
% options = optimset('Display','iter','TolX', 1e-6, 'TolFun', 1e-6, 'MaxIter', 10000, 'MaxFunEvals', 10000);
% 'quick and dirty':
options = optimset('Display','iter','TolX', 1e-1, 'TolFun', 1e-1, 'MaxIter', 10, 'MaxFunEvals', 10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if gridSearch
    iterationsSearch = 30;
    gridDesign = lhsdesign(iterationsSearch, numParams);
    gridDesign(:,1) = gridDesign(:,1) .* 20; % Drift Rate 
    gridDesign(:,2) = gridDesign(:,2) .* -1.5 - .025; % start position
    gridDesign(:,3) = gridDesign(:,3) .* .400; % mean non-DT time
    gridDesign(:,4) = gridDesign(:,4) .* 2; % theta
    % Fminbound boundaries
    LB = [0    -1.8  .001 0];
    UB = [20   -.025 .450 10];

    % Do we need extra params?
    extraParamsName = [];
    if useUrgency
        extraParams = extraParams+1;
        gridDesign(:,4+extraParams) = gridDesign(:,4+extraParams) .* 8; % urgency signal (For linear use 6, for Urgency use something much higher ~10)
        LB = [LB 0];
        UB = [UB 20]; 
        extraParamsName{extraParams} = 'urgencySignal';
    end
    if useUrgencyHalfMax
        extraParams = extraParams+1;
        gridDesign(:,4+extraParams) = gridDesign(:,4+extraParams) .* 8; % urgency signal
        LB = [LB eps];
        UB = [UB 20]; 
        extraParamsName{extraParams} = 'useUrgencyHalfMax';
    end
    if asymmetrical % For rt distribution
        extraParams = extraParams+1;
        gridDesign(:,4+extraParams) = gridDesign(:,3); % mean non-DT time
        LB = [LB 0];
        UB = [UB .450]; 
        extraParamsName{extraParams} = 'asymmetrical';
    end
    if wagerOffset
        extraParams = extraParams+1;
        gridDesign(:,4+extraParams) = gridDesign(:,4+extraParams) .* .3; % mean non-DT time
        LB = [LB 0];
        UB = [UB .3]; 
        extraParamsName{extraParams} = 'wagerOffset';
    end
    if thetaTimeCutoff
        extraParams = extraParams+1;
        gridDesign(:,4+extraParams) = gridDesign(:,4+extraParams) .* .01 - .005; % mean non-DT time
        LB = [LB -.01];
        UB = [UB .01]; 
        extraParamsName{extraParams} = 'thetaTimeCutoff';
    end
else 
    % Design search space for only wager theta
    vectorDesign = lhsdesign(iterationsSearch, 1);
    vectorDesign = vectorDesign .* 5;
end

% Variables for memory storage 
allErrorLL_wTheta = nan(1, iterationsSearch);
outcomeParamsMat_wTheta = nan(iterationsSearch, numParams); % the 8 can change
%

whichSolve = 3; % How to solve it?
for i = 1:iterationsSearch
    tic

    if gridSearch
        params0_wTheta = gridDesign(i,:);
    else
        % Add some noise to best parameters from fitting only Choice and RT
        noiseOutcomeParams = normrnd(outcomeParams, sqrt((outcomeParams.*.1).^2));
        params0_wTheta = [noiseOutcomeParams vectorDesign(i)];
    end
    
    % Map with as high as Urgency signal of 5 still works well 
    [tempOutcomeParams_wTheta, errorLL] = fminsearchbnd(@(params0_wTheta) Parallel_Fitting(monkeyData, k, totalTime, deltaT, bound, whichSolve, extraParamsName, params0_wTheta), params0_wTheta, LB, UB, options);
    
    % Save some parameters (Error and outcome parametersof that error)
    allErrorLL_wTheta(i) = errorLL; 
    outcomeParamsMat_wTheta(i,:) = tempOutcomeParams_wTheta;
    cd ../results
    save(['genji_Fitting_k=' num2str(k) '_Joint_' char(datetime("today")) '.mat'], 'outcomeParamsMat_wTheta', 'allErrorLL_wTheta', 'gridDesign');
    cd ../code

    % Display time information and iteration value
    display(['Iteration Number =' num2str(i)])
    datetime 
end

[outcomeError_wTheta, b] = min(allErrorLL_wTheta);
outcomeParams = outcomeParamsMat_wTheta(b, :);

%% Now that we have fits, use those to fit the data. 

[totalError, logOddMaps_FromFit] = Parallel_Fitting(monkeyData, k, totalTime, deltaT, bound, whichSolve, extraParamsName, outcomeParams);

% Equations needed
%twoDGauss = @(x, s, mu, sigma, t) 1./(2.*pi.*sqrt(det(t.*sigma))) .* exp(-1./(2).* sum(((x - s - mu.*t)/(sigma.*t)).*(x - s - mu.*t), 2)); 
sj_even = @(j, alpha, k, s0) (1/sin(pi/k))*([sin(j*alpha + pi/k)    sin(j*alpha);           -sin(j*alpha)           -sin(j*alpha - pi/k)])  * s0';
sj_odd = @(j, alpha, k, s0) (1/sin(pi/k))* ([sin(j*alpha)           sin(j*alpha - pi/k);    -sin(j*alpha + pi/k)    -sin(j*alpha)])         * s0';
alpha =@(k) (k-1)/k * pi;
weightj =@(j,mu,sigma,sj,s0) (-1)^j * exp(mu/(sigma) * (sj-s0)');

% Inputs: Conditions (Working in 3rd quadrant so the entire grid is negative -x,-y)
k = k; % parameter for correlation and number of images needed
rho = -cos(pi/k); %correlation
sigma=[1 rho;rho 1]; % covariance matrix (You can plug in correlations for CovXY when Variance is 1 for X and Y)
senK = outcomeParams(1); %10; %coherence and sensitivity
s0 = [outcomeParams(2), outcomeParams(2)]; %Start of propagation (Notice negative; this can change depending on your need for distance to upper limit/bound of 0)
deltaT = .01; startT = 0; endT = 3; %Time is in seconds % delta normally .001
deltaX = .025; deltaY = .025; % These should always be the same
nonDT_std = .034; % Deduce using van den Berg et al., 2016

nonDT = outcomeParams(3);
theta = outcomeParams(4);
timeTheta = 0; %inf; %outcomeParams(6); %0 or inf; %round(outcomeParams(5));

% Urgency add on
if length(outcomeParams) > 4
    urgencyMax = outcomeParams(5); %*1000;
else
    urgencyMax = 0;
end

% Asymmetrical Non-DT
if any(find(contains(extraParamsName, 'asymmetrical')))
    indexAdd = find(contains(extraParamsName, 'asymmetrical'));
    nonDT_Other = outcomeParams(4+indexAdd);
else
    nonDT_Other = nonDT;
end
% P high
if any(find(contains(extraParamsName, 'wagerOffset')))
   %bias p high term (for hanzo)
   indexAdd = find(contains(extraParamsName, 'wagerOffset'));
   biasPHigh = outcomeParams(4+indexAdd); %.05;
else
   biasPHigh = 0;
end
% Urgency half life
if any(find(contains(extraParamsName, 'useUrgencyHalfMax')))
   %bias p high term (for hanzo)
   indexAdd = find(contains(extraParamsName, 'useUrgencyHalfMax'));
   urgencyTauHalf = outcomeParams(4+indexAdd); %.05;
else
   urgencyTauHalf = 0;
end

multiplyBy = 3;
s0_fitted = round(s0/deltaX) * deltaX;

yi = 0:-deltaY:s0_fitted(1)*multiplyBy; % X-Mesh Grid
xi = s0_fitted(2)*multiplyBy:deltaX:0; % Y-Mesh Grid
ti = startT:deltaT:endT; %startT:deltaT:endT; % Time Grid (should always go from 0 -> T)
[x, y] = meshgrid(xi, yi); %make mesh
subInd = [x(:) y(:)]; %make 2D points
Alpha = alpha(k);
reflectingBound = 1;

% For choice do we calculate a Fokker Planck Map for all coherences of
% interest
modelCoherences = -.52:.02:.52; % unique(monkeyData.cohs); % -.52:.02:.52; %try to make this even (this does not have to range from negative to positive since its a symmetrical at 0, 0:.52)
% modelCoherences(modelCoherences == 0) = -eps;
% modelCoherences = [modelCoherences(1:27) eps modelCoherences(28:end)];

totalPassBoundProb = nan(size(modelCoherences));
averageRT_model = nan(size(modelCoherences)); averageRT_High_model = nan(size(modelCoherences)); averageRT_Low_model = nan(size(modelCoherences)); 
probRight_model = nan(1, length(modelCoherences)); probRight_High_model = nan(1, length(modelCoherences)); probRight_Low_model = nan(1, length(modelCoherences));
probCorrectxTime_model = nan(3/deltaT, length(modelCoherences));
probRight_model_2 = nan(1, length(modelCoherences));
 
probHigh_model = nan(1, length(modelCoherences)); probHigh_model_Temp = nan(1, length(modelCoherences)); 
probHighxTime_model = nan(3/deltaT, length(modelCoherences));
probHigh_Correct_model = nan(1, length(modelCoherences)); probHigh_Correct_model_Temp = nan(1, length(modelCoherences));
probHigh_Incorrect_model = nan(1, length(modelCoherences)); probHigh_Incorrect_model_Temp = nan(1, length(modelCoherences));
probHigh_Correct_Error = nan(1, length(modelCoherences));
probHigh_Incorrect_Error = nan(1, length(modelCoherences));
probHigh_temp = nan(1, length(modelCoherences));

if 1
    tic
    for c=1:length(modelCoherences) % (c = 1:length(modelCoherences), 8) 
    
        coh = modelCoherences(c);  %coherence and sensitivity
        mu = [senK*coh, -senK*coh]; %Drift rate  
        
        [xytFPMat_FromC] = FP_Drugowitsch(subInd, deltaX, length(yi)*length(xi), mu, k, ti(end), deltaT, s0_fitted, urgencyMax, urgencyTauHalf); %you might then have to multiply by deltaT to actually get the p(x,y,t|parameters)
        % Make all values below eps -> eps
        xytFPMat_FromC(xytFPMat_FromC<eps) = eps;
        xytFPCell = reshape(xytFPMat_FromC, [length(xi), length(yi), size(xytFPMat_FromC,2)]); %re-arrange 

        % Make it a reflecting bound
        if reflectingBound == 1
%             tempXYT = xytFPCell;
            % Flip at MidPoint
            pointOfReflection = size(xytFPCell);
            if pointOfReflection(1) == pointOfReflection(2)
                pointOfReflection = floor(pointOfReflection(1)/2); %priority takes X-axis
            else
                pointOfReflection = floor(min(pointOfReflection(1:2))/2);
            end
            % Do the geometry needed to reflect on the bounds
            horizontalExcess = xytFPCell(:,1:pointOfReflection,:);
            verticalExcess = xytFPCell((pointOfReflection+1+rem(size(xytFPCell,1), 2)):end,:,:); % Must add plus 1!
            % Loop through time now
            for tt = 2:length(pointOfReflection)
                % Reflect Horizontal part
                horizontalExcess(:,:,tt) = fliplr(horizontalExcess(:,:,tt));
                xytFPCell(:,(pointOfReflection+1):pointOfReflection*2, tt) = xytFPCell(:,(pointOfReflection+1):pointOfReflection*2, tt) + horizontalExcess(:,:,tt);
                % Reflect Horizontal part
                verticalExcess(:,:,tt) = flipud(verticalExcess(:,:,tt));
                xytFPCell(1+rem(size(xytFPCell,1), 2):(pointOfReflection+rem(size(xytFPCell,1), 2)), :, tt) = xytFPCell(1+rem(size(xytFPCell,1), 2):(pointOfReflection+rem(size(xytFPCell,1),2)), :, tt) + verticalExcess(:,:,tt);
            end
            % After everything remove the segments that were reflected
            xytFPCell(:,1:pointOfReflection,:) = [];
            xytFPCell((pointOfReflection+1+rem(size(xytFPCell,1), 2)):end,:,:) = [];
%             xytFPCell = tempXYT;
        end

        %Calculate Probability of Betting High (logOddsMap)
        atBound1 = flipud(squeeze(xytFPCell(yi==bound,:, 2:end))); 
        if reflectingBound == 0
            atBound2 = squeeze(xytFPCell(: ,xi == bound, 2:end)); 
        else
            atBound2 = squeeze(xytFPCell(: ,end-1, 2:end)); 
        end
    
        % Calulcate threshold values 
        [indexHighWager] = calculating_MoreFlexibleTheta(logOddMaps_FromFit, theta, timeTheta); % Allows for flat cutoff and changes in timePoint 
        indexHighWager = indexHighWager(:, 2:end);
        
        % Any Lapses to P High (Only for Hanzo)
        % Probability of High should be the probabilities that higher than theta when either bound is hit 
        probHigh_model(c) = (nansum(atBound1(indexHighWager), 'all') + nansum(atBound2(indexHighWager), 'all'))/ (nansum(atBound1, 'all') + nansum(atBound2, 'all')); %Adding a bias term (.05)

        % Also Calculate P(Confidence=High|{Correct, Incorrect})
        if modelCoherences(c) < 0
            probHigh_Correct_model(c) = nansum(atBound1(indexHighWager), 'all')/ nansum(atBound1, 'all') - biasPHigh;
            probHigh_Incorrect_model(c) = nansum(atBound2(indexHighWager), 'all')/ nansum(atBound2, 'all') - biasPHigh;
        else
            probHigh_Correct_model(c) = nansum(atBound2(indexHighWager), 'all')/ nansum(atBound2, 'all') - biasPHigh;
            probHigh_Incorrect_model(c) = nansum(atBound1(indexHighWager), 'all')/ nansum(atBound1, 'all') - biasPHigh;
        end
        % IF YOU CAN COMPUTE THE LOG LIKEHOOD FOR P(HIGH|CORRECT) AND
        % P(HIGH|INCORRECT)

        % Calculate Confidence Through Time (This is wrong, you need to
        % multiply by the P(Correct={0,1}|Coh. Actually, is this true? I think both equations end up being equivalent)
        % This is true for probHigh_model
        probHighxTime_model(:, c) = (sum(atBound1 .* indexHighWager,1, 'omitnan') + sum(atBound2 .* indexHighWager,1, 'omitnan'))./ ...
            (sum(atBound1,1, 'omitnan') + sum(atBound2,1, 'omitnan'));

        % Calculate Average RT
        survival = squeeze(sum(sum(xytFPCell))); %sum(xytFPMat_FromC);
        flux = diff(1-survival);
        normFlux = flux./nansum(flux);
        normFlux(1) = eps;
        if modelCoherences(c) > 0 
            totalTimeFlux = conv(normFlux, normpdf(ti(2:end), nonDT, nonDT_std)); %include the non-decision time
        elseif modelCoherences(c) <= 0 
            totalTimeFlux = conv(normFlux, normpdf(ti(2:end), nonDT_Other, nonDT_std)); %include the non-decision time
        end

        % Normalize
        totalTimeFlux = totalTimeFlux(1:(3/deltaT)) / sum(totalTimeFlux(1:(3/deltaT)));
        averageRT_model(c) = sum(totalTimeFlux(1:(3/deltaT)) .* (1:(3/deltaT)));
    
        % Calculate prob Right
        probRight_model(c) = sum(atBound2, 'all') ./ (sum(atBound1, 'all') + sum(atBound2, 'all'));

        % Calculate Prob{High|Right) vs Prob(Low|Right)
        probHigh_Right = nansum(atBound2(indexHighWager), 'all') / nansum(atBound2, 'all'); %THIS COULD BE WRONG, IT COULD BE BOUND2 AND NOT BOUND1
        probHigh_Left = nansum(atBound1(indexHighWager), 'all') / nansum(atBound1, 'all'); 
        % P[Right|High]
        probRight_High_model(c) = probHigh_Right * probRight_model(c) / probHigh_model(c);
        probRight_Low_model(c) = (1-probHigh_Right) * probRight_model(c) / (1-probHigh_model(c));

        % Calculate Correct Through Time
        if modelCoherences(c) >= 0 
            probCorrectxTime_model(:, c) = sum(atBound2, 1, 'omitnan') ./ (sum(atBound1, 1, 'omitnan') + sum(atBound2, 1, 'omitnan'));
        else
            probCorrectxTime_model(:, c) = sum(atBound1, 1, 'omitnan') ./ (sum(atBound1, 1, 'omitnan') + sum(atBound2, 1, 'omitnan'));
        end

        % Now try solving for E[P(RT|Wager)] 
        RTdist = flux(2:end)./sum(flux(2:end));
        
        if modelCoherences(c) > 0 
            totalTimeFlux = conv(RTdist, normpdf([deltaT:deltaT:ti(end)], nonDT, nonDT_std)); %include the non-decision time
        elseif modelCoherences(c) <= 0 
            totalTimeFlux = conv(RTdist, normpdf([deltaT:deltaT:ti(end)], nonDT_Other, nonDT_std)); %include the non-decision time
        end
        
        probHigh_Time_model = (nansum(atBound1 .* indexHighWager, 1) + nansum(atBound2 .* indexHighWager, 1)) ./ (nansum(atBound1, 1) + nansum(atBound2,1)); %distribution of P(high wager|RT)
        probHigh_Time_model(1:2) = 1;
        
        distribution_RThigh_model = (probHigh_Time_model(2:end) .* totalTimeFlux(2:ti(end)*(1/deltaT))) ./ (probHigh_model(c)+eps); % bayes for p(RT|wager==1)
        distribution_RTlow_model = ((1-probHigh_Time_model(2:end)) .* totalTimeFlux(2:ti(end)*(1/deltaT))) ./ (1-probHigh_model(c) + eps); % bayes for p(RT|wager==0)
        
        averageRT_High_model(c) = sum( distribution_RThigh_model .* ti(3:end) ./ sum(distribution_RThigh_model) );
        averageRT_Low_model(c) = sum( distribution_RTlow_model .* ti(3:end) ./ sum(distribution_RTlow_model) );

        % Lastly, update PHigh based on the bias term if it exist
        probHigh_model(c) = probHigh_model(c) - biasPHigh;
    end 
end


%% Plot fits with actual data
[avgChoiceData, ~, stdChoiceData]   = behavioralAverages(monkeyData.choice, monkeyData.cohs); % choice
[avgHighChoiceData, ~, stdHighChoiceData]   = behavioralAverages(monkeyData.choice(monkeyData.wager==1), monkeyData.cohs(monkeyData.wager==1)); % choice
[avgLowChoiceData, ~, stdLowChoiceData]   = behavioralAverages(monkeyData.choice(monkeyData.wager==0), monkeyData.cohs(monkeyData.wager==0)); % choice

[avgRTData, ~, stdRTData]           = behavioralAverages(monkeyData.RT, monkeyData.cohs); % RT
[avgRTHighData, ~, stdRTHighData]    = behavioralAverages(monkeyData.RT(monkeyData.wager==1), monkeyData.cohs(monkeyData.wager==1)); % RT High
[avgRTLowData, ~, stdRTLowData]      = behavioralAverages(monkeyData.RT(monkeyData.wager==0), monkeyData.cohs(monkeyData.wager==0)); % RT Low

[avgWagerData, ~, stdWagerData]             = behavioralAverages(monkeyData.wager, monkeyData.cohs); % Wager
[avgWagerCorrData, ~, stdWagerCorrData]     = behavioralAverages(monkeyData.wager(monkeyData.correct==1), monkeyData.cohs(monkeyData.correct==1)); % wager | Correct
[avgWagerIncorrData, ~, stdWagerIncorrData] = behavioralAverages(monkeyData.wager(monkeyData.correct==0), monkeyData.cohs(monkeyData.correct==0)); % wager | incorrect

f = figure;
f.Position = [0 500 1400 800];

subplot(2,3,1)
uniqMonkeyCoh = unique(monkeyData.cohs);
errorbar(uniqMonkeyCoh, avgChoiceData, stdChoiceData, 'k.',  'Markersize', 25); hold on;
plot(modelCoherences , probRight_model, '-', 'Color', [.5 .5 .5], 'Linewidth', 2.5); hold off;
% labels
% ylabel('Proportion rightward choice')
% xlabel('Motion strength (%coh)')
xticks([-.512 -.256 0 .256 .512])
% xticklabels({'', '', '', '', ''})
xticklabels({'-51.2', '-25.6', '0', '25.6', '51.2'})
xlim([-.512 .512]);
% y-label
ylabel('Proportion rightward choices')
% Pretty graphs
% prettifyGraphs(20, 'normal', 1)

subplot(2,3,2);
errorbar(uniqMonkeyCoh, avgRTData./1000, stdRTData./1000, 'k.',  'Markersize', 25); hold on; %Simulation calculations
% plot([-flip(modelCoherences(2:end)) modelCoherences], [flip(averageRT_model(1:end-1)) averageRT_model], 'ks--', 'Linewidth', 2); hold off;
if length(outcomeParams)>6
    plot(modelCoherences(1:27), averageRT_model(1:27)./(1/deltaT), '-', 'Color', [.5 .5 .5], 'Linewidth', 2.5);
    plot(modelCoherences(28:end), averageRT_model(28:end)./(1/deltaT), '-', 'Color', [.5 .5 .5], 'Linewidth', 2.5); hold off
else
    plot(modelCoherences, averageRT_model./(1/deltaT), '-', 'Color', [.5 .5 .5], 'Linewidth', 2.5); hold off;
end
% labels
% ylabel('Reaction time (s)')
% xlabel('Motion strength (%coh)')
xticks([-.512 -.256 0 .256 .512])
% xticklabels({'', '', '', '', ''})
xticklabels({'-51.2', '-25.6', '0', '25.6', '51.2'})
xlim([-.512 .512]);
xlabel('Motion strength (% coh)')
% y-label
ylabel('Reaction time (s)')
% Pretty graphs
% prettifyGraphs(20, 'normal', 1)

subplot(2,3,3);
errorbar(uniqMonkeyCoh, avgWagerData, stdWagerData, 'k.',  'Markersize', 25); hold on; %Simulation calculations
% plot([flip(-modelCoherences(2:end)) modelCoherences], [flip(probHigh_model(2:end)) probHigh_model], 'k--', 'Linewidth', 2); hold off;  %Model fit
plot(modelCoherences, probHigh_model, '-', 'Color', [.5 .5 .5], 'Linewidth', 2.5); hold off;  %Model fit
% labels
% ylabel('Proportion high bet')
% xlabel('Motion strength (%coh)')
xticks([-.512 -.256 0 .256 .512])
%xticklabels({'', '', '', '', ''})
xticklabels({'-51.2', '-25.6', '0', '25.6', '51.2'})
xlim([-.512 .512]);
ylim([.5 1])
% y-label
ylabel('Proportion high bets')
% Pretty graphs
% prettifyGraphs(20, 'normal', 1)


%% Plot Conditional Results
subplot(2,3,5)
errorbar(unique(monkeyData.cohs), avgRTHighData./1000, stdRTHighData./1000, 'r.', 'Markersize', 25); hold on;
errorbar(unique(monkeyData.cohs), avgRTLowData./1000, stdRTLowData./1000, 'b.', 'Markersize', 25); 
if length(outcomeParams)>6
    plot(modelCoherences(1:27), averageRT_High_model(1:27), '-', 'Color', [255 127 127]./255, 'Linewidth', 2.5);
    plot(modelCoherences(28:end), averageRT_High_model(28:end), '-', 'Color', [255 127 127]./255, 'Linewidth', 2.5); 
    plot(modelCoherences(1:27) , averageRT_Low_model(1:27), '-', 'Color', [127 127 255]./255, 'Linewidth', 2.5);
    plot(modelCoherences(28:end) , averageRT_Low_model(28:end), '-', 'Color', [127 127 255]./255, 'Linewidth', 2.5); hold off
else
    plot(modelCoherences, averageRT_High_model, 'r-', 'Linewidth', 2.5); 
    plot(modelCoherences , averageRT_Low_model, 'b-', 'Linewidth', 2.5); hold off;
end
hold off;
ylabel('Reaction time (s)')
xlabel('Motion strength (% coh)')
xticks([-.512 -.256 0 .256 .512])
xticklabels({'-51.2', '-25.6', '0', '25.6', '51.2'})
xlim([-.512 .512]);
ylim([.3 .8])
% Pretty graphs
% prettifyGraphs(20, 'normal', 1)


subplot(2,3,4)
errorbar(uniqMonkeyCoh, avgHighChoiceData, stdHighChoiceData, 'r.', 'Markersize', 25); hold on;
errorbar(uniqMonkeyCoh, avgLowChoiceData, stdLowChoiceData, 'b.', 'Markersize', 25);
plot(modelCoherences , probRight_High_model, '-', 'Color', [255 127 127]./255, 'Linewidth', 2.5); 
plot(modelCoherences , probRight_Low_model, '-', 'Color', [127 127 255]./255, 'Linewidth', 2.5); hold off;
% labels
ylabel('Proportion rightward choices')
xticks([-.512 -.256 0 .256 .512])
xticklabels({'-51.2', '-25.6', '0', '25.6', '51.2'})
xlim([-.512 .512]);
ylim([0 1])
legend('high bet', 'low bet', 'Location', 'Northwest')
legend box off
% Pretty graphs
% prettifyGraphs(20, 'normal', 1)

subplot(2,3,6)
errorbar(uniqMonkeyCoh, avgWagerCorrData, stdWagerCorrData, '.', 'Color', [0.4940 0.1840 0.5560], 'Markersize', 25); hold on;
errorbar(uniqMonkeyCoh(2:end-1), avgWagerIncorrData(2:end-1), stdWagerIncorrData(2:end-1), '.', 'Color', [0.4660 0.6740 0.1880], 'Markersize', 25);
plot(modelCoherences , probHigh_Correct_model, '-', 'Color', [0.4940 0.1840 0.5560], 'Linewidth', 2.5); 
plot(modelCoherences , probHigh_Incorrect_model, '-', 'Color', [0.4660 0.6740 0.1880], 'Linewidth', 2.5); hold off;
ylabel('Proportion high bets')
xticks([-.512 -.256 0 .256 .512])
xticklabels({'-51.2', '-25.6', '0', '25.6', '51.2'})
xlim([-.512 .512]);
ylim([.5 1])
legend('correct', 'incorrect', 'Location', 'Southwest')
legend box off
% Pretty graphs
% prettifyGraphs(20, 'normal', 1)

%% Comparing Models Using Both BIC and AIC
if 1
    indexModel = 0;
    aicVector = nan(1,10);
    bicVector = nan(1,10);
    nameOfModel = cell(1,10);
end

% AIC & BIC
[aic, bic] = aicbic(-outcomeError_wTheta , length(outcomeParams)+0, length(monkeyData.choice)); 

indexModel = indexModel+1;
aicVector(indexModel) = aic;
bicVector(indexModel) = bic;
nameOfModel{indexModel} = 'Genji 2D_Accumulator';

%% Last thing to do is compute the Log-Likelihood error for P(High|Correct) and P(High|Incorrect)
% Only works if modelCoherences == monkeydata.coh
errorWagerHigh = nan(1,length(modelCoherences));
errorWagerLow = nan(1,length(modelCoherences));

for i = 1:length(modelCoherences)
    % Find the empirical data
    empData_High  = sum(monkeyData.wager == 1 & monkeyData.correct == 1 & monkeyData.cohs == modelCoherences(i));
    empData_Low  = sum(monkeyData.wager == 0 & monkeyData.correct == 1 & monkeyData.cohs == modelCoherences(i));
    errorWagerHigh(i) = empData_High*log(probHigh_Correct_model(i)) + empData_Low * log(1-probHigh_Correct_model(i));
    % Low
    empData_High  = sum(monkeyData.wager == 1 & monkeyData.correct == 0 & monkeyData.cohs == modelCoherences(i));
    empData_Low  = sum(monkeyData.wager == 0 & monkeyData.correct == 0 & monkeyData.cohs == modelCoherences(i));
    errorWagerLow(i) = empData_High*log(probHigh_Incorrect_model(i)) + empData_Low * log(1-probHigh_Incorrect_model(i));

end

% Final values
errorWagerHigh = -sum(errorWagerHigh);
errorWagerLow = -sum(errorWagerLow);


%% Save the figure and workspace just in case
cd ../results
saveas(gcf,'data_wFit_Parallel','fig')
clear f
save workspace_Parallel.mat
cd ../code
