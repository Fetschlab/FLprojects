function [totalError, logOddsMap, posteriorProb_LatentState_GivenChoice, firstRaceLoss] = Parallel_Fitting(data, numOfImages, endTi, deltaT, bound, whichSolve, extraParamsName, params)% Bl, Bu, T0, sigma0, sensitivity)
% Description: Fit choice, RT, and wager data using FP_Drugowitsch.c,
% which is based on:

% Shan, Moreno-Bote, and Drugowitsch (2019) Family of closed-form solutions
% for two-dimensional correlated diffusion processes. Phys Rev E100, 032132

% Also useful for fitting any collapsing bound model. This type of model
% requires finding an error that includes not just RT but also Wager.

% Code provided for peer review of:
% M. Vivar-Lazo & C.R. Fetsch, Neural basis of concurrent deliberation
% toward a choice and degree of confidence (2024)

% run using FP_Fitting_Wrapper_NN

%%
% CF: for legacy reasons the data are expected to be row vectors, otherwise 
% the error calculations (which use elementwise multiplication) give n x n
% matrices instead of vectors, which then do not sum to a scalar.
if iscolumn(data.choice)
    data.direction = data.direction';
    data.coherence = data.coherence';
    data.choice = data.choice';
    data.RT = data.RT';
    data.correct = data.correct';
    data.wager = data.wager';
    data.cohs = data.cohs'; 
end


%% Parameters of interest (only 1):

% whichSolve: 1--all covariates independent, 2--RT independent but Choice &
% Wager multinomial, 3--all together as one big joint probabiity

% Get Parameters of Interest
sensitivity = params(1); % Sensitivty K value
s0 = [params(2), params(2)]; % Starting point (assuming symmetry)
s0_1D = params(2);
mu_nonDT = params(3); % Non-decision time average
var_nonDT = .034^2; % Non-decision time variance (For Genji = .023 for Hanzo = .032)
if length(params) > 3
    theta = params(4); % wager cutoff
    timeTheta = 0; %0 inf; % time that moderates the theta cutoff when constructing LogOdds Map
    urgencyTauHalf = 0;
    % Next thing to try is having a time dependent theta (linear or
    % non-linear -> up to you)
end
% Load extra params. If needed these params will change based on needed
% from wrapper
urgencyMax = 0;
mu_nonDT_null = mu_nonDT;
biasWager = 0;

if ~isempty(extraParamsName)
    extraParams = 0;
    for i = 1:length(extraParamsName)
        if strcmp(extraParamsName{i}, 'urgencySignal')
            extraParams = extraParams + 1;
            urgencyMax = params(4+extraParams); % * 1000; % Highest value the urgency signal can take (theres a good reason to multiply by 1000)
        elseif strcmp(extraParamsName{i}, 'asymmetrical')
            extraParams = extraParams + 1;
            mu_nonDT_null = params(4+extraParams);
        elseif strcmp(extraParamsName{i}, 'wagerOffset')
            extraParams = extraParams + 1;
            biasWager = params(4+extraParams);
        elseif strcmp(extraParamsName{i}, 'thetaTimeCutoff')
            extraParams = extraParams + 1;
            timeTheta = params(4+extraParams); %round(params(4+extraParams));
        elseif strcmp(extraParamsName{i}, 'useUrgencyHalfMax')
            extraParams = extraParams + 1;
            urgencyTauHalf = params(4+extraParams);
        end
    end
end

% Any variables needed
k = numOfImages; 
% bound = -.5; %temporary for now
reflectingBound = 1; % is it reflecting bound?

% Set up the mesh
deltaX = .025; deltaY = .025; % These should always be the same
deltaT = deltaT; startT = 0; endT = endTi; %Time is in seconds (normally delta is .001)
multiplyGrid = 3;
% Formulation
% You have to double check to make sure the starting point is divisible by
% the deltaX/deltaY
if rem(s0_1D, deltaX) ~= 0 
    s0_1D = deltaX * round(s0_1D/deltaX);
end
yi = 0:-deltaY:s0_1D*multiplyGrid; % X-Mesh Grid
xi = s0_1D*multiplyGrid:deltaX:0; % Y-Mesh Grid
ti = startT:deltaT:endT; %startT:deltaT:endT; % Time Grid (should always go from 0 -> T)
[x, y] = meshgrid(xi, yi); %make mesh
subInd = [x(:) y(:)]; %make 2D points
s0 = [s0_1D, s0_1D];

%modelCoherences =  [-.512 -.256 -.128 -.064 -.032 0 .032 .064 .128 .256 .512]; 
modelCoherences = unique(data.cohs);

% Error varialbes
errorRT = nan(1, length(modelCoherences));
errorRT_High = nan(1, length(modelCoherences));
errorRT_Low = nan(1, length(modelCoherences));
errorWager = nan(1, length(modelCoherences));

probRight_Coh = nan(1,length(modelCoherences));

%

if exist('theta', 'var')
    % Variables to make LogOddsMap
    if reflectingBound == 0
        correctFinalFPMat = zeros(length(xi), length(ti));
        errorFinalFPMat = zeros(length(yi), length(ti));
    else
        correctFinalFPMat = zeros(ceil(length(xi)/2), length(ti));
        errorFinalFPMat = zeros(ceil(length(yi)/2), length(ti));
    end

    % Number of Parallels
    if length(modelCoherences) < 8
        numOfParallels = length(modelCoherences);
    else 
        numOfParallels = 8;
    end
    % Algorthm run
%     parfor (c = 1:length(modelCoherences), numOfParallels) %reduce the number of workers because of memory usage (macbook has 10)
    for c = 1:length(modelCoherences)
        coh = modelCoherences(c);  %coherence and sensitivity
        mu = [sensitivity*coh, -sensitivity*coh]; %Drift rate
        % You need pdf to calculate things correctly (Wager)
        % Run Drugowitsch Method
        [xytFPMat_FromC] = FP_Drugowitsch(subInd, deltaX, length(yi)*length(xi), mu, k, endTi, deltaT, s0, urgencyMax, urgencyTauHalf); %you might then have to multiply by deltaT to actually get the p(x,y,t|parameters)
        
        % Make all values below eps -> eps
        xytFPMat_FromC(xytFPMat_FromC<eps) = eps;
        xytFPCell = reshape(xytFPMat_FromC, [length(xi), length(yi), size(xytFPMat_FromC,2)]); % .* 1/3000; %re-arrange 
        % Construct 1D Propogation for Time-x-X and Time-x-Y (For Choice and
        % RT)

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
            verticalExcess = xytFPCell((pointOfReflection+1+rem(size(xytFPCell,1), 2)):end,:,:); % Must add plus 1 for even xmesh but 2 for odd xmesh!
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
    
        % Calculate how much has passed and flux
        survival = squeeze(sum(sum(xytFPCell)));
        flux = diff(survival);
        normalizeFlux = flux./sum(flux, 'omitnan'); % normalize
        normalizeFlux(1) = eps;
        % Save flux incases needed for later
        flux(1) = eps;
        laterFlux{c} = flux; %normalizeFlux
    
        % Now that you have probability flux, calculate RT average (Add
        % non-decision time)
        
        if modelCoherences(c) > 0 
            totalTimeFlux = conv(normalizeFlux, normpdf([deltaT:deltaT:endTi], mu_nonDT, sqrt(var_nonDT))); %include the non-decision time (var_nonDT)
        elseif modelCoherences(c) <= 0 
            totalTimeFlux = conv(normalizeFlux, normpdf([deltaT:deltaT:endTi], mu_nonDT_null, sqrt(var_nonDT))); %include the non-decision time (var_nonDT_null)
        end
        
        % Normalize totalTimeFlux
        % totalTimeFlux = totalTimeFlux(1:3000) / sum(totalTimeFlux(1:3000));
        totalTimeFlux = totalTimeFlux/sum(totalTimeFlux);
        totalTimeFlux(totalTimeFlux > 1-eps) = 1-eps; 
        totalTimeFlux(totalTimeFlux<eps) = eps; 
        % Calculate Estimate Error 
        tempRT = data.RT(data.cohs == modelCoherences(c));
        if deltaT == .01
            tempRT = round(tempRT(~isnan(tempRT))/10);
        elseif deltaT == .025
            tempRT = round(tempRT(~isnan(tempRT))/25);
        elseif deltaT == .001
            tempRT = round(tempRT(~isnan(tempRT)));
        else
            warning('Your time specificity is not supported')
        end        
        errorRT(c) = sum(log(totalTimeFlux(tempRT)), 'omitnan');
    
        % Construct 1D propogation for X and Y at bound (For calculation of
        % Wager Error)
        firstRaceLoss{c} = flipud(squeeze(xytFPCell(yi == bound, :, :))); %Find the p(vR1, t | vR2==Bound, C) (For Left Choices)
        if reflectingBound == 0
            secondRaceLoss{c} = squeeze(xytFPCell(:, xi == bound, :)); %Find the p(vR2, t | vR1==Bound, C) (For Right Choices)
        else
            secondRaceLoss{c} = squeeze(xytFPCell(:, end-1, :));
        end
        % Place Boundaries in Correct or Incorrect Matrix
        if modelCoherences(c) > 0 %Positive Coherence
            correctFinalFPMat = correctFinalFPMat + secondRaceLoss{c}  .* (1/(length(modelCoherences)+1)); % Integral( p(DV|coh) * p(coh) dCoh)
            errorFinalFPMat = errorFinalFPMat + firstRaceLoss{c} .* (1/(length(modelCoherences)+1));
        elseif modelCoherences(c) < 0 % negative and 0 coherence
            correctFinalFPMat = correctFinalFPMat + firstRaceLoss{c}  .* (1/(length(modelCoherences)+1)); % Integral( p(DV|coh) * p(coh) dCoh)
            errorFinalFPMat = errorFinalFPMat + secondRaceLoss{c} .* (1/(length(modelCoherences)+1));
        else
            correctFinalFPMat = correctFinalFPMat + (firstRaceLoss{c} + secondRaceLoss{c})  .* (1/(length(modelCoherences)+1)); % Integral( p(DV|coh) * p(coh) dCoh)
            errorFinalFPMat = errorFinalFPMat + (firstRaceLoss{c} + secondRaceLoss{c}) .* (1/(length(modelCoherences)+1));
        end

        % Also save the amount of probabilites that will exit each time
        % point (Values at the bound at every t)
        leftChoice = sum(squeeze(xytFPCell(yi == bound, :, :)), 'all', 'omitnan');
        if reflectingBound == 0
            rightChoice = sum(xytFPCell(:, xi == bound, :), 'all', 'omitnan');
        else
            rightChoice = sum(xytFPCell(:, end-1, :), 'all', 'omitnan');
        end
        % You should in theyr have to do P(:,xi=b,t) * P(RT=t) / (P(:,xi=b,t) * P(RT=t) + P(yi=b,:,t) * P(RT=t)) 
        % but no need cause P(RT) cancels out from top and bottom
        probRight_Coh(c) = rightChoice ./ (rightChoice + leftChoice);
    end

    % Now find the LogOddsMap 
    logOddsMap = log(correctFinalFPMat ./ errorFinalFPMat);
    [indexHighWager] = calculating_MoreFlexibleTheta(logOddsMap, theta, timeTheta); % Allows for flat cutoff and changes in timePoint 
    % This is a work in progress, calculating the Posterior in correct trials
    posteriorProb_LatentState_GivenChoice = correctFinalFPMat./(correctFinalFPMat + errorFinalFPMat); % Part of Equation in Kiani 2009 Science which is equivalent to Figure 1 equation for Confidence in Pouget, Drugowitsch, and Kepecs 2016  

    % After having P(Right|Coherence) & P(RT|Coherence) decide whether
    % you're going to find the likelihood assuming independence between
    % Choice and Wager or if there is some relationship
    if whichSolve == 1 % Independent Choice and Wager
        
        % Now find P(High|Coherence)
        probHigh_model = nan(1, length(modelCoherences));
         % Probability of High should be the probabilities that higher than theta when either bound is hit 
        for c = 1:length(modelCoherences)
            % Calculate prob high (Is this right?). Shouldn't we be
            % calculating P(High|Acc=Corr, Coh=C)*P(Acc=Corr|Coh=C) + P(High|Acc=Error, Coh=C)*P(Acc=Error|Coh=C)
            % If so then this is weighted incorrectly. When you fix this,
            % also fix it for the rest of the fitting procedures and the
            % portion of this in the wrapper. Its equivalent.
            probHigh_model(c) = (nansum(firstRaceLoss{c}(indexHighWager), 'all') + nansum(secondRaceLoss{c}(indexHighWager), 'all'))/ (nansum(firstRaceLoss{c}, 'all') + nansum(secondRaceLoss{c}, 'all')) - biasWager;
        end

        % Find the logLL 
        % Avoid making -inf or inf with Log
        probHigh_model(probHigh_model<eps) = eps;
        probHigh_model(probHigh_model>1-eps) = 1-eps;
        % Have everytrial take a prob of making its wager based on coherence
        % Eventually you can just use the symmetry to save run time
        probHigh_trial = nan(1, length(modelCoherences));
        for z = 1:length(modelCoherences)
            probHigh_trial(data.cohs == modelCoherences(z)) = probHigh_model(z); 
        end

        % Binary Wager Outcome 
        % Also fit for conditional reaction times (High vs Low)
        errorWager = sum((data.wager == 1) .* log(probHigh_trial) + (data.wager ==0) .* log(1-probHigh_trial), 'omitnan');

        % calculate the Probability Right
        probRight_trial = nan(1, length(modelCoherences));
        % Calculate error
        for z = 1:length(modelCoherences)
            probRight_trial(data.cohs == modelCoherences(z)) = probRight_Coh(z); 
        end

        % LogLL for Choice
        errorChoice = sum((data.choice == 1) .* log(probRight_trial) + (data.choice ==0) .* log(1-probRight_trial), 'omitnan');
    
        % total error (Binary outcomes)
        totalError = -sum(errorRT) + -errorWager + -errorChoice;

    
    elseif whichSolve == 2 %RT independent and Choice,Wager multinomial
        
        % RUN INSTEAD OF P(HIGH|RIGHT)*P(RIGHT) DO P(RIGHT|HIGH)*P(HIGH)

        % Now find P(High|Coherence)
        ProbHigh_ConditionedRightChoice_model = nan(1, length(modelCoherences));
        ProbHigh_ConditionedLeftChoice_model = nan(1, length(modelCoherences));
        % Probability of High should be the probabilities that higher than theta when either bound is hit 
        for c = 1:length(modelCoherences)
            % [All the Probs above the theta @ Bound]/[All Probs @ Bound]
            ProbHigh_ConditionedRightChoice_model(c) = (nansum(secondRaceLoss{c}(indexHighWager), 'all'))/ (nansum(secondRaceLoss{c}, 'all')) - biasWager;
            ProbHigh_ConditionedLeftChoice_model(c) = nansum(firstRaceLoss{c}(indexHighWager), 'all')/ nansum(firstRaceLoss{c}, 'all') - biasWager;
            % Now do P(Choice|Confidence)*P(Confidence)
            % Now find P(Choice|Confidence)*P(Confidence)
            probRight_High_model(c) = nansum(secondRaceLoss{c}(indexHighWager), 'all')/ (nansum(secondRaceLoss{c}(indexHighWager), 'all') + nansum(firstRaceLoss{c}(indexHighWager), 'all'));
            probRight_Low_model(c) = nansum(secondRaceLoss{c}(indexHighWager==0), 'all')/ (nansum(secondRaceLoss{c}(indexHighWager==0), 'all') + nansum(firstRaceLoss{c}(indexHighWager==0), 'all'));
        end
        % most likely no high confidence exist
        probRight_High_model(isnan(probRight_High_model)) = eps;
        probRight_Low_model(isnan(probRight_Low_model)) = eps;

        % Now find P(High|Coherence)
        probHigh_model = nan(1, length(modelCoherences));
         % Probability of High should be the probabilities that higher than theta when either bound is hit 
        for c = 1:length(modelCoherences)
            % Calculate prob high
            probHigh_model(c) = (nansum(firstRaceLoss{c}(indexHighWager), 'all') + nansum(secondRaceLoss{c}(indexHighWager), 'all'))/ (nansum(firstRaceLoss{c}, 'all') + nansum(secondRaceLoss{c}, 'all')) - biasWager;
        end
        % If anything below 0 because of the bias then rezero
        probHigh_model(probHigh_model < 0) = eps;

        % total error (multinomial for choice and wager
        % Calculate distributions
        errorMultinomial = 0;
        for z = 1:length(modelCoherences)
            % Could do it like this P(Conf|Choice)*P(Choice)
% %             probRightHighModel = ProbHigh_ConditionedRightChoice_model(z) * probRight_Coh(z);
% %             probRightLowModel = (1-ProbHigh_ConditionedRightChoice_model(z)) * probRight_Coh(z);
% %             probLeftHighModel = ProbHigh_ConditionedLeftChoice_model(z) * (1-probRight_Coh(z));
% %             probLeftLowModel = (1-ProbHigh_ConditionedLeftChoice_model(z)) * (1-probRight_Coh(z));
% %             totalProb = probRightHighModel + probRightLowModel + probLeftHighModel + probLeftLowModel; 
            % Or do it like P(Choice|Conf)*P(Conf)
            probRightHighModel = probRight_High_model(z) * (probHigh_model(z)); % Bias is already included
            probRightLowModel = probRight_Low_model(z) * (1-(probHigh_model(z)));
            probLeftHighModel = (1-probRight_High_model(z)) * (probHigh_model(z));
            probLeftLowModel = (1-probRight_Low_model(z)) * (1-(probHigh_model(z) ));
            totalProb = probRightHighModel + probRightLowModel + probLeftHighModel + probLeftLowModel;
            
            % normalize incase any descripencies
            probRightHighModel = probRightHighModel / totalProb;
            probRightLowModel = probRightLowModel / totalProb;
            probLeftHighModel = probLeftHighModel / totalProb;
            probLeftLowModel = probLeftLowModel / totalProb;
    
            % Make sure to turn any 0s to eps
            probRightHighModel(probRightHighModel <= 0) = eps;
            probRightLowModel(probRightLowModel <= 0) = eps;
            probLeftHighModel(probLeftHighModel <= 0) = eps;
            probLeftLowModel(probLeftLowModel <= 0) = eps;
            % any 1 to 1-eps
            probRightHighModel(probRightHighModel >= 1) = 1-eps;
            probRightLowModel(probRightLowModel >= 1) = 1-eps;
            probLeftHighModel(probLeftHighModel >= 1) = 1-eps;
            probLeftLowModel(probLeftLowModel >= 1) = 1-eps;
    
            numRightHigh    = sum(data.cohs == modelCoherences(z) & data.choice == 1 & data.wager == 1);
            numRightLow     = sum(data.cohs == modelCoherences(z) & data.choice == 1 & data.wager == 0);
            numLeftHigh     = sum(data.cohs == modelCoherences(z) & data.choice == 0 & data.wager == 1);
            numLeftLow      = sum(data.cohs == modelCoherences(z) & data.choice == 0 & data.wager == 0);
            % Calculate log likelihood (for each coherence)
            errorMultinomial = errorMultinomial + numRightHigh*log(probRightHighModel) + numRightLow*log(probRightLowModel) + numLeftHigh*log(probLeftHighModel) + numLeftLow*log(probLeftLowModel);
        end
      
        % total Erroor
        totalError = -sum(errorRT) + -errorMultinomial;
    elseif whichSolve == 3 % Last way to fit, make everything attached

        % First calculate wager
         % Now find P(High|Coherence, t)
        probHigh_model = nan(1, length(modelCoherences));
         % Probability of High should be the probabilities that higher than theta when either bound is hit 
        for c = 1:length(modelCoherences)
            % Calculate prob high
            probHigh_model(c) = (nansum(firstRaceLoss{c}(indexHighWager), 'all') + nansum(secondRaceLoss{c}(indexHighWager), 'all'))/ (nansum(firstRaceLoss{c}, 'all') + nansum(secondRaceLoss{c}, 'all')) - biasWager;
        end

        % Find the logLL 
        % Avoid making -inf or inf with Log
        probHigh_model(probHigh_model<eps) = eps;
        probHigh_model(probHigh_model>1-eps) = 1-eps;
        % Have everytrial take a prob of making its wager based on coherence
        % Eventually you can just use the symmetry to save run time
        probHigh_trial = nan(1, length(modelCoherences));
        for z = 1:length(modelCoherences)
            probHigh_trial(data.cohs == modelCoherences(z)) = probHigh_model(z); 
        end
        % Binary Wager Outcome 
        % Also fit for conditional reaction times (High vs Low)
        errorWager = sum((data.wager == 1) .* log(probHigh_trial) + (data.wager ==0) .* log(1-probHigh_trial), 'omitnan');

        % Calculate P(Reaction Time|Wager).
        % Not as simple, you have to for each trial calculate only one Prob
        % and added to the others
        
%         rt_Temp = nan(1,length(modelCoherences)); 
%         rt_High_Temp = nan(1,length(modelCoherences)); 
%         rt_Low_Temp = nan(1,length(modelCoherences));
        errorRT_Wager = nan(1,length(modelCoherences));
        errorChoice_HighWagerRT = nan(1,length(modelCoherences));
        errorChoice_LowWagerRT = nan(1,length(modelCoherences));
        for c = 1:length(modelCoherences)
            RTdist = laterFlux{c};
            if modelCoherences(c) > 0 
                totalTimeFlux = conv(RTdist, normpdf([deltaT:deltaT:ti(end)], mu_nonDT, sqrt(var_nonDT))); %include the non-decision time
            elseif modelCoherences(c) <= 0 
                totalTimeFlux = conv(RTdist, normpdf([deltaT:deltaT:ti(end)], mu_nonDT_null, sqrt(var_nonDT))); %include the non-decision time
            end
            
            probHigh_Time_model = (nansum(firstRaceLoss{c} .* indexHighWager, 1) + nansum(secondRaceLoss{c} .* indexHighWager, 1)) ./ (nansum(firstRaceLoss{c}, 1) + nansum(secondRaceLoss{c},1)); %distribution of P(high wager|RT)
            probHigh_Time_model(1:2) = 1; % Correct for the weird anamoly early on.
            % Old way, but works.
            distribution_RThigh_model = (probHigh_Time_model(2:end) .* totalTimeFlux(1:ti(end)*(1/deltaT))) ./ (probHigh_model(c)+eps); % bayes for p(RT|wager==1)
            distribution_RTlow_model = ((1-probHigh_Time_model(2:end)) .* totalTimeFlux(1:ti(end)*(1/deltaT))) ./ (1-(probHigh_model(c) + eps)); % bayes for p(RT|wager==0)
           
            % normalize distribution
            distribution_RThigh_model = distribution_RThigh_model ./ sum(distribution_RThigh_model);
            distribution_RTlow_model = distribution_RTlow_model ./ sum(distribution_RTlow_model);
            
            % Eliminate 0s
            distribution_RThigh_model(distribution_RThigh_model < eps) = eps;
            distribution_RTlow_model(distribution_RTlow_model < eps) = eps; 

            % Now calculate the trial to trial error
            % First find the RTs for High and RTs for low
            highRTs = data.RT(data.wager==1 & data.cohs == modelCoherences(c));
            lowRTs = data.RT(data.wager==0 & data.cohs == modelCoherences(c));
            % Now map onto the distributions you have here
            if deltaT == .01
                highRTs = round(highRTs(~isnan(highRTs))/10);
                lowRTs = round(lowRTs(~isnan(lowRTs))/10);
            elseif deltaT == .001
                highRTs = round(highRTs(~isnan(highRTs)));
                lowRTs = round(lowRTs(~isnan(lowRTs)));
            else
                warning('Your time specificity is not supported')
            end
            errorRT_Wager(c) = sum(log(distribution_RThigh_model(highRTs)), 'omitnan') + sum(log(distribution_RTlow_model(lowRTs)), 'omitnan');

            % Calculate RT and RT|Wager={0,1}
%             rt_Temp(c) = sum(RTdist'/sum(RTdist) .* ti(1:300));
%             rt_High_Temp(c) = sum(distribution_RThigh_model .* ti(1:300));
%             rt_Low_Temp(c) = sum(distribution_RTlow_model .* ti(1:300));
        
            % Now calculate P(Right|Reaction Time, Wager)
            % Hardest one yet
            probRight_HighTime_model = nansum(secondRaceLoss{c} .* indexHighWager, 1) ./ (nansum(secondRaceLoss{c} .* indexHighWager, 1) + nansum(firstRaceLoss{c} .* indexHighWager, 1));
            probRight_LowTime_model = nansum(secondRaceLoss{c} .* indexHighWager==0, 1) ./ (nansum(secondRaceLoss{c} .* indexHighWager==0, 1) + nansum(firstRaceLoss{c} .* indexHighWager==0, 1));
            probRight_HighTime_model=probRight_HighTime_model(2:end); %probRight_HighTime_model(isnan(probRight_HighTime_model)) = 0;
            probRight_LowTime_model=probRight_LowTime_model(2:end); %probRight_LowTime_model(isnan(probRight_LowTime_model)) = 0;
            
            % Now calculate the trial to trial error
            % First find the RTs and Wager for each trial and map it back to
            % the distributions you have here
            probRightConditionedHighRT = probRight_HighTime_model(highRTs); probRightConditionedHighRT(probRightConditionedHighRT<eps) = eps; probRightConditionedHighRT(probRightConditionedHighRT>1-eps) = 1-eps;
            probRightConditionedLowRT = probRight_LowTime_model(lowRTs); probRightConditionedLowRT(probRightConditionedLowRT<eps) = eps; probRightConditionedLowRT(probRightConditionedLowRT>1-eps) = 1-eps;
                        
            % now calculate the error
            errorChoice_HighWagerRT(c) = sum((data.choice(data.wager == 1 & data.cohs==modelCoherences(c))) .* log(probRightConditionedHighRT) + ...
                abs(data.choice(data.wager == 1 & data.cohs==modelCoherences(c)) - 1) .* log(1-probRightConditionedHighRT), 'omitnan');

            errorChoice_LowWagerRT(c) = sum((data.choice(data.wager == 0 & data.cohs==modelCoherences(c))) .* log(probRightConditionedLowRT) + ...
                abs(data.choice(data.wager == 0 & data.cohs==modelCoherences(c)) - 1) .* log(1-probRightConditionedLowRT), 'omitnan');
        end
        totalError = -errorWager + -sum(errorRT_Wager) + -sum(errorChoice_HighWagerRT + errorChoice_LowWagerRT);

    elseif whichSolve == 4 % Solve when Confidence is not binary
        % Because confidence is not binary but instead continous for now
        % make it an independent fit between all variables. 
        % Additionally, fit correct and not prob right. 
        % 1) Correct/Error - binary
        % 2) RT - distribution from model
        % 3) Confidence Rating - distribution from model

        % Variables to save
        errorCorrect = 0;
        errorConf = 0;

        % Calculate error for each coherence
        probCorrect_model = nan(1, length(modelCoherences));
        for c = 1:length(modelCoherences)
            % If the assumption is that Accuracy is based on a subset of
            % the First and Second Race Loss, than that should be
            % considered. Ex: secondRaceLoss{c}(:,1:4000) &
            % firstRaceLoss{c}(:,1:4000), Only consider the first 4 seconds
            % for of Prob(Stimulus|DV, Correct/Incorrect). 
            cutOffDT = 4000/(deltaT*1000);

            % Correct vs Error
            if modelCoherences(c) > 0
                probCorrect_model(c) = nansum(secondRaceLoss{c}(:,1:cutOffDT+1), 'all') ./ (nansum(firstRaceLoss{c}(:,1:cutOffDT+1), 'all') + nansum(secondRaceLoss{c}(:,1:cutOffDT+1), 'all'));
            else
                probCorrect_model(c) = nansum(firstRaceLoss{c}(:,1:cutOffDT+1), 'all') ./ (nansum(firstRaceLoss{c}(:,1:cutOffDT+1), 'all') + nansum(secondRaceLoss{c}(:,1:cutOffDT+1), 'all'));
            end
            % Calculate error
            % Avoid making -inf or inf with Log
            probCorrect_model(probCorrect_model<eps) = eps;
            probCorrect_model(probCorrect_model>1-eps) = 1-eps;
   
            % Binary Wager Outcome 
            % Also fit for conditional reaction times (High vs Low)
            errorCorrect = errorCorrect + sum((data.correct(data.cohs == modelCoherences(c)) == 1) .* log(probCorrect_model(c)) + (data.correct(data.cohs == modelCoherences(c))==0) .* log(1-probCorrect_model(c)), 'omitnan');
            
            % Now reaction time, this is easy just copy and paste the
            % previous method. Even better its been solve already.

            % So really all that is left is confidence :D, given time. This
            % is the best way to do it. 
%             posteriorProb_LatentState_GivenChoice = posteriorProb_LatentState_GivenChoice(~isnan(posteriorProb_LatentState_GivenChoice));
            posteriorProb_LatentState_GivenChoice = round(posteriorProb_LatentState_GivenChoice,2); % Posterior that tells you Conf Rating (round to .01)
            % You should somewhat bin the confidence ratings like by steps
            % of 5. I think this would make things easier. 
            % Built distribution like this 
            confRatings = .5:.01:1; % Play around with this
            distOfConf = nan(1, length(.5:.01:1));
            for cRat = 1:length(confRatings)
                indexForDistOfConf = posteriorProb_LatentState_GivenChoice == confRatings(cRat);
                distOfConf(cRat) = (nansum(firstRaceLoss{c}(indexForDistOfConf), 'all') + nansum(secondRaceLoss{c}(indexForDistOfConf), 'all'))/ (nansum(firstRaceLoss{c}, 'all') + nansum(secondRaceLoss{c}, 'all'));
            end
            % Should be 1
%             distOfConf = distOfConf./sum(distOfConf);
%             keyboard
            % Now that you have the distribution simply make the error
            distOfConf(distOfConf<eps) = eps;
            distOfConf(distOfConf>1-eps) = 1-eps;
            errorConf = errorConf + sum(log(distOfConf(data.wager(data.cohs==modelCoherences(c))-49))); %This only works if Confidence is round to 50-100, no decimals    
        end

        % total error (Binary outcomes)
        totalError = -sum(errorRT) + -sum(errorCorrect) + -errorConf;

    end


else %if theta does not exist 
    % calculate the Probability Right
    choiceParams = [params(1) params(2) urgencyMax];
    probRight_trial = nan(1, length(data.cohs));
    probRight_model = nan(1, length(modelCoherences));
    % Correct and Incorrect Heat Map for LogOdds
    correctMap = cell(1, length(modelCoherences));
    incorrectMap = cell(1, length(modelCoherences));
    tic
    parfor z = 1:length(modelCoherences)
        [temp, rtDist, corr, incorr] = FP_Fitting_Choice(modelCoherences(z), numOfImages, ti, bound, choiceParams);
        probRight_model(z) = temp;
        % calculate RT probability
        rtDist = rtDist(2:end);
        normRTDist = rtDist/sum(rtDist, 'omitnan');
        totalTimeFlux = conv(normRTDist, normpdf(ti(2:end), mu_nonDT, sqrt(var_nonDT)));
        % Normalize totalTimeFlux
        totalTimeFlux = totalTimeFlux/sum(totalTimeFlux);
        totalTimeFlux(totalTimeFlux > 1-eps) = 1-eps; 
        totalTimeFlux(totalTimeFlux<eps) = eps; 
        % Calculate Error 
        tempRT = data.RT(data.cohs == modelCoherences(z));
        tempRT = round(tempRT(~isnan(tempRT)));
        errorRT(z) = sum(log(totalTimeFlux(tempRT)), 'omitnan'); 

        % save 
        correctMap{z} = corr{1};
        incorrectMap{z} = incorr{1};
    end
    toc
    
    probRight_model(probRight_model ==0) = eps;
    % Calculate error
    for z = 1:length(modelCoherences)
        probRight_trial(data.cohs == modelCoherences(z)) = probRight_model(z); 
    end

    % LogLL for Choice
    errorChoice = sum((data.choice == 1) .* log(probRight_trial) + (data.choice ==0) .* log(1-probRight_trial), 'omitnan');

    % Calculate Logodds Correct
    correctAddMap = zeros(size(correctMap{1}));
    incorrectAddMap = zeros(size(correctMap{1}));
    for z = 1:length(modelCoherences)
        correctAddMap = correctAddMap + correctMap{z} .* (1/length(modelCoherences));
        incorrectAddMap = incorrectAddMap + (incorrectMap{z}) .* (1/length(modelCoherences));
    end

    % total error
    totalError = -sum(errorRT) + -errorChoice;
end

%% Calculate the distributions for RT|high and RT|low
% % for m = 1:length(modelCoherences)
% %     probHigh_Time_model = (nansum(firstRaceLoss{m} .* indexHighWager, 1) + nansum(secondRaceLoss{m} .* indexHighWager, 1)) ./ (nansum(firstRaceLoss{m} , 1) + nansum(secondRaceLoss{m},1)); %distribution of P(high wager|RT)
% %     distribution_RThigh_model = (probHigh_Time_model(2:end) .* totalTimeFlux_Cell{m}(2:length(ti))) ./ (probHigh_model(m) + eps); % bayes for p(RT|wager==1)
% %     distribution_RTlow_model = ((1-probHigh_Time_model(2:end)) .* totalTimeFlux_Cell{m}(2:length(ti))) ./ (1-probHigh_model(m) + eps); % bayes for p(RT|wager==0)
% %     % Now calculate the error
% %     highRT = round(data.RT(data.cohs == modelCoherences(m) & data.wager == 1));
% %     lowRT  = round(data.RT(data.cohs == modelCoherences(m) & data.wager == 0)); 
% %     errorRT_High(m) = sum(log(distribution_RThigh_model(highRT)));
% %     errorRT_Low(m) = sum(log(distribution_RTlow_model(lowRT)));
% % end

% total error
% totalError = -sum(errorRT_High) + -sum(errorRT_Low) + -errorWager;
