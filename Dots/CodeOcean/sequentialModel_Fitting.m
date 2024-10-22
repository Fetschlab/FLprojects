function [totalError, probRight, avgRT, probHigh, probRightHigh, probRightLow, avgRTHigh, avgRTLow, probHigh_Corr, probHigh_Inco] ... 
    = sequentialModel_Fitting(data, cohs, tEnd, howToCalculateError, params, extraParamsName)

% Description: Fit choice, RT, and wager data to "serial" model. 
% Computes numerical solution to 1D Fokker-Planck equation,
% to calculate RT distributions and confidence in a 1D bounded accumulator 
% model (aka DDM) using FP4.c by Roozbeh Kiani,
% based on the method of Chang & Cooper, 1970
%
% Originally used in:
% Representation of Confidence Associated with a Decision by Neurons in the
% Parietal Cortex, Science 324, 759 (2009)
%
% Code provided for peer review of:
% M. Vivar-Lazo & C.R. Fetsch, Neural basis of concurrent deliberation
% toward a choice and degree of confidence (2024)

% run using Sequential_Wrapper_NN.m


%%

% (For now only urgency signal for confidence bounds and not choice.) 

 
% Choice
% cohs = [-.512 -.256 -.128 -.064 -.032 0 .032 .064 .128 .256 .512];
K = params(1); % free parameter (transfer to confidence)
uvect = K.*cohs;
t = 0:.01:tEnd;
sigma = 1;
Bup = (params(2)).*ones(size(t));  % Free Parameter (and symmetric for choice)
Blo = -Bup;
y0 = 0;
nonDT_Mean = params(5);
nonDT_Std = .023; % Non-decision time variance (For Genji = .023 for Hanzo = .032)
nonDT_Mean_Asym = nonDT_Mean;

% Any extra parameters
biasWager = 0;
urgencyMax = 0;
if nargin > 5 && ~isempty(extraParamsName)
    extraParams = 0;
    for i = 1:length(extraParamsName)
        if strcmp(extraParamsName{i}, 'urgencySignal') % Should urgency signal be added to both choice or pdw
            extraParams = extraParams + 1;
            urgencyMax = params(5+extraParams); % Highest value the urgency signal can take (theres a good reason to multiply by 1000)
        elseif strcmp(extraParamsName{i}, 'asymmetrical')
            extraParams = extraParams + 1;
            nonDT_Mean_Asym = params(5+extraParams);
        elseif strcmp(extraParamsName{i}, 'wagerOffset')
            extraParams = extraParams + 1;
            biasWager = params(5+extraParams);
        end
    end
end


% Confidence
t = 0:.01:3; %this CAN be a free parameter in a different form of this model
Bup_Conf = (params(3))*ones(size(t)); % Free parameter (asymmetric)
Blo_Conf = (params(4)).*ones(size(t)); % Free Parameter


% Main algorithm from Shadlen lab
% What do we need?
% For 1st Propogation we need only Choice and Reaction Time distribution
[pRightAbs, ~, ~, rightDist, leftDist,pLeftAbs,~,~,~] = runFPonSeries(uvect,t,Bup,Blo,y0,sigma); % Could add an urgency max here 

% For 2nd propogation we need Confidence
% I believe the urgency signal gets added to the bound, as a changing vector.
[pUpAbs_R, ~, ~, upDist_R, loDist_R,pLoAbs_R,~,~,~] = runFPonSeries(uvect,t,Bup_Conf,Blo_Conf,y0,sigma,urgencyMax); % If right bound is hit
[pUpAbs_L, ~, ~, loDist_L, upDist_L,pLoAbs_L,~,~,~] = runFPonSeries(uvect,t,-Blo_Conf,-Bup_Conf,y0,sigma,urgencyMax); % If left bound is hit

% Calculate simple stuff
probRight = pRightAbs./(pRightAbs + pLeftAbs);

% Now this is a little less straight forward because if you hit the bottom
% bound things are flipped. How do you account for that? To flip it you
% might just need to subtract by 1.
% probHigh_R = [pLoAbs_R(cohs<0)./(pUpAbs_R(cohs<0) + pLoAbs_R(cohs<0)); pUpAbs_R(cohs>=0)./(pUpAbs_R(cohs>=0) + pLoAbs_R(cohs>=0))]';
% probHigh_L = [pLoAbs_L(cohs<0)./(pUpAbs_L(cohs<0) + pLoAbs_L(cohs<0)); pUpAbs_L(cohs>=0)./(pUpAbs_L(cohs>=0) + pLoAbs_L(cohs>=0))]';
% or
probHigh_R = (pUpAbs_R./(pUpAbs_R + pLoAbs_R))';
probHigh_L = (pLoAbs_L./(pUpAbs_L + pLoAbs_L))'; %pLoAbs_L is for High Left choices

% Why am I not using "probRight" here? Well that is because P(Bet=High) is
% dependent on the choice bound being absorb so that the subsequent
% accumulator can start. However, P(Right) is a comparison between Right
% and Left which could overestimate the probability of hitting a bound and
% subsequently starting the Wager outcome. So using the Prob of Absorbtion
% is correct. This is correct. However, the other way might be right too.
probHigh = probHigh_R .* pRightAbs' + probHigh_L .* pLeftAbs'; %(1-pRightAbs)'; 
probHigh(probHigh>=1) = 1-eps;
probHigh(probHigh<=0) = eps;

% Do the same for Distributions for RT
upDist = cell(size(cohs)); loDist = cell(size(cohs));
rightUpDist = cell(size(cohs)); rightloDist = cell(size(cohs));
leftUpDist = cell(size(cohs)); leftloDist = cell(size(cohs));
for i = 1:length(cohs) % calculate correct RT distributions for confidence 
    upDist{i} = upDist_R{i} .* pUpAbs_R(i) + upDist_L{i} .* pLoAbs_L(i); % P(Crossing|Conf=High)
    loDist{i} = loDist_R{i} .* pLoAbs_R(i) + loDist_L{i} .* pUpAbs_L(i); % P(Crossing|Conf=Low)
    % Need this for joint calculation
    rightUpDist{i} = (upDist_R{i} .* pUpAbs_R(i)) ./ (upDist_R{i} .* pUpAbs_R(i) + loDist_R{i} .* pLoAbs_R(i)); % P(High|RT)
    rightloDist{i} = (loDist_R{i} .* pLoAbs_R(i)) ./ (loDist_R{i} .* pLoAbs_R(i) + upDist_R{i} .* pUpAbs_R(i));
    % Left
    leftUpDist{i} = (upDist_L{i} .* pLoAbs_L(i)) ./ (upDist_L{i} .* pLoAbs_L(i) + loDist_L{i} .* pUpAbs_L(i));
    leftloDist{i} = (loDist_L{i} .* pUpAbs_L(i)) ./ (upDist_L{i} .* pLoAbs_L(i) + loDist_L{i} .* pUpAbs_L(i));
end

% Compute P(High|Correct) & P(High|Incorrect)
probHigh_Corr = [pLoAbs_L(cohs<0)./(pUpAbs_L(cohs<0) + pLoAbs_L(cohs<0)); pUpAbs_R(cohs>=0)./(pUpAbs_R(cohs>=0) + pLoAbs_R(cohs>=0))]' - biasWager;
probHigh_Inco = [pUpAbs_R(cohs<0)./(pUpAbs_R(cohs<0) + pLoAbs_R(cohs<0)); pLoAbs_L(cohs>=0)./(pUpAbs_L(cohs>=0) + pLoAbs_L(cohs>=0))]' - biasWager;

% Place holder for future use
avgRT = nan(1,length(probHigh));
avgRTHigh = nan(1,length(probHigh));
avgRTLow = nan(1,length(probHigh));

% Calculate Choice|Confidence
probRightHigh = (probHigh_R .* probRight')./probHigh;
probRightLow = ((1-probHigh_R) .* probRight')./(1-probHigh);

% Using Mutinomial fitting?
errorMultinomial = 0;

% In this particular case it might be best to fit on the Conditioned
% Responses but we will try both 
if howToCalculateError == 0 % non conditioned
    % First thing is to find the error in Choices since its a simple
    % binomial
    uniqCoh = cohs;
    errorConf = 0;
    errorChoice = 0;
    errorRT = 0;
    for i = 1:length(uniqCoh)
        if 0 % if doing binomial do this but instead we are running the multinomial fits 
            % choice
            trials = sum(data.cohs == uniqCoh(i));
            n = sum(data.choice==1 & data.cohs == uniqCoh(i));
            tempPRight = probRight(i); tempPRight(tempPRight == 1) = 1-eps; tempPRight(tempPRight==0) = eps;
            errorChoice = errorChoice + n*log(tempPRight) + (trials-n)*log(1-tempPRight);
            % Confidence
            trials = sum(data.cohs == uniqCoh(i));
            n = sum(data.wager==1 & data.cohs == uniqCoh(i));
            tempPHigh = probHigh(i) - biasWager; tempPHigh(tempPHigh == 1) = 1-eps; tempPHigh(tempPHigh==0) = eps;
            errorConf = errorConf + n*log(tempPHigh) + (trials-n)*log(1-tempPHigh);
            errorMultinomial=0;
        else
            % Calculate the joint probabilities
            probRightHighModel  = probRightHigh(i) * (probHigh(i) - biasWager);
            probRightLowModel   = probRightLow(i) * (1-(probHigh(i) - biasWager));
            probLeftHighModel   = (1-probRightHigh(i)) * (probHigh(i) - biasWager);
            probLeftLowModel    = (1-probRightLow(i)) * (1-(probHigh(i) - biasWager));
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

            numRightHigh    = sum(data.cohs == uniqCoh(i) & data.choice == 1 & data.wager == 1);
            numRightLow     = sum(data.cohs == uniqCoh(i) & data.choice == 1 & data.wager == 0);
            numLeftHigh     = sum(data.cohs == uniqCoh(i) & data.choice == 0 & data.wager == 1);
            numLeftLow      = sum(data.cohs == uniqCoh(i) & data.choice == 0 & data.wager == 0);
            % Calculate log likelihood (for each coherence)
            errorMultinomial = errorMultinomial + numRightHigh*log(probRightHighModel) + numRightLow*log(probRightLowModel) + numLeftHigh*log(probLeftHighModel) + numLeftLow*log(probLeftLowModel);
            errorChoice=0;
            errorConf=0;
        end

        % Lastly do the Reaction time (include: Choice RT, Conf RT, and Non-DT
        % First calculate the correct Prob(Crossing Bound at T|Choice)
        probChoiceFlux = rightDist{i}.*probRight(i) + leftDist{i}.*(1-probRight(i)); % there might be a need to normalize the leftDist and righDist
        probConfFlux = upDist{i}.*probHigh(i) + loDist{i}.*(1-probHigh(i));
        totalChoiceConfFlux = conv(probChoiceFlux, probConfFlux); % first do choice and conf time
        % normalize things
        totalChoiceConfFlux = totalChoiceConfFlux/sum(totalChoiceConfFlux);
        totalChoiceConfFlux(totalChoiceConfFlux > 1-eps) = 1-eps; 
        totalChoiceConfFlux(totalChoiceConfFlux<eps) = eps; 
        % Now include non-DT
        if uniqCoh(i) <= 0 
            nonDTFlux = normpdf(t, nonDT_Mean, nonDT_Std);
            totalFlux = conv(totalChoiceConfFlux, nonDTFlux); %include the non-decision time (var_nonDT)
            totalFlux = totalFlux/sum(totalFlux);
            totalFlux(totalFlux > 1-eps) = 1-eps; 
            totalFlux(totalFlux<eps) = eps; 
        else
            nonDTFlux = normpdf(t, nonDT_Mean_Asym, nonDT_Std);
            totalFlux = conv(totalChoiceConfFlux, nonDTFlux); %include the non-decision time (var_nonDT)
            totalFlux = totalFlux/sum(totalFlux);
            totalFlux(totalFlux > 1-eps) = 1-eps; 
            totalFlux(totalFlux<eps) = eps; 
        end

        % Finally calulcate the error
        rtOfInterest = round(data.RT(data.cohs == uniqCoh(i))/10); % ms/time
        errorRT = errorRT + sum(log(totalFlux(rtOfInterest)));
    end

    totalError = -errorRT + -errorChoice + -errorConf + -errorMultinomial; % total Error 
elseif howToCalculateError == 1 % conditioned
    % Calculate Choice|Confidence (I calculated this outside)
%     probRightHigh = (probHigh_R .* probRight')./probHigh;
%     probRightLow = ((1-probHigh_R) .* probRight')./(1-probHigh);
    % Now calculate error that goes with this
    uniqCoh = cohs;
    errorChoiceHigh = 0;
    errorChoiceLow = 0;
    errorRT_High = 0;
    errorRT_Low = 0;
    errorConf = 0;
    for i = 1:length(uniqCoh)
        % high
        ntrials = sum(data.cohs == uniqCoh(i) & data.wager == 1);
        nTrialRight = sum(data.choice==1 & data.wager == 1 & data.cohs == uniqCoh(i));
        tempPRightHigh = probRightHigh(i); tempPRightHigh(tempPRightHigh > 1-eps) = 1-eps; tempPRightHigh(tempPRightHigh<eps) = eps;
        errorChoiceHigh = errorChoiceHigh + nTrialRight*log(tempPRightHigh) + (ntrials-nTrialRight)*log(1-tempPRightHigh);
        % low
        ntrials = sum(data.cohs == uniqCoh(i) & data.wager == 0);
        nTrialRight = sum(data.choice==1 & data.wager == 0 & data.cohs == uniqCoh(i));
        tempPRightLow  = probRightLow(i); tempPRightLow(tempPRightLow > 1-eps) = 1-eps; tempPRightLow(tempPRightLow<eps) = eps;
        errorChoiceLow = errorChoiceLow + nTrialRight*log(tempPRightLow) + (ntrials-nTrialRight)*log(1-tempPRightLow);

        % Calculate RT|Confidence
        probChoiceFlux = rightDist{i}.*probRight(i) + leftDist{i}.*(1-probRight(i)); % there might be a need to normalize the leftDist and righDist
        if uniqCoh(i) <= 0 
            nonDTFlux = normpdf(t, nonDT_Mean, nonDT_Std);
        else
            nonDTFlux = normpdf(t, nonDT_Mean_Asym, nonDT_Std);
        end
%         nonDTFlux = normpdf(t, nonDT_Mean, nonDT_Std);
        
        % Reaction time conditionals (FIX THIS, ITS WRONG)
        totalChoiceConfFlux = conv(probChoiceFlux, upDist{i}); % first do choice and conf time
        % normalize things
        totalChoiceConfFlux = totalChoiceConfFlux/sum(totalChoiceConfFlux);
        totalChoiceConfFlux(totalChoiceConfFlux > 1-eps) = 1-eps; 
        totalChoiceConfFlux(totalChoiceConfFlux<eps) = eps; 
        % Now include non-DT
        totalFlux = conv(totalChoiceConfFlux, nonDTFlux); %include the non-decision time (var_nonDT)
        totalFlux = totalFlux/sum(totalFlux);
        totalFlux(totalFlux > 1-eps) = 1-eps; 
        totalFlux(totalFlux<eps) = eps; 

        % Finally calulcate the error
        rtOfInterest = round(data.RT(data.cohs == uniqCoh(i) & data.wager == 1)/10);
        errorRT_High = errorRT_High + sum(log(totalFlux(rtOfInterest)));

        totalChoiceConfFlux = conv(probChoiceFlux, loDist{i}); % first do choice and conf time
        % normalize things
        totalChoiceConfFlux = totalChoiceConfFlux/sum(totalChoiceConfFlux);
        totalChoiceConfFlux(totalChoiceConfFlux > 1-eps) = 1-eps; 
        totalChoiceConfFlux(totalChoiceConfFlux<eps) = eps; 
        % Now include non-DT
        totalFlux = conv(totalChoiceConfFlux, nonDTFlux); %include the non-decision time (var_nonDT)
        totalFlux = totalFlux/sum(totalFlux);
        totalFlux(totalFlux > 1-eps) = 1-eps; 
        totalFlux(totalFlux<eps) = eps; 

        % Finally calulcate the error
        rtOfInterest = round(data.RT(data.cohs == uniqCoh(i) & data.wager == 0)/10);
        errorRT_Low = errorRT_Low + sum(log(totalFlux(rtOfInterest)));

        % Finally error for wager (which i think might be redundant)
        trials = sum(data.cohs == uniqCoh(i));
        n = sum(data.wager==1 & data.cohs == uniqCoh(i));
        tempPHigh = probHigh(i) - biasWager; tempPHigh(tempPHigh > 1-eps) = 1-eps; tempPHigh(tempPHigh<eps) = eps;
        errorConf = errorConf + n*log(tempPHigh) + (trials-n)*log(1-tempPHigh);
    end

    totalError = -errorChoiceHigh + -errorChoiceLow + -errorRT_High + -errorRT_Low + -errorConf;

elseif howToCalculateError == 2 %P(Choice, RT, Wager) , this is what I used for the results of paper (Vivar-Lazo & Fetsch 2024, BioArchive)
    uniqCoh = cohs;
    errorRT_High = 0;
    errorRT_Low = 0;
    errorConf = 0;
    errorChoice_RTHigh = 0;
    errorChoice_RTLow = 0;
    for i = 1:length(uniqCoh)
        % First wager
        trials = sum(data.cohs == uniqCoh(i));
        n = sum(data.wager==1 & data.cohs == uniqCoh(i));
        tempPHigh = probHigh(i) - biasWager; tempPHigh(tempPHigh >= 1) = 1-eps; tempPHigh(tempPHigh<=0) = eps;
        errorConf = errorConf + n*log(tempPHigh) + (trials-n)*log(1-tempPHigh);

        % Next P(RT|Wager)
        % Calculate RT|Confidence
        probChoiceFlux = rightDist{i}.*probRight(i) + leftDist{i}.*(1-probRight(i)); % there might be a need to normalize the leftDist and righDist
        if uniqCoh(i) <= 0 
            nonDTFlux = normpdf(t, nonDT_Mean, nonDT_Std);
        else
            nonDTFlux = normpdf(t, nonDT_Mean_Asym, nonDT_Std);
        end
        
        % Reaction time conditionals (FIX THIS, ITS WRONG)
        totalChoiceConfFlux = conv(probChoiceFlux, upDist{i}); % first do choice and conf time
        % normalize things
        totalChoiceConfFlux = totalChoiceConfFlux/sum(totalChoiceConfFlux);
        totalChoiceConfFlux(totalChoiceConfFlux > 1-eps) = 1-eps; 
        totalChoiceConfFlux(totalChoiceConfFlux<eps) = eps; 
        % Now include non-DT
        totalFlux = conv(totalChoiceConfFlux, nonDTFlux); %include the non-decision time (var_nonDT)
        totalFlux = totalFlux/sum(totalFlux); % Normalize here
        totalFlux(totalFlux > 1-eps) = 1-eps; 
        totalFlux(totalFlux<eps) = eps; 

        % Finally calulcate the error
        rtOfInterest = round(data.RT(data.cohs == uniqCoh(i) & data.wager == 1)/10);
        errorRT_High = errorRT_High + sum(log(totalFlux(rtOfInterest)));

        totalChoiceConfFlux = conv(probChoiceFlux, loDist{i}); % first do choice and conf time
        % normalize things
        totalChoiceConfFlux = totalChoiceConfFlux/sum(totalChoiceConfFlux);
        totalChoiceConfFlux(totalChoiceConfFlux > 1-eps) = 1-eps; 
        totalChoiceConfFlux(totalChoiceConfFlux<eps) = eps; 
        % Now include non-DT
        totalFlux = conv(totalChoiceConfFlux, nonDTFlux); %include the non-decision time (var_nonDT)
        totalFlux = totalFlux/sum(totalFlux); % Normalize again. 
        totalFlux(totalFlux > 1-eps) = 1-eps; 
        totalFlux(totalFlux<eps) = eps; 

        % Finally calulcate the error
        rtOfInterest = round(data.RT(data.cohs == uniqCoh(i) & data.wager == 0)/10);
        errorRT_Low = errorRT_Low + sum(log(totalFlux(rtOfInterest)));

        % Last P(Right|RT,Wager)
        probRightTime = (rightDist{i}.*probRight(i)) ./ (rightDist{i}.*probRight(i) + leftDist{i}.*(1-probRight(i)));
        probLeftTime = (leftDist{i}.*(1-probRight(i))) ./ (rightDist{i}.*probRight(i) + leftDist{i}.*(1-probRight(i)));
        % Right 
        probRightHighFlux = conv(probRightTime, rightUpDist{i}) ./ probHigh(i); % P(Right|RT, Wager=High)
        probRightLowFlux = conv(probRightTime, rightloDist{i}) ./ probHigh(i); % P(Right|RT, Wager=Low)
        % Left 
        probLeftHighFlux = conv(probLeftTime, leftUpDist{i}) ./ probHigh(i) ; % P(Right|RT, Wager=High)
        probLeftLowFlux = conv(probLeftTime, leftloDist{i}) ./ probHigh(i); % P(Right|RT, Wager=Low)
        % Normalize?
        probRightHighFlux = probRightHighFlux./(probRightHighFlux+probLeftHighFlux);
        probRightLowFlux = probRightLowFlux./(probRightLowFlux+probLeftLowFlux);

        % Convolute the Non-DT
        nonDTFlux = nonDTFlux ./ sum(nonDTFlux);

        probRightHighFlux = conv(probRightHighFlux, nonDTFlux);
        probRightLowFlux = conv(probRightLowFlux, nonDTFlux);
        % Left
        probLeftHighFlux = conv(probLeftHighFlux, nonDTFlux);
        probLeftLowFlux = conv(probLeftLowFlux, nonDTFlux);
        % Normalize right/(right+left)
        probRightHighFlux = probRightHighFlux./(probRightHighFlux+probLeftHighFlux);
        probRightLowFlux = probRightLowFlux./(probRightLowFlux+probLeftLowFlux);

        % Finally calulcate the error (a bit more difficult than you think)
        rtOfInterest_Right = round(data.RT(data.choice == 1 & data.cohs == uniqCoh(i) & data.wager == 1)/10);
        rtOfInterest_Left = round(data.RT(data.choice == 0 & data.cohs == uniqCoh(i) & data.wager == 1)/10);
        probRightHighFlux(probRightHighFlux>=1) = 1-eps;
        probRightHighFlux(probRightHighFlux<=0) = eps;
        errorChoice_RTHigh = errorChoice_RTHigh + sum(log(probRightHighFlux(rtOfInterest_Right))) + sum(log(1-probRightHighFlux(rtOfInterest_Left)));

        rtOfInterest_Right = round(data.RT(data.choice == 1 & data.cohs == uniqCoh(i) & data.wager == 0)/10);
        rtOfInterest_Left = round(data.RT(data.choice == 0 & data.cohs == uniqCoh(i) & data.wager == 0)/10);
        % Control Probs
        probRightLowFlux(probRightLowFlux>=1) = 1-eps;
        probRightLowFlux(probRightLowFlux<=0) = eps;
        errorChoice_RTLow = errorChoice_RTLow + sum(log(probRightLowFlux(rtOfInterest_Right))) + sum(log(1-probRightLowFlux(rtOfInterest_Left)));

        % NOTE: SHOULD YOU CONVOLUTE THE P(CHOICE=RIGHT| RT, WAGER) BY THE
        % NON-DT THEREBY SHIFTING THE DISTRIBUTION??? THIS WOULD ALTERED
        % THE PROBABILITY VALUES FOR EACH RT AND GIVE DIFFERENT
        % LIKELIHOODS,RIGHT?
        % I'm starting to think that you do have to convolute the
        % distribution
    end
    % ErrorChoice_RTLow seems to be wrong, you've to check it that its not
    totalError = -errorConf + -errorRT_High + -errorRT_Low + -errorChoice_RTHigh + -errorChoice_RTLow; 

elseif howToCalculateError == -1 % this is not for fitting but instead calculating averages
    % Lastly if you are not fitting but instead trying to generate average
    % outcomes then you can compute it here
    % You already have the average choice and conf. You just need the
    % average RT.
    for i = 1:length(cohs)
%         if i == 101
%             keyboard
%         end
        probChoiceFlux = rightDist{i}.*probRight(i) + leftDist{i}.*(1-probRight(i)); % there might be a need to normalize the leftDist and righDist

        probConfFlux = upDist{i}.*probHigh(i) + loDist{i}.*(1-probHigh(i));
        totalChoiceConfFlux = conv(probChoiceFlux, probConfFlux); % first do choice and conf time
        % normalize things
        totalChoiceConfFlux = totalChoiceConfFlux/sum(totalChoiceConfFlux);
        totalChoiceConfFlux(totalChoiceConfFlux > 1-eps) = 1-eps; 
        totalChoiceConfFlux(totalChoiceConfFlux<eps) = eps; 
        % Now include non-DT
        if cohs(i) <= 0 
            nonDTFlux = normpdf(t, nonDT_Mean, nonDT_Std);
            totalFlux = conv(totalChoiceConfFlux, nonDTFlux); %include the non-decision time (var_nonDT)
            totalFlux = totalFlux/sum(totalFlux);
            totalFlux(totalFlux > 1-eps) = 1-eps; 
            totalFlux(totalFlux<eps) = eps; 
        else
            nonDTFlux = normpdf(t, nonDT_Mean_Asym, nonDT_Std);
            totalFlux = conv(totalChoiceConfFlux, nonDTFlux); %include the non-decision time (var_nonDT)
            totalFlux = totalFlux/sum(totalFlux);
            totalFlux(totalFlux > 1-eps) = 1-eps; 
            totalFlux(totalFlux<eps) = eps; 
        end
        
        avgRT(i) = sum(totalFlux(1:301)' .* t); % mean RT

        % Reaction time conditionals (FIX THIS, ITS WRONG, is it wrong?)
        totalChoiceConfFlux = conv(probChoiceFlux, upDist{i}); % first do choice and conf time
        % normalize things
        totalChoiceConfFlux = totalChoiceConfFlux/sum(totalChoiceConfFlux);
        totalChoiceConfFlux(totalChoiceConfFlux > 1-eps) = 1-eps; 
        totalChoiceConfFlux(totalChoiceConfFlux<eps) = eps; 
        % Now include non-DT
        if cohs(i) <= 0 
            nonDTFlux = normpdf(t, nonDT_Mean, nonDT_Std);
            totalFlux = conv(totalChoiceConfFlux, nonDTFlux); %include the non-decision time (var_nonDT)
            totalFlux = totalFlux/sum(totalFlux);
            totalFlux(totalFlux > 1-eps) = 1-eps; 
            totalFlux(totalFlux<eps) = eps; 
        else
            nonDTFlux = normpdf(t, nonDT_Mean_Asym, nonDT_Std);
            totalFlux = conv(totalChoiceConfFlux, nonDTFlux); %include the non-decision time (var_nonDT)
            totalFlux = totalFlux/sum(totalFlux);
            totalFlux(totalFlux > 1-eps) = 1-eps; 
            totalFlux(totalFlux<eps) = eps; 
        end

        avgRTHigh(i) = sum(totalFlux(1:301)' .* t); % mean RT

        totalChoiceConfFlux = conv(probChoiceFlux, loDist{i}); % first do choice and conf time
        % normalize things
        totalChoiceConfFlux = totalChoiceConfFlux/sum(totalChoiceConfFlux);
        totalChoiceConfFlux(totalChoiceConfFlux > 1-eps) = 1-eps; 
        totalChoiceConfFlux(totalChoiceConfFlux<eps) = eps; 
        % Now include non-DT
        if cohs(i) <= 0 
            nonDTFlux = normpdf(t, nonDT_Mean, nonDT_Std);
            totalFlux = conv(totalChoiceConfFlux, nonDTFlux); %include the non-decision time (var_nonDT)
            totalFlux = totalFlux/sum(totalFlux);
            totalFlux(totalFlux > 1-eps) = 1-eps; 
            totalFlux(totalFlux<eps) = eps; 
        else
            nonDTFlux = normpdf(t, nonDT_Mean_Asym, nonDT_Std);
            totalFlux = conv(totalChoiceConfFlux, nonDTFlux); %include the non-decision time (var_nonDT)
            totalFlux = totalFlux/sum(totalFlux);
            totalFlux(totalFlux > 1-eps) = 1-eps; 
            totalFlux(totalFlux<eps) = eps; 
        end

        avgRTLow(i) = sum(totalFlux(1:301)' .* t); % mean RT
    end

    % Add bias part
    probHigh = probHigh - biasWager;
    
    % Also calculate Conditionals (Already did this outside)
    % Prob Right | Confidence = {1,0}
%     probRightHigh = (probHigh_R .* probRight')./probHigh;
%     probRightLow = ((1-probHigh_R) .* probRight')./(1-probHigh);
    totalError = 0;

end
