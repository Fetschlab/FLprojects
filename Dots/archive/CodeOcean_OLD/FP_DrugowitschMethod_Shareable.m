%Miguel Vivar-Lazo
%01/21/2022
%Fetsch Lab
%Solving Fokker Planck Equation Analytical Method: Drugowitsch 2019
%method. Its limited but might be what we just need at this time!

%Note: Authors claim and extension for higher dimensions should be possible, although at this moment I don't know how. 
%Note: This process only works with NON-CHANGING bounds.
%Note: One might be able to add an urgency signal based on a modification
%to the t or mu (Increase mu,senK,or time). Since decreasing the bound is
%analogous to increasing drift rate.

clear

%Plot flags
oneDPlot = 1;
twoDPlot = 1;

%Equations Needed (from paper: Family of closed-form solutions for
%two-dimensional correlated diffusion processes)
%Image Rotataion Formalism (MoI=Method of Images)
twoDGauss = @(x, s, mu, sigma, t) 1./(2.*pi.*sqrt(det(t.*sigma))) .* exp(-1./(2).* sum(((x - s - mu.*t)/(sigma.*t)).*(x - s - mu.*t), 2)); 
sj_even = @(j, alpha, k, s0) (1/sin(pi/k))*([sin(j*alpha + pi/k)    sin(j*alpha);           -sin(j*alpha)           -sin(j*alpha - pi/k)])  * s0';
sj_odd = @(j, alpha, k, s0) (1/sin(pi/k))* ([sin(j*alpha)           sin(j*alpha - pi/k);    -sin(j*alpha + pi/k)    -sin(j*alpha)])         * s0';
alpha =@(k) (k-1)/k * pi;
weightj =@(j,mu,sigma,sj,s0) (-1)^j * exp(mu/(sigma) * (sj-s0)');

% Conditions (Working in 3rd quadrant so the entire grid is negative -x,-y)
k = 4; % parameter for correlation and number of images needed
rho = -cos(pi/k); %correlation
cohs = [-.512 -.256 -.128 -.064 -.032 -eps eps .032 .064 .128 .256 .512]; 
senK = .65; %coherence and sensitivity
mu = [senK*cohs, -senK*cohs]; %Drift rate
sigma=[1 rho;rho 1]; % covariance matrix (You can plug in correlations for CovXY when Variance is 1 for X and Y)
limU = 0; %is always zero
limL = -3; %lowest part of the grid
s0 = [-2, -2]; %Start of propagation (Notice negative; this can change depending on your need for distance to upper limit/bound of 0)
gridLimits = limU + limL; %Grid limits (Think of 3rd Quadrant on a Grid);
deltaX = .01;
deltaY = .01;
deltaT = .01; % 0.001
startT = .001;
endT = 2; %3 %Time is in seconds

%Formulation
Alpha = alpha(k); % Solve for alpha variable (In paper)
yi = 0:-deltaY:gridLimits; % X-Mesh Grid
xi = gridLimits:deltaX:0; % Y-Mesh Grid
ti = startT:deltaT:endT; % Time Grid
xytFPMat = zeros(length(yi), length(xi), length(ti)); % Construct X-Y-T Grid
margCorr= zeros(length(yi), length(ti));
margError= zeros(length(yi), length(ti));

%%
indexY = 0; indexX=0;% Values for Indexing and Placing Prob. Values
% The free parameters should be (senK, rho, Upper Bound, Lower Bound,
% NonDecision Mean and Var)
tic
[x, y] = meshgrid(xi, yi);
subInd = [x(:) y(:)];
%This is running funky
probCoh = 1/length(cohs);
logOddsPlot = 1;
for c = cohs
% for c = cohs(1) %TEMP
    coh = c;  %coherence and sensitivity
    mu = [senK*coh, -senK*coh]; %Drift rate
    tic
    for t = 1:length(ti)
        indexT=ti(t);
        %twoDGauss gives the PDF, but you must turn it to Probability 
        addRest = twoDGauss(subInd,  s0, mu, sigma,  indexT); %Compute 3rd Quadrant (or First) Value or Image
        for i = 1:((2*k)-1) %iterate through Reflective Points Needed (Method Of Images)
            if round(i/2) == i/2 %is it even? (Diff Functions for Even or Odd)
                s = sj_even(i, Alpha, k, s0);
                restAdd = twoDGauss(subInd,  s', mu, sigma, indexT);
                weight = weightj(i, mu, sigma, s', s0);
            else
                s = sj_odd(i, Alpha, k, s0);
                restAdd = twoDGauss(subInd,  s', mu, sigma, indexT);
                weight = weightj(i, mu, sigma, s', s0);
            end
            addRest = addRest + weight.*restAdd; % Add the weighted image to the original image
        end
        addRest = addRest.*(deltaX*deltaY); %addRest/sum(addRest);%normalize the PDF (IT ACTUALLY MIGHT NOT BE GOOD TO NORMALIZE IN THIS MANNER, probabilities should be leaving as time goes on);
        xytFPMat(:, :, t) = reshape(addRest, [length(yi), length(xi)]); % Add Probability to the Grid (Fills colmns first)
    end
    toc
    xytFPMat(xytFPMat<eps) = eps; %Remove any small negative values, and zeros (to apply Log function)
    %Form LogOdds (P[Correct Motion]/P[Incorrect Motion]))
    if c > 0 %Marginalization for Log Odds: LET ME KNOW IF THIS IS INCORRECT
        margCorr = margCorr + flipud(squeeze(sum(xytFPMat, 1))) .* probCoh; %correct
        margError = margError + squeeze(sum(xytFPMat, 2)) .* probCoh; % incorrect
    else
        margCorr = margCorr + squeeze(sum(xytFPMat, 2)) .* probCoh; %correct
        margError = margError + flipud(squeeze(sum(xytFPMat, 1))) .* probCoh; %incorrect
    end
end
toc % this takes about 7 mins for dt=0.001, 30s for 0.01

%% Convert to 1D for x and y, and plot a 2D
% function plot2D()
if oneDPlot == 1
    firstAccum  = flipud(squeeze(sum(xytFPMat, 1)));
    secondAccum = squeeze(sum(xytFPMat, 2));
    figure; set(gcf,'position',[558         539        1196         409])
    subplot(1,2,1);
    imagesc(log(firstAccum)); colormap jet; colorbar;
    yticks([1:100:length(yi)]); yticklabels({'20', '10','0','-10'});
    xlabel('Time (ms)')
    ylabel('Accumulated Evidence')
    title(['Correct Accumulator: Coh = ' num2str(c)])
    subplot(1,2,2);
    imagesc(log(secondAccum)); colormap jet; colorbar;
    yticks([1:100:length(yi)]); yticklabels({'20', '10','0','-10'});
    xlabel('Time (ms)')
    ylabel('Accumulated Evidence')
    title(['Incorrect Accumulator: Coh = ' num2str(c)])
end
% %% Plot Map Through Space
if twoDPlot == 1
    for z = length(ti):-round(length(ti)/10):1
        figure;
        imagesc(log(xytFPMat(:,:,z))); colormap jet; colorbar();
        xlabel('Accumulator 1')
        ylabel('Accumulator 2')
        yticks([1:99:300]); yticklabels({'20', '10','0','-10'});
        xticks([1:99:300]); xticklabels(fliplr({'20', '10','0','-10'}));
        title(num2str(c));
    end
end


%% Log Odd Plot, which doesn't look like Kaini's :,(
if logOddsPlot == 1
    logOdds = log(margError./margCorr); %./margError;
    figure; imagesc(logOdds); colormap jet; colorbar;
end


%% Now try simulating data on the model (How will you add urgency signal?)
nTrials = 1e5;
simCoh = [eps .032 .064 .128 .256 .512]; simDir = [0 180];
simCohN = randsample(simCoh, nTrials, 'true');
simDirN = randsample(simDir, nTrials, 'true');
simSignCoh = simCohN .* cosd(simDirN);

simChoice = nan(1, nTrials); simRT = nan(1,nTrials); simPDW = nan(1,nTrials);
% Generate two random accumulators with some correlation
K=.07; %senK;
weightTheta = .85;
theta = weightTheta*max(logOdds(:)) + (1-weightTheta)*min(logOdds(:)); %Is distance from boundary high here?
limL = -3;
for i=1:nTrials
    momEvid = mvnrnd([simSignCoh(i) , -simSignCoh(i)], sigma, 3000);
    %Whose Bound was hit first?
    dv = s0;%based off model
    for t = 1:length(momEvid)
        dv = dv + momEvid(t,:).*K; %Turn to DV
        %Make sure to not go over reflective bound
        if dv(1) < limL || dv(2) < limL
            dv(dv < limL) = limL;
        end
        %Has a bound been hit?
        if any(dv >= 0)
            accumWon = find(dv >= 0);
            if length(accumWon) > 1
                [~, accumWon] = max(dv);
            end
            break
        end
    end
    simChoice(i) = abs(accumWon-2);
    %Reaction Time
    simRT(i) = t;
    %High or Low bet? Here we have to use Heat Map made above
    % First which is the losing accumulator?
    loserAcc = abs(accumWon-3);
    dvLoser = round(dv(:, loserAcc) / deltaX); %Loser DV, convert to appropriate index value using grid partition value (deltaX or deltaY)
    %Is this DV higher than the probability located here
    logOddPDW = logOdds(abs(dvLoser)+1, simRT(i));
    simPDW(i) = logOddPDW > theta; %Is it greater than theta or critereon?
end
%% %Map results from Simulation (Simple Choice and RT conditioned on Bet)
plotSimSignCoh = simSignCoh; plotSimSignCoh(plotSimSignCoh == eps) = 0; plotSimSignCoh(plotSimSignCoh==-eps)=0;
[avgChoiceSim, ~, stdChoiceSim]   = behavioralAverages(simChoice, plotSimSignCoh); % choice
[avgChoiceHighSim, ~, stdChoiceHighSim]   = behavioralAverages(simChoice(simPDW==1), plotSimSignCoh(simPDW==1)); % choice Low
[avgChoiceLowSim, ~, stdChoiceLowSim]   = behavioralAverages(simChoice(simPDW==0), plotSimSignCoh(simPDW==0)); % choice High

[avgRTSim, ~, stdRTSim]   = behavioralAverages(simRT, plotSimSignCoh); % RT
[avgRTHighSim, ~, stdRTHighSim]   = behavioralAverages(simRT(simPDW==1), plotSimSignCoh(simPDW==1)); % RT High
[avgRTLowSim, ~, stdRTLowSim]   = behavioralAverages(simRT(simPDW==0), plotSimSignCoh(simPDW==0)); % RT Low

[avgPDWSim, ~, stdPDWSim]   = behavioralAverages(simPDW, plotSimSignCoh); % PDW

figure(100); 
errorbar(unique(plotSimSignCoh), avgChoiceSim, stdChoiceSim, 'k', 'Linewidth', 2); hold on;
errorbar(unique(plotSimSignCoh), avgChoiceHighSim, stdChoiceHighSim, 'r', 'Linewidth', 2);
errorbar(unique( plotSimSignCoh(simPDW==0)), avgChoiceLowSim, stdChoiceLowSim, 'b', 'Linewidth', 2); hold off
xlabel('Coherence'); ylabel('Proportion Right'); legend('All', 'High Bet', 'Low Bet');


figure(101);
errorbar(unique(plotSimSignCoh), avgRTSim, stdRTSim, 'k', 'Linewidth', 2); hold on;
errorbar(unique(plotSimSignCoh), avgRTHighSim, stdRTHighSim, 'r', 'Linewidth', 2);
errorbar(unique(plotSimSignCoh(simPDW==0)), avgRTLowSim, stdRTLowSim, 'b', 'Linewidth', 2); hold off;
xlabel('Coherence'); ylabel('Reaction Time'); legend('All', 'High Bet', 'Low Bet');


figure(103);
errorbar(unique(plotSimSignCoh), avgPDWSim, stdPDWSim, 'k', 'Linewidth', 2);
xlabel('Coherence'); ylabel('Wager');
%% Look At PDW Avg based On Error and Correct
correct = simChoice == cosd(simDirN);
condition = correct == 1;
[avgPDWCorrSim, ~, stdPDWCorrSim]   = behavioralAverages(simPDW(condition), abs(plotSimSignCoh(condition))); % PDW for Correct
condition = correct == 0;
[avgPDWErrSim, ~, stdPDWErrSim]   = behavioralAverages(simPDW(condition), abs(plotSimSignCoh(condition))); % PDW for Correct

figure(104);
errorbar(unique(abs(plotSimSignCoh)), avgPDWCorrSim, stdPDWCorrSim, 'Linewidth', 2); hold on;
errorbar(unique(abs(plotSimSignCoh)), avgPDWErrSim, stdPDWErrSim, 'Linewidth', 2); hold off;
xlabel('Coherence'); ylabel('Wager'); legend('Correct', 'Error');
%%
function [avgDepVar, N, stderr, uniqIndVar] = behavioralAverages(depVariable, indVariable)
uniqIndVar  = unique(indVariable);
avgDepVar   = nan(1,length(uniqIndVar));
N           = nan(1,length(uniqIndVar));
stderr      = nan(1,length(uniqIndVar));
for i = 1:length(uniqIndVar)
    avgDepVar(i)    = mean(depVariable(uniqIndVar(i) == indVariable));
    N(i)            = sum(uniqIndVar(i) == indVariable);
    stderr(i)       = std(depVariable(uniqIndVar(i) == indVariable))/sqrt(N(i));
end
end