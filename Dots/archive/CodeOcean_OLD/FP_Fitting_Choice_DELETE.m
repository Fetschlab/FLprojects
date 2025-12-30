function [probRight_model, rtDist] = FP_Fitting_Choice(modelCoherences, numOfImages, ti, bound, params)
% Drugowitsch Method for Fitting
% Only Choice
% Copied from "FP_Drugowitsch_Fitting", keeping only what is necessary for
% choice
% 1) You need to make the bound and starting point much lower than the
% original vales because for some reason this algorithm gives funky values
% if you use the original values. And when far enough you don't have use
% Inf you can use 0 to get the amount of DV that has left only one
% boundary.

% Equations needed for Drugowitsch Method 2019
% twoDGauss = @(x, s, mu, sigma, t) 1./(2.*pi.*sqrt(det(t.*sigma))) .* exp(-1./(2).* sum(((x - s - mu.*t)/(sigma.*t)).*(x - s - mu.*t), 2)); 
sj_even = @(j, alpha, k, s0) (1/sin(pi/k))*([sin(j*alpha + pi/k)    sin(j*alpha);           -sin(j*alpha)           -sin(j*alpha - pi/k)])  * s0';
sj_odd = @(j, alpha, k, s0) (1/sin(pi/k))* ([sin(j*alpha)           sin(j*alpha - pi/k);    -sin(j*alpha + pi/k)    -sin(j*alpha)])         * s0';
alpha =@(k) (k-1)/k * pi;
weightj =@(j,mu,sigma,sj,s0) (-1)^j * exp(mu/(sigma) * (sj-s0)');
% Gaussian equation not needed because we just need the CDF

% Introduce Parameters that we are guessing for
sensitivity = params(1);
startingPoint = [params(2), params(2)];

%Extra params for urgency signal (which depend on which kind of urgency)

if length(params) == 4 %This doesnt work yet
    urgencyMax = params(3);
    tauHalf = params(4);
    urgencyVector = abs(urgencyMax) .* ti(1:end)./(ti(1:end)+tauHalf);
elseif length(params) == 3 
    urgencyMax = params(3);
    urgencyVector = ones(1, length(ti)-1) .* abs(urgencyMax)/(length(ti)-1);
    urgencyVector = [0 urgencyVector];
else
    urgencyVector = zeros(1,length(ti));
end

%Set everything else up
k = numOfImages; % Number of Images (Controls the correlation value)
Alpha = alpha(k); 
rho = -cos(pi/k);
sigma=[1 rho;rho 1]; % covariance matrix (You can plug in correlations for CovXY when Variance is 1 for X and Y)
% bound = -.025; %temporary for now (Change this to 0 eventually)

boundForMarginal = inf;

% Run the Algorithm
for c = 1:length(modelCoherences) %parfor
    coh = modelCoherences(c);  %coherence and sensitivity
    mu = [sensitivity*coh, -sensitivity*coh]; %Drift rate
    %Need to add this in here for Parfor loop reasons!
    survivalProb = zeros(1, length(ti)); % Survival Prob for 2D
    flux_RaceOne = zeros(1, length(ti));
    flux_RaceTwo = zeros(1, length(ti));
    
    % Don't runti=0 (that is why t starts @ 2 instead of 1)
    for t = 2:length(ti) %is running parfor loop here faster or at Coh?!?!?!?!?!?!
        indexT=ti(t);
        cdfRest = mvncdf([0 0], startingPoint + (mu + urgencyVector(t)).*indexT, sigma.*indexT); % Probability outside of the Boundaries
        cdfRest_1D = mvncdf([-boundForMarginal bound], [0 0], startingPoint + (mu + urgencyVector(t)).*indexT, sigma.*indexT); 
        cdfRest_2D = mvncdf([bound -boundForMarginal], [0 0], startingPoint + (mu + urgencyVector(t)).*indexT, sigma.*indexT);
        for i = 1:((2*k)-1) %iterate through Reflective Points Needed (Method Of Images)
            if round(i/2) == i/2 %is it even? (Diff Functions for Even or Odd)
                s = sj_even(i, Alpha, k, startingPoint);
                cdfAdd = mvncdf([0 0], s' + (mu + urgencyVector(t)).*indexT, sigma.*indexT);
                cdfAdd_1D = mvncdf([-boundForMarginal bound], [0 0], s' + (mu + urgencyVector(t)).*indexT, sigma.*indexT); 
                cdfAdd_2D = mvncdf([bound -boundForMarginal], [0 0], s' + (mu + urgencyVector(t)).*indexT, sigma.*indexT);
                weight = weightj(i, (mu + urgencyVector(t)), sigma, s', startingPoint);
            else
                s = sj_odd(i, Alpha, k, startingPoint);
                cdfAdd = mvncdf([0 0], s' + (mu + urgencyVector(t)).*indexT, sigma.*indexT);
                cdfAdd_1D = mvncdf([-boundForMarginal bound], [0 0], s' + (mu + urgencyVector(t)).*indexT, sigma.*indexT); 
                cdfAdd_2D = mvncdf([bound -boundForMarginal], [0 0], s' + (mu + urgencyVector(t)).*indexT, sigma.*indexT);
                weight = weightj(i, (mu + urgencyVector(t)), sigma, s', startingPoint);
            end
            % Sometimes weight will be "inf" which will break things, so just make it a high number
            if weight == inf
                weight = realmax;
            elseif weight == -inf
                weight = realmin;
            end
            cdfRest = cdfRest + weight * cdfAdd;
            cdfRest_1D = cdfRest_1D + weight * cdfAdd_1D; % Winner for -Coh
            cdfRest_2D = cdfRest_2D + weight * cdfAdd_2D; % winner for +Coh
        end
        survivalProb(t)     = cdfRest; % P
        flux_RaceOne(t)  = cdfRest_1D; % when CDF for one variable (y) goes to inf its like finding the pdf of the marginal cdf of P(x)
        flux_RaceTwo(t)  = cdfRest_2D; % other cdf variable X
    end
    % Calculate Probability Right 
    probRight_model(c) = sum(flux_RaceTwo) / (sum(flux_RaceTwo) + sum(flux_RaceOne));

    % RT distribution
    death = 1 - survivalProb;
    rtDist = diff(death);
end

% Avoid getting inf or -inf when using log()
probRight_model(probRight_model<eps) = eps;
probRight_model(probRight_model>1-eps) = 1-eps;