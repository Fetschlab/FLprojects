% Miguel Vivar-Lazo
% 03/17/2023
% Fetsch Lab: Explain Fokker Plank Drugowitsch Code

%% I Only use two functions 
% 1) Both assume the bound is 0. I sometimes make the bound -.025 for
% reasons I cannot really explain.
% 2) Urgency Signal is implemented but only works for linear.  I have some
% written for a hyperbolic functions buuuut requires some work.
% 3) I have only use the C code with mesh grid sizes of .025 and time as
% .001 . Ideally any change to that should work so if it doesn't please let
% me know. 

%% FP_Fitting_Choice: Uses Gaussian CDFs to (faster) to find the first
% passage-time distributions to calculate Choice and Reaction time. This is
% matlab code so much easier to dissect if needed.

% Input
% params = [driftRate startingPoint UrgencyMax]
%        = [0:20      0:-inf        0:10000] Urgency Max should be in the
%        1000s, at least for linear
%        making "params" on a vector of two makes the function run with no
%        urgency.
% modelCoherences   = coherences for the model. Although it does loop in the
%                   function I tend to do the looping outside. 
% numOfImages       = normally 4 for the typical -.7 corr
% ti                = vector with length of time (0:.001:3). I always start with 0.
% bound             = -.025. I don't make it 0 but close as possible to 0. We can talk
%                       about this if you like. It's not really a bound.
%                       Named it incorrectly.

% Output
% probRight_model = Probability of Making a right choice (Going to +Drift
%                   Rate bound)
% rtDist = Distribution of Reaction time, you can them find the mean or fit
% using this.

%[probRight_model, rtDist] = FP_Fitting_Choice(modelCoherences, numOfImages, ti, bound, params)

%% FP_Drugowitsch: C code that constructs heatmap
% A little more complicated

% subInd = [x(:) y(:)]; 2D vector with all cartesian mesh points. I
%           normally use meshgrid to make this.
% deltaX = distance between mesh points. I make sure to make X and Y mesh
%           points the same distance.
% length(yi)*length(xi) = is just the number of meshpoints at any Time (t)
% mu = [driftRate*coh, -driftRate*coh]; vector 
% k = number of images (default I use is 4)
% endTi = last time value is seconds 
% deltaT = timesteps in seconds
% s0 = [s0_1D, s0_1D]; Vector with X and Y starting point
% optional: urgencyMax = should be in the 1000s, at least for linear. If no
%                       input the function will work as if no urgency
%                       signal.

% [xytFPMat_FromC] = FP_Drugowitsch(subInd, deltaX, length(yi)*length(xi), mu, k, endTi, deltaT, s0, urgencyMax); %you might then have to multiply by deltaT to actually get the p(x,y,t|parameters)

% Output 
% xytFPMat_FromC = 2D heatmap

% USE THIS LINES AFTER TO MAKE THE 2D HEATMAP INTO 3D HEATMAP (x-y-time)
% xytFPMat_FromC(xytFPMat_FromC<eps) = eps;
% xytFPCell = reshape(xytFPMat_FromC, [length(xi), length(yi), size(xytFPMat_FromC,2)]); 
% Now it should be what you expect. 