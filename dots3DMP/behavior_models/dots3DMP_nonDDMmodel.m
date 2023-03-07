

% build simulation of dots3DMP behavioral computation (non-DDM)


% analytic solution


% heading is either L or R (relevant dimension)
% L or R are equally likely

% observer makes noisy measurement of heading (in either/both modalities)
% model as gaussian centered on true heading, with sensory noise
% or model based on readout of poisson population...

% p(C|h) = P(h|C)*P(C) / P(h)
% P(h) is exponential around 0? lambda*exp(-lambda*x)


% stimulus conditions
mods  = [1 2 3];        % stimulus modalities: ves, vis, comb
cohs  = [0.4 0.8];      % visual coherence levels (these are really just labels, since k's are set manually)
hdgs  = [-12 -6 -3 -1.5 0 1.5 3 6 12];
% deltas = [-3 0 3];    % conflict angle; positive means vis to the right
deltas  = 0;

nreps = 20;

[hdg, modality, coh, delta, ntrials] = dots3DMP_create_trial_list(hdgs,mods,cohs,deltas,nreps,0); % don't shuffle

nstep = 200;


% need to calculate likelihoods for vis and ves separately (in combined
% case), then combine them


% flat prior, or decaying exponential around 0...

hdg_vec = linspace(hdgs(1),hdgs(end),nstep);
hdgMat  = kron(hdg_vec',ones(simntrl*ntrials),1);
hdgMat  = reshape(hdgMat,simntrl,ntrials);

simntrl = 10000; 

% need to create hdgMat and est_h

sigma_p = sqrt(sigma_stim^2 + sigma_meas^2);
sim_lh  = 1./(2*pi*sigma_p^2) .* exp(-(est_h - hdgMat).^2) ./ (2*sigma_p^2);

sim_post = bsxfun(@rdivide,sim_lh,sum(sim_lh,3));





% params to fit

% sigma_scale


% mu_vis = cohs * k
% mu_ves = mean(kvis)


% mu_comb and sigma_comb are Bayes calculated from vis + ves


% MAP or area comparison to zero tells choice, also compare to simulated
% zero trials

% could add beta distribution parameter for decision noise

% how to calculate conf?
% MAP or area then comparison to criterion
% what about heuristic method?

% what about RT? separately? some kind of cut-off for integration time based on acc/vel, +
% decision-noise
% maybe this would fail to fit well because RT would be independent of sampling (i.e.
% wouldn't see any average difference between correct and errors in distribution of
% posterior

% noisy estimate of peak of acc/vel for ves/vis, weighted comb for comb
% also wouldn't have any relationship between accuracy and RT?
%
% skewed distribution estimate of peak with sigma also dependent on sigma
% 





% exponential prior, a la Fetsch 2009
x = linspace(0,12,100);
mu = 0.16;

prior = mu*exp(-mu.*x);
prior = [fliplr(prior) prior];
x = [-fliplr(x) x];

prior = prior ./ sum(prior);

figure; plot(x, prior)