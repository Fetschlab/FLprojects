

% build simulation of dots3DMP behavioral computation (non-DDM)


% analytic solution


% heading is either L or R (relevant dimension)
% L or R are equally likely

% observer makes noisy measurement of heading (in either/both modalities)
% model as gaussian centered on true heading, with variance sigma

% p(C|h) = P(h|C)*P(C) / P(h)
% P(C) & P(h) are uniform, so normalized likelihood P(h|C) is equivalent to
% posterior p(C|h)

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