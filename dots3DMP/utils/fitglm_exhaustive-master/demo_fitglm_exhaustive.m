% Demo for fitglm_exhaustive().
%
% WARNING:
% Returning mdls (all models) can be memory intensive. 
% When size(X) is about 1500 x 17 and 2^17 models are fitted,
% mdls can take up >50GB.
% Use it only when you have enough memory.
%
% NOTE 1:
% Fitting all possible models can be impractical when size(X,2) > 25.
% - Try reducing the dimensionality using PCA.
% - If you have a priori reasons to believe some columns should always be
%   included, use the 'must_include' option.
%
% NOTE 2:
% Estimate the time and memory expenditure first by using a small subset of
% columns, like:
%
%     tic;
%     [mdl, info, mdls] = fitglm_exhaustive(X(:,1:8), ...)
%     toc;
%
%     whos mdls
%
% Then estimate the time and memory needed by multiplying 
% the elapsed time and mdls's size in the memory (Bytes) by 2^(size(X,2)-8).

% 2015-2016 (c) Yul Kang. hk2699 at columbia dot edu.

%%
addpath(genpath('lib'));

%% Preparation
n = 1e4;
X = rand(n,3) - 0.5;
y = X(:,1) * 0.4 + X(:,3) * 0.6 + 0.1;
y = max(min(y, 10), -10);
p = exp(y) ./ (1 + exp(y));
ch = rand(size(p)) < p;

%% BIC
[mdl, info, mdls] = fitglm_exhaustive(X, ch, {'Distribution', 'binomial'}, ...
    'model_criterion', 'BIC')

%% BIC without mdls output - use when you have many columns in X.
tic;
[mdl, info] = fitglm_exhaustive(X, ch, {'Distribution', 'binomial'}, ...
    'model_criterion', 'BIC')
toc;

%% BIC, Using Table
% Use dataset2table() to convert datasets.
column1 = X(:,1);
column2 = X(:,2);
column_unused = X(:,3);
tbl = table(column1, column2, column_unused, y);

[mdl, info] = fitglm_exhaustive(tbl, 'y', {'Distribution', 'binomial'}, ...
    'model_criterion', 'BIC', ...
    'UseParallel', 'none')

%% Under construction

% %% cross-validation with default settings
% [mdl, info, mdls] = fitglm_exhaustive(X, ch, {'Distribution', 'binomial'}, ...
%     'model_criterion', 'crossval')
% 
% %% cross-validation with 20 simulations of 50% holdout
% [mdl, info, mdls] = fitglm_exhaustive(X, ch, {'Distribution', 'binomial'}, ...
%     'model_criterion', 'crossval', ...
%     'crossval_method', 'HoldOut', ...
%     'crossval_args', {0.5}, ...
%     'n_sim', 20)
% 
% %% 5-fold cross-validation
% [mdl, info, mdls] = fitglm_exhaustive(X, ch, {'Distribution', 'binomial'}, ...
%     'model_criterion', 'crossval', ...
%     'crossval_method', 'Kfold', ...
%     'crossval_args', {5})
