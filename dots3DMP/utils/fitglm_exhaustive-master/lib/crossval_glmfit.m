function [loglik, loglik0] = crossval_glmfit(X, y, glm_args, varargin)
% [loglik, loglik0, cvix] = crossval_glmfit(X, y, glm_args, varargin)
%
% OPTION
% ------
% 'n_sim', 1e2
% 'crossval_method', 'HoldOut'
% 'crossval_args', {0.1}
% 'group', []

% 2016 (c) Yul Kang. hk2699 at columbia dot edu.

S = varargin2S(varargin, {
    'n_sim', 1e2
    'crossval_method', 'HoldOut'
    'crossval_args', {0.1}
    'group', []
    });

n = size(X, 1);
assert(iscolumn(y));
assert(n == length(y));

if isempty(S.group)
    S.group = ones(n,1);
%     [~, ~, S.group] = uique(X); 
end

if verLessThan('matlab', '8.6')
    % Use MATLAB's crossval
    loglik0 = crossval(@glm_train_test, X, y, ...
        lower(S.crossval_method), S.crossval_args{1}, ...
        'stratify', S.group, ...
        'mcreps', S.n_sim); % , ...
%         'options', statset('UseParallel', 'always'));

else
    [incl_train0, incl_test0] = crossvaltf( ...
        S.crossval_method, S.n_sim, S.group, S.crossval_args{:});
    loglik0 = zeros(S.n_sim, 1);

    glm_args = varargin2C(glm_args, {
        'Distribution', 'binomial'
        });
    glm_opt = varargin2S(glm_args);
    assert(strcmp(glm_opt.Distribution, 'binomial'));

    for i_sim = 1:S.n_sim
        incl_test = incl_test0(:, i_sim);
        incl_train = incl_train0(:, i_sim);


        mdl1 = fitglm(X(incl_train,:), y(incl_train), glm_args{:});
        y1 = predict(mdl1, X(incl_test, :));
        y0 = y(incl_test);

        loglik0(i_sim) = glmlik(X(incl_test, :), y0, y1, ...
            glm_opt.Distribution);
    end
end

loglik = mean(loglik0);
end