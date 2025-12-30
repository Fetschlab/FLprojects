function [mdl, info, mdls] = fitglm_exhaustive(X, y, glm_args, varargin)
% Picks the best model among all 2^n_param possible models.
%
% [mdl, info, mdls] = fitglm_exhaustive(X, y, glm_args, varargin)
% : X is a numeric matrix. y is a numeric vector.
%   glm_args is a cell array of name-value pairs given to fitglm,
%   e.g., {'Distribution', 'binomial', ...
%          'CategoricalVars', {'column1', 'column2'}}
%
% [...] = fitglm_exhaustive(tbl, [], glm_args, varargin)
% : tbl is a table. Use dataset2table to convert datasets.
%   The last column is taken as the response variable.
%
% [...] = fitglm_exhaustive(tbl, ResposeVar, glm_args, varargin)
% : tbl is a table. Use dataset2table to convert datasets.
%   ResponseVar is the name of the column in tbl.
%
% Options
% -------
% ... % 'model_criterion'
% ... % : 'crossval' : cross validates using negative log likelihood
% ... % : 'AIC', 'AICc', 'BIC', 'BICc', 'CAIC' : see mdl.ModelCriterion
% 'model_criterion', 'BIC' 
% 'must_include',    [] % Indices of columns of X to include.
% 'must_exclude',    [] % Indices of columns of X to exclude.
% 'crossval_args',   {}
% 'UseParallel',     'model' % 'model'|'none'
% 'verbose',         true
% 'return_mdls',     (nargout >= 3)
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
%
% 2015 (c) Yul Kang. hk2699 at columbia dot edu.

    if ~exist('glm_args', 'var')
        glm_args = {};
    end
    S = varargin2S(varargin, {
        ... % 'model_criterion'
        ... % : 'crossval' : cross validates using negative log likelihood
        ... % : 'AIC', 'AICc', 'BIC', 'BICc', 'CAIC' : see mdl.ModelCriterion
        'model_criterion', 'BIC' 
        'must_include',    [] % Indices of columns of X to include.
        'must_exclude',    [] % Indices of columns of X to exclude.
        'crossval_args',   {}
        ...
        ... % 'group'
        ... % : Used in cross-validation.
        ... %   Leave empty to treat each unique row as its own group.
        'group',           [] 
        ...
        ... % 'UseParallel'
        ... % : 'auto'|'model'|'none'
        ... %   'auto' by default. Chooses 'model' if n_model > 1000
        'UseParallel',     'auto' 
        'verbose',         true
        'return_mdls',     (nargout >= 3)
        });
    
    if islogical(S.must_include)
        S.must_include = find(S.must_include);
    end
    if islogical(S.must_exclude)
        S.must_exclude = find(S.must_exclude);
    end

    % Construct param_incl_all
    if istable(X)
        assert(ischar(y) || isempty(y));
        var_names = X.Properties.VariableNames;
        if isempty(y)
            y = table2array(X(:,end));
            y_name = var_names{end};
        else
            y_name = y;
            y = X.(y_name);
        end
        X = table2array(X(:, setdiff(var_names, y_name, 'stable')));
        
        glm_args = varargin2C({
            'VarNames', var_names
            'ResponseVar', y_name
            }, glm_args);
    else
        assert(isnumeric(X));
        assert(ismatrix(X));
        
        assert(isnumeric(y) || islogical(y));
        assert(isvector(y));
    end
    
    n_param = size(X, 2);
    n_model = 2 ^ n_param;
    param_incl_all = false(n_model, n_param);
    model_incl = true(n_model, 1);

    for i_model = n_model:-1:1
        param_incl = dec2bin(i_model - 1, n_param) == '1';
        if ~isempty(S.must_include) && any(~param_incl(S.must_include))
            model_incl(i_model) = false;
        end
        if ~isempty(S.must_exclude) && any(param_incl(S.must_exclude))
            model_incl(i_model) = false;
        end
        param_incl_all(i_model, :) = param_incl;
    end

    param_incl_all = param_incl_all(model_incl, :);
    
    if S.verbose
        t_st = tic;
        fprintf('Choosing among %d models began at %s\n', ...
            nnz(model_incl), datestr(now, 30));
    end
    
    if strcmp(S.UseParallel, 'auto')
        if nnz(model_incl) > 1000
            S.UseParallel = 'model';
        else
            S.UseParallel = 'none';
        end
    end

    % Fit
    [ic_all, ic_all0, param_incl_all, mdls] = ...
        fitglm_all(X, y, glm_args, param_incl_all, S);

    if S.verbose
        t_el = toc(t_st);
        fprintf('Chose the best among %d models in %1.1f sec at %s\n', ...
            nnz(model_incl), t_el, datestr(now, 30));
    end
    
    % Output
    ic_all_se = cellfun(@sem, ic_all0);

    [ic_min, ic_min_ix] = min(ic_all);
    
    if S.return_mdls
        mdl = mdls{ic_min_ix};
    else % Estimate it again
        [~,~,mdl] = fitglm_unit( ...
            X, y, glm_args, param_incl_all(ic_min_ix,:), ...
            'none', {}, []);
    end

    param_incl = param_incl_all(ic_min_ix, :);

    info = packStruct(param_incl, ic_min, ic_min_ix, ...
        ic_all, param_incl_all, ...
        ic_all0, ic_all_se);
    info = copyFields(info, S);
end
function [ic_all, ic_all0, param_incl_all, mdls] = ...
        fitglm_all(X, y, glm_args, param_incl_all, S)
    
    n_model = size(param_incl_all, 1);

    ic_all = zeros(n_model, 1);
    ic_all0 = cell(n_model, 1);

    mdls = cell(n_model, 1);
    return_mdls = S.return_mdls;
    
    model_criterion = S.model_criterion;
    crossval_args = S.crossval_args;
    if isempty(S.group)
        [~, ~, group] = unique(X, 'rows');
    else
        group = S.group;
    end

    switch S.UseParallel
        case 'model'
            parfor i_model = 1:n_model
                param_incl = param_incl_all(i_model, :);

                [c_ic, c_ic0, c_mdl] = fitglm_unit( ...
                    X, y, glm_args, param_incl, ...
                    model_criterion, crossval_args, group);

                ic_all(i_model) = c_ic;
                ic_all0{i_model} = c_ic0;

                if return_mdls
                    mdls{i_model} = c_mdl;
                end
            end
        otherwise
            for i_model = 1:n_model
                param_incl = param_incl_all(i_model, :);

                [c_ic, c_ic0, c_mdl] = fitglm_unit( ...
                    X, y, glm_args, param_incl, ...
                    model_criterion, crossval_args, group);

                ic_all(i_model) = c_ic;
                ic_all0{i_model} = c_ic0;

                if return_mdls
                    mdls{i_model} = c_mdl;
                end
            end
    end
end
function [c_ic, c_ic0, c_mdl] = fitglm_unit(X, y, glm_args, param_incl, ...
    model_criterion, crossval_args, group)

    c_mdl = fitglm(X, y, glm_args{:}, ...
        'PredictorVars', find(param_incl));

    switch model_criterion
        case 'crossval'
            if verLessThan('matlab', '8.6')
                glm_args1 = [glm_args(:)', ...
                    {'PredictorVars', find(param_incl)}];
                [c_ic, c_ic0] = bml.stat.crossval_glmfit(X, y, glm_args1, ...
                    'group', group, crossval_args{:});

                % Take negative log likelihood
                c_ic = -c_ic;
                c_ic0 = -c_ic0;
            else
                glm_args1 = [glm_args(:)', ...
                    {'PredictorVars',find(param_incl)}];
                [c_ic, c_ic0] = bml.stat.crossval_glmfit(X, y, glm_args1, ...
                    'group', group, crossval_args{:});

                % Take negative log likelihood
                c_ic = -c_ic;
                c_ic0 = -c_ic0;
            end
            
        case 'none'
            c_ic = nan;
            c_ic0 = nan;

        otherwise
            c_ic = c_mdl.ModelCriterion.(model_criterion);
            c_ic0 = c_ic;
    end
end