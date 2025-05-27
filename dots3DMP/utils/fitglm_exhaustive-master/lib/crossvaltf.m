function [train, test] = crossvaltf(meth, n_sim, group, varargin)
% [train, test] = crossvalind(meth, n_sim, group, varargin)
%
% train(TR, SIM) : true if TR is included in the training set
%                  on the SIM-th simulation.
% test(TR, SIM) : true if TR is included in the test set
%                 on the SIM-th simulation.

% 2016 (c) Yul Kang. hk2699 at columbia dot edu.

if isnumeric(group) && isscalar(group)
    n_tr = group;
else
    n_tr = size(group, 1);
end

train = false(n_tr, n_sim);
test = false(n_tr, n_sim);

if isempty(group)
    group = ones(n_tr, 1);
end

switch meth
    case 'Kfold'
        error('Not implemented yet!');
        
    case 'HoldOut'
        for i_sim = 1:n_sim
            [train(:, i_sim), test(:, i_sim)] = crossvalind('HoldOut', ...
                group, varargin{:});
        end
        
    otherwise
        error('Not implemented yet!');
end