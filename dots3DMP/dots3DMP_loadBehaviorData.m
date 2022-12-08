function data = dots3DMP_loadBehaviorData(subject,conftask,RTtask)

% clean means no breakfix trials, only 'good' subjects, coherence and headings standardized, normalized conf
% see dots3DMP_cleanHumanData, dots3DMP_cleanMonkeyData
% also keeping raw cohs, headings and conf values as well, since these can vary across sessions/subjects)
% and keeping 1-target trials in monkey data, for comparison in PDW
% analyses

cleaned = 1; % load cleaned version as an option, not implemented yet
if nargin<3 || isempty(RTtask), RTtask = 1; end


if isempty(conftask)
    switch subject
        case 'human', conftask = 1;
        otherwise,    conftask = 2;
    end
end

% now load the data
switch subject

    case 'lucio'
        d = load('lucio_20220301-20221006_clean.mat'); % recent lucio data, PDW + RT


    case 'zarya'
        d = load('zarya_20220301-20220920_clean.mat');

    case 'human'

        if RTtask
            d = load('human_20200213-20220317_RT_clean_Apr2022.mat');    % human RT
        else
            d = load('human_20190625-20191231_nonRT_clean_Apr2022.mat');% human non-RT
        end

    case 'simul'
        error('No simulation dataStruct chosen')
end

data = d.data;