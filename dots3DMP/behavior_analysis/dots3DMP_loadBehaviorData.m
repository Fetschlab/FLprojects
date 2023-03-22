function data = dots3DMP_loadBehaviorData(subject,datapath,conftask,RTtask)

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
%         filename = 'lucio_20220301-20221006_clean.mat'; % recent lucio data, PDW + RT (SfN 2022)
        filename = 'lucio_20220301-20230303_clean.mat'; 


    case 'zarya'
        filename = 'zarya_20220301-20220920_clean.mat';

    case 'human'

        if RTtask
%             filename = 'human_20200213-20220317_RT_clean_Apr2022.mat';    % human RT
            filename = 'human_20200213-20211020_RT_clean.mat';
        else
            filename = 'human_20190625-20191231_nonRT_clean_Apr2022.mat';% human non-RT
        end

    case 'simul'
        error('No simulation dataStruct chosen')
end

try
    d = load(fullfile(datapath,filename));
catch 
    fprintf('could not find file %s\n', filename)
end
data = d.data;

fprintf('loaded file %s, %d trials\n', filename, length(data.heading))