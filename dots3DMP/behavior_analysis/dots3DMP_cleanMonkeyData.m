% clean monkey data structure 

clear all; close all

conftask = 2;
RTtask = 1;

subject = 'lucio';
% paradigm = 'dots3DMP';
% dateRange = 20210315:20210805; % RT
% dateRange = 20220512:20230605;
% dateRange = 20221201:20230601;

%%
% folder = '/Users/chris/Documents/MATLAB/PLDAPS_data/';
[file, path] = uigetfile;
load(fullfile(path, file));

%% new 04/2020: selecting 'good' data (esp RT)

% some new useful vars
% for k = 1:length(data.filename)
%     data.date(k,1) = str2double(data.filename{k}(6:13));
%     data.subjDate{k,:} = [data.subj{k} data.filename{k}(6:13)];
% end

% SJ 10-2023
data.date = cellfun(@(x) str2double(x(6:13)), data.filename);
data.subjDate = cellfun(@(x, y) [x y(6:13)], data.subj, data.filename, 'UniformOutput', false);


%% Some manual excludes e.g. of bad sessions, non-RT/PDW

% TODO convert these to functions

excludes_filename = {};
excludes_date = [];

% remove fixation breaks (indicated by nans), or manually excluded filenames
removethese = isnan(data.choice) | isnan(data.RT) | isinf(data.RT) | isnan(data.PDW);
removethese = removethese | ismember(data.filename,excludes_filename) | ismember(data.date,excludes_date);
fnames = fieldnames(data);
for F = 1:length(fnames)
    data.(fnames{F})(:, removethese) = [];
end

% do the trial number based exclusion here, once brfixes are removed

blocks = unique(data.filename);
nTrialsByBlock = nan(length(blocks),1);
for u = 1:length(blocks)
    nTrialsByBlock(u) = sum(ismember(data.filename,blocks(u)));
end

% discard blocks with <N good trials
N = 50;
removethese = ismember(data.filename,blocks(nTrialsByBlock<N));
for F = 1:length(fnames)
    data.(fnames{F})(:, removethese) = [];
end

% quick look at blocks again
blocks = unique(data.filename);
nTrialsByBlock = nan(length(blocks),1);
for u = 1:length(blocks)
    nTrialsByBlock(u) = sum(ismember(data.filename,blocks(u)));
end


%% cull data

mods = unique(data.modality); 

data.delta(data.delta<0) = -3;
data.delta(data.delta>0) = 3;
deltas = unique(data.delta); % aka conflict angle

% simplify cohs (collapse similar ones)
data.coherence(data.coherence<=0.5) = 0.2;
data.coherence(data.coherence>0.5) = 0.7;
cohs = unique(data.coherence);

% the coh assigned to vestib trials (as a placeholder) depends on which
% cohs were shown in a particular block, so we need to standardize it:
data.coherence(data.modality==1) = cohs(1);

data.heading(abs(data.heading)<0.01) = 0;
hdgs = unique(data.heading);


% 05/2023
% re-assign the slightly different headings from cue conflict blocks
% (because I switched from 9 headings overall to 7...)

fnames = fieldnames(data);
if strcmp(subject, 'zarya')
    removethese = abs(data.heading) > 4.5 & abs(data.heading) < 8;
    for F = 1:length(fnames)
        data.(fnames{F})(:, removethese) = [];
    end
end

%%
% keep the raw heading values
data.heading_original = data.heading;

switch subject
    case 'lucio'
        hdg_ranges = [0 1 2 3 4 5 8 10 12];
        newhdgvals = [0 1.5 1.5 3 3 6 6 12 12]; % match length of hdg_ranges

    case 'zarya'
        hdg_ranges = [0 1.5 3.873 10 30];
        newhdgvals = [0 1.5 3.873 10 10]; % match length of hdg_ranges
end

hdg_index = interp1(hdg_ranges, 1:numel(hdg_ranges), abs(data.heading), 'nearest');
data.heading = newhdgvals(hdg_index) .* sign(data.heading);


% remove the rest
% removethese = ~ismember(data.heading,hdgs) | ~ismember(data.coherence,cohs);
% for F = 1:length(fnames)
%     data.(fnames{F})(:, removethese) = [];
% end


% fix data.correct (see function description for issue!)
% irrelevant after 04/2022
% data.correct = dots3DMPCorrectTrials(data.choice,data.heading,data.delta);

% remove one target trials
% leaving in, these can be removed for analysis later if desired
% removethese = data.oneTargChoice | data.oneTargConf;
% for F = 1:length(fnames)
%     data.(fnames{F})(:, removethese) = [];
% end

try data = rmfield(data,'reward'); catch, end
% try data = rmfield(data,'subj'); catch, end
% try data = rmfield(data,'oneTargChoice'); catch, end
% try data = rmfield(data,'oneTargConf'); catch, end
try data = rmfield(data,'TargMissed'); catch, end
try data = rmfield(data,'subjDate'); catch, end
try data = rmfield(data,'insertTrial'); catch, end
try data = rmfield(data,'confRT'); catch, end
try data = rmfield(data,'amountRewardHighConf'); catch, end

% sorted_fnames = {'filename','subj','date','heading','modality','coherence','delta','choice','RT','PDW','correct'};
% data = orderfields(data,sorted_fnames);


%% check block and day trial counts
blocks = unique(data.filename);
nTrialsByBlock = nan(length(blocks),1);
for u = 1:length(blocks)
    nTrialsByBlock(u) = sum(ismember(data.filename,blocks(u)));
end

dates = unique(data.date);
nTrialsByDate = nan(length(dates),1);
for u = 1:length(dates)
    nTrialsByDate(u) = sum(ismember(data.date,dates(u)));
end



%% save it
save(fullfile(path,[file(1:end-4) '_clean.mat']),'data','nTrialsByDate','nTrialsByBlock', '-v7.3');




