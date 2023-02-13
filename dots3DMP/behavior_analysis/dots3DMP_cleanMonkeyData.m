% clean monkey data structure 

clear all; close all

conftask = 2;
RTtask = 1;

subject = 'zarya';
paradigm = 'dots3DMP';
% dateRange = 20210315:20210805; % RT
dateRange = 20220301:20221006;   % for SfN2022

dateRange = 20220512:20230222;

dateRange = 20221207:20221222;


%%
% folder = '/Users/chris/Documents/MATLAB/PLDAPS_data/';
folder = '/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/dataStructs/';

file = [subject '_' num2str(dateRange(1)) '-' num2str(dateRange(end)) '.mat'];
load([folder file], 'data');

try data = rmfield(data,'amountRewardLowConf'); catch, end
try data = rmfield(data,'amountRewardHighConf'); catch, end


% load('/Users/chris/Downloads/lucio_20210315-20210707_clean.mat')

%% temp: CF needs these for now
for k = 1:length(data.filename)
    data.subj{k,:} = subject;
%     data.oneTargChoice(k,:) = 0;
%     data.oneTargConf(k,:) = 0;
end

%%
% struct data has fields:
% filename
% subj: subject code
% choice: 1=left, 2=right, nan = fixation break or otherwise invalid
% heading: angle in deg relative to straight ahead, positive = rightward
% coherence: dot motion coherence aka visual reliability
% modality: stimulus modality: 1=ves, 2=vis, 3=comb
% delta: conflict angle in deg, positive means visual right, vestib left
% correct: was the choice correct, 1 or 0
% conf: confidence rating via saccadic end-point to color-bar target
%       (or in other datasets, PDW)


%% new 04/2020: selecting 'good' data (esp RT)

% some new useful vars
for k = 1:length(data.filename)
    data.date(k,1) = str2double(data.filename{k}(6:13));
    data.subjDate{k,:} = [data.subj{k} data.filename{k}(6:13)];
end


%% Some manual excludes e.g. of bad sessions, non-RT/PDW

excludes_filename = {};
excludes_date = [];

% remove fixation breaks (indicated by nans), or manually excluded filenames
removethese = isnan(data.choice) | isnan(data.RT) | isinf(data.RT) | isnan(data.PDW);
removethese = removethese | ismember(data.filename,excludes_filename) | ismember(data.date,excludes_date);
fnames = fieldnames(data);
for F = 1:length(fnames)
    eval(['data.' fnames{F} '(removethese) = [];']);
end

% should do the trial number based exclusion here, once brfixes are removed

% quick look at blocks, for when some need to be excluded
blocks = unique(data.filename);
nTrialsByBlock = nan(length(blocks),1);
for u = 1:length(blocks)
    nTrialsByBlock(u) = sum(ismember(data.filename,blocks(u)));
end

% discard blocks with <N good trials
N = 50;
removethese = ismember(data.filename,blocks(nTrialsByBlock<N));
for F = 1:length(fnames)
    eval(['data.' fnames{F} '(removethese) = [];']);
end

% quick look at blocks again
blocks = unique(data.filename);
nTrialsByBlock = nan(length(blocks),1);
for u = 1:length(blocks)
    nTrialsByBlock(u) = sum(ismember(data.filename,blocks(u)));
end


%% cull data

mods = unique(data.modality); 

data.delta(data.delta==-2) = -3;
data.delta(data.delta==2) = 3;
deltas = unique(data.delta); % aka conflict angle

% simplify cohs (collapse similar ones)
data.coherence(data.coherence<=0.5) = 0.2;
data.coherence(data.coherence>0.5) = 0.7;
cohs = unique(data.coherence);

% the coh assigned to vestib trials (as a placeholder) depends on which
% cohs were shown in a particular block, so we need to standardize it:
data.coherence(data.modality==1) = cohs(1);

data.heading(abs(data.heading)<0.01) = 0;
% hdgs = unique(data.heading);

switch subject
    case 'zarya'
        hdgs = [-12 -6 -3 -1.5 0 1.5 3 6 12];
    case 'lucio'
        hdgs = [-12 -6 -3 -1.5 0 1.5 3 6 12];
end

% remove the rest
removethese = ~ismember(data.heading,hdgs) | ~ismember(data.coherence,cohs);
for F = 1:length(fnames)
    eval(['data.' fnames{F} '(removethese) = [];']);
end


% fix data.correct (see function description for issue!)
data.correct = dots3DMPCorrectTrials(data.choice,data.heading,data.delta);

% remove one target trials
% removethese = data.oneTargChoice | data.oneTargConf;
% for F = 1:length(fnames)
%     eval(['data.' fnames{F} '(removethese) = [];']);
% end

% try data = rmfield(data,'reward'); catch, end
% try data = rmfield(data,'subj'); catch, end
% try data = rmfield(data,'oneTargChoice'); catch, end
% try data = rmfield(data,'oneTargConf'); catch, end
try data = rmfield(data,'TargMissed'); catch, end
try data = rmfield(data,'subjDate'); catch, end
try data = rmfield(data,'insertTrial'); catch, end
try data = rmfield(data,'confRT'); catch, end


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
save(fullfile(folder,[file(1:end-4) '_clean.mat']),'data');


% CF temp:
% parsedData = dots3DMP_parseData(data,mods,cohs,deltas,hdgs,conftask,RTtask);

%%
% dots3DMP_plots_func_forAumena(parsedData,mods,cohs,deltas,hdgs,conftask,RTtask);

%%
% gfit = dots3DMP_fit_cgauss(data,mods,cohs,deltas,conftask,RTtask);

%%
% dots3DMP_plots_cgauss_func(gfit,parsedData,mods,cohs,deltas,hdgs,conftask,RTtask)
% dots3DMP_plots_cgauss_byCoh(gfit,parsedData,mods,cohs,deltas,hdgs,conftask,RTtask)















