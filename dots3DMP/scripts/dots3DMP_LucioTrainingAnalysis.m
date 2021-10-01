% offline analysis wrapper for PLDAPS data, dots3DMP paradigm
% lucio training data, SJ 07-16-2020

clear; close all

RTtask = 1;
conftask = 2; % 1=colorbars, 2=PDW

subject = 'lucio';
paradigm = 'dots3DMP';
dateRange = 20210315:20210512;


%%
% change to whatever local folder data is stored in
% folder = '/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/';
folder = '/Users/chris/Documents/MATLAB/PLDAPS_data/';

file = [subject '_' num2str(dateRange(1)) '-' num2str(dateRange(end)) '.mat'];
load([folder file], 'data');

% some new useful vars
for k = 1:length(data.filename)
    data.date(k,1) = str2double(data.filename{k}(6:13));
    data.subjDate{k,:} = [data.subj{k} data.filename{k}(6:13)];
end

%% manual excludes

% block noted as bad, e.g. always betting high, technical errors, or bias
% excludes_filename = {'lucio20200727dots3DMP0946','lucio20200727dots3DMP1024',...
%     'lucio20200819dots3DMP1200','lucio20200819dots3DMP1220'};
% excludes_subjDate = {'lucio20200625','lucio20200724','lucio20200729', 'lucio20200820',...
%     'lucio20200827','lucio20200902','lucio20200910'};

% excludes_filename = {'lucio20210315dots3DMP1112_basic', 'lucio20210319dots3DMP1332_basic'};
% excludes_subjDate = {};
% 
% cutoff = 20210401; % quickly exclude everything before this date
% removethese = ismember(data.filename,excludes_filename) | ...
%     ismember(data.subjDate,excludes_subjDate) | data.date<cutoff; %#ok<NASGU>
% fnames = fieldnames(data);
% for F = 1:length(fnames)
%     eval(['data.' fnames{F} '(removethese) = [];']);
% end

%% remove invalid trials (fixation breaks (which gives nans))
removethese = isnan(data.choice) | isnan(data.PDW) | data.oneConfTargTrial;
% removethese = isnan(data.choice); % keep nan PDWs for choice curves?
fnames = fieldnames(data);
for F = 1:length(fnames)
    eval(['data.' fnames{F} '(removethese) = [];']);
end

blocks = unique(data.filename);
nTrialsByBlock = nan(length(blocks),1);
for u = 1:length(blocks)
    nTrialsByBlock(u) = sum(ismember(data.filename,blocks(u)));
end

% we can be pretty sure blocks with <N trials (say, 50) are to be discarded
removethese = ismember(data.filename,blocks(nTrialsByBlock<80));
for F = 1:length(fnames)
    eval(['data.' fnames{F} '(removethese) = [];']);
end

% recompute blocks and counts
blocks = unique(data.filename);
nTrialsByBlock = nan(length(blocks),1);
for u = 1:length(blocks)
    nTrialsByBlock(u) = sum(ismember(data.filename,blocks(u)));
end
% 
% figure;
% hist(nTrialsByBlock,15);
% xlabel('Session');
% ylabel('NumTrials');

%% manual exclude specific blocks
% 
i  = [1 4:10 12:13];
removethese = ismember(data.filename,blocks(i));
for F = 1:length(fnames)
    eval(['data.' fnames{F} '(removethese) = [];']);
end
%% manual selection of specific session/s

% i  = 23;
% removethese = ~ismember(data.filename,blocks(i));
% for F = 1:length(fnames)
%     eval(['data.' fnames{F} '(removethese) = [];']);
% end

% keep_SubjDates = {'lucio20200727','lucio20200803','lucio20200805','lucio20200807','lucio20200819',...
%     'lucio20200820','lucio20200821','lucio20200901','lucio20200922','lucio20200923'};
% removethese = ~ismember(data.subjDate,keep_SubjDates);

% keep_Subjfiles = {'lucio20200715dots3DMP1005'};
% removethese = ~ismember(data.filename,keep_Subjfiles);

% for F = 1:length(fnames)
%     eval(['data.' fnames{F} '(removethese) = [];']);
% end

%% cull data

mods = unique(data.modality); 

% only a few 0.4 trials, so reassign to 0.3
data.coherence(data.coherence<=0.4) = 0.4; 
data.coherence(data.coherence>=0.7) = 0.8; 

cohs = unique(data.coherence);

% remove the rest
removethese = ~ismember(data.coherence,cohs) & data.modality~=1;
for F = 1:length(fnames)
    eval(['data.' fnames{F} '(removethese) = [];']);
end
    
% the coh assigned to vestib trials (as a placeholder) depends on which
% cohs were shown in a particular block, so we need to standardize it:
data.coherence(data.modality==1) = cohs(1);

deltas = unique(data.delta); % aka conflict angle

% group some of the hdgs together to simplify things...
% data.heading(data.heading>-5 & data.heading<=-4) = -4;
% data.heading(data.heading<5 & data.heading>=4) = 4;

% data.heading(abs(data.heading)<=eps*1.5) = 0;

% hdgs = unique(data.heading);
% hdgs = [-12 -6 -4 -3 -2 -1.5 1.5 2 3 4 6 12];
hdgs = [-12 -6 -3 -1.5 1.5 3 6 12];

% remove the rest
removethese = ~ismember(data.heading,hdgs);
for F = 1:length(fnames)
    eval(['data.' fnames{F} '(removethese) = [];']);
end

    
%% summary data and plot

dots3DMP_parseData
dots3DMP_plots

%% split choice plots into high and low conf trials

if 0
    
%% fit cumulative gaussians (needed for weights calculation and thresholds)
% and plot them

dots3DMP_fit_cgauss
dots3DMP_plots_cgauss

%% split choice data by high and low confidence

dots3DMP_parseData_splitConf
dots3DMP_plots_splitConf


%% summary statistic for each session

% %correct - overall, and split by high and low confidence
hiCor = data.correct & data.PDW==1 & data.oneConfTargTrial==0;
loCor = data.correct & data.PDW==0 & data.oneConfTargTrial==0;

    
for i=1:length(blocks)
    blocktrs = ismember(data.filename,blocks{i});
    
    PrcCorr(i) = sum(data.correct(blocktrs)) / sum(blocktrs);
    HiCorr(i)  = sum(hiCor(blocktrs)) / sum(data.PDW==1 & blocktrs & data.oneConfTargTrial==0);
    LoCorr(i)  = sum(loCor(blocktrs)) / sum(data.PDW==0 & blocktrs & data.oneConfTargTrial==0);
end
    



%% DDM fitting

options.fitMethod = 'fms';
% options.fitMethod = 'global';
% options.fitMethod = 'multi';
% options.fitMethod = 'pattern';
% options.fitMethod = 'bads';

fixed = [0 1 0 1 1];


% initial guess (or hand-tuned params)
% Tnd is irrelevant for non-RT task
ks = 25;
sigma = 0.025;
B  = 0.8;
Tnd = 250; 
theta = 0.6; % PDW only

guess = [ks sigma B Tnd theta];

% ************************************
% set all fixed to 1 for hand-tuning:
fixed(:)=1;
% (can be used to fix some params and not others)
% ************************************

% plot error trajectory (prob doesn't work with parallel fit methods)
options.ploterr = 0;

options.RTtask = RTtask; % 0 - no need to fit RTs
options.conftask = conftask; % 1 - sacc endpoint, 2 - PDW

[X, err_final, fit, fitInterp] = dots3DMP_fitDDM(data,options,guess,fixed);

% plot it!
dots3DMP_plots_fit(data,fitInterp,conftask,RTtask)
    
end