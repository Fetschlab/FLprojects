% offline analysis wrapper for PLDAPS data, dots3DMP paradigm
% CF started it 10-3-18

clear all; close all

conftask = 1; % 1=colorbars, 2=PDW

normalize = 0;

% specify which (previously saved) mat file to load
 
% %% pre-RT
% subject = 'human';
% paradigm = 'dots3DMP';
% dateRange = 20190612:20191231;
% RTtask = 0;

%% RT

subject = 'human';
paradigm = 'dots3DMP';
dateRange = 20200213:20210526; % RT
RTtask = 1;


%%
folder = '/Users/chris/Documents/MATLAB/PLDAPS_data/';
% folder = '/Users/stevenjerjian/Desktop/FetschLab/Analysis/PLDAPS_data/';
file = [subject '_' num2str(dateRange(1)) '-' num2str(dateRange(end)) '.mat'];
load([folder file], 'data');

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
    data.date(k,1) = str2double(data.filename{k}(9:16));
    data.subjDate{k,:} = [data.subj{k} data.filename{k}(9:16)];
end

%% Some manual excludes e.g. of bad sessions, RT trainings

% consider beter ways of cleaning or labelling sessions beforehand?
% based on notes? add field

% other manual excludes (e.g., RT training)
excludes_filename = {'humanIPQ20200227dots3DMP0904_basic','humanVZC20200229dots3DMP1239_basic'};
excludes_subjDate = {'FRK20200216','FRK20200223','NKT20200215','VZC20200222'};
excludes_subj = {'ASQ', 'HXL', 'XRJ', 'EMF','NEX'};

% excludes_filename = {};
% excludes_subjDate = {};
% excludes_subj = {};

removethese = ismember(data.filename,excludes_filename) | ismember(data.subjDate,excludes_subjDate) | ismember(data.subj,excludes_subj); %#ok<NASGU>
fnames = fieldnames(data);
for F = 1:length(fnames)
    eval(['data.' fnames{F} '(removethese) = [];']);
end
% now this should reflect only good data, per spreadsheet:
blocks = unique(data.filename);


%% choose subjects to include based on 3-letter code

subjs = unique(data.subj); % all

% subjs = {'AAW'}; % the best single subject
% subjs = {'LLV'}; % good shifts at low coh
% subjs = {'CXD'}; % fine, should maybe normalize
% subjs = {'IWT'}; % conf doesn't match large shift at high coh!
% subjs = {'EMF'}; % conf flat, prob discard
% subjs = {'AAW' 'LLV' 'CXD' 'IWT' 'EMF'};
    % ^ the First Five
    
% subjs = {'ASQ'}; % only 62 tr, discard
% subjs = {'HXL'}; % did not do full paradigm, discard
% subjs = {'XRJ'}; % only 31 tr, discard

% subjs = {'AAW' 'LLV' 'CXD'};
% subjs = {'AAW' 'LLV' 'IWT' 'CXD' 'EMF'};

% RT:
% subjs = {'DRH'};
% subjs = {'IPQ'}
% subjs = {'LLV'}
% subjs = {'SJJ'}
% subjs = {'VZC'}

subjs = {'DRH','SJJ','LLV','IPQ','FRK'};

% subjs = {'AAW' 'LLV' 'CXD' 'DRH' 'IPQ' 'SJJ' 'VZC'}; % all 'good' data (pre and post RT)


% remove invalid trials (fixation breaks (which gives nans), excluded subj,
% and obvious testing trials, signaled by very large confidence (saccade
% endpoint) values
removethese = isnan(data.choice) | ~ismember(data.subj,subjs) | data.conf>3 | isnan(data.conf);
fnames = fieldnames(data);
for F = 1:length(fnames)
    eval(['data.' fnames{F} '(removethese) = [];']);
end

% should do the trial number based exclusion here, earlier on we are
% counting fixation breaks

% quick look at blocks, for when some need to be excluded
blocks = unique(data.filename);
nTrialsByBlock = nan(length(blocks),1);
for u = 1:length(blocks)
    nTrialsByBlock(u) = sum(ismember(data.filename,blocks(u)));
end

% we can be pretty sure blocks with <N trials (say, 30) are to be discarded
removethese = ismember(data.filename,blocks(nTrialsByBlock<30));
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

% [C,IA,IC] = unique(data.coherence)
% N = hist(data.coherence,unique(data.coherence))

% cohs = unique(data.coherence); 
    % currently too many cohs in the dataset, so...

% N = hist(data.coherence,unique(data.coherence))'
% unique(data.coherence)


% ...lump coherences together (work in progress)

% % avery
% data.coherence(data.coherence<0.2) = 0.1;
% data.coherence(data.coherence>=0.2 & data.coherence<0.6) = 0.5;
% data.coherence(data.coherence>=0.6) = 0.9;
% cohs = [0.1 0.5 0.9];

% everyone else
data.coherence(data.coherence<=0.5) = 0.4;
% data.coherence(data.coherence>0.1 & data.coherence<=0.4) = 0.2;
data.coherence(data.coherence>0.5) = 0.7;

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
hdgs = unique(data.heading);
hdgs(hdgs==0) = [];

%%

% if RTtask==0
    % same here.  the 1.5-12 range was only used rarely, and in fact is a
      % good signature of warmup or testing-mode trials to be excluded
    %hdgs = [-10 -5 -2.5 -1.25 0 1.25 2.5 5 10]';
    % some zero values were stored as +/- eps in an older version of the gui
    data.heading(abs(data.heading)<0.01) = 0;
% end
    

% remove the rest
removethese = ~ismember(data.heading,hdgs);
for F = 1:length(fnames)
    eval(['data.' fnames{F} '(removethese) = [];']);
end

%% normalize confidence ratings, *within subject*

if normalize

data_orig = data;
usubj = unique(data.subj);
for s = 1:length(usubj)
    data = data_orig;
    removethese = ~strcmp(data.subj,usubj{s});
    for F = 1:length(fnames)
        eval(['data.' fnames{F} '(removethese) = [];']);
    end    
    
    % subtract min and divide by max
    % data.conf = (data.conf - min(data.conf)) / max((data.conf - min(data.conf)));

    % OR
    
    % subtract/divide by *means*
    dots3DMP_parseData
    minOfMeans = nanmin(nanmin(nanmin(nanmin(confMean))));
    data.conf = data.conf - minOfMeans;
    dots3DMP_parseData
    maxOfMeans = nanmax(nanmax(nanmax(nanmax(confMean))));
    data.conf = data.conf / maxOfMeans;

    % OR
    
    % simply cap at 1/0
    % data.conf(data.conf>1) = 1;
    % data.conf(data.conf<0) = 0;

    % append each subj to a new data struct
    if s==1
        data_new = data;
    else
        for F = 1:length(fnames)
            eval(['data_new.' fnames{F} '(end+1:end+length(data.date)) = data.' fnames{F} ';']);
        end
    end
end
data = data_new;

end


    
%% after settling on the above, run this script to generate summary data

dots3DMP_parseData_OLD


%% some plots

% dots3DMP_plots(parsedData,mods,cohs,deltas,hdgs,conftask,RTtask,D)
dots3DMP_plots_OLD


%% fit cumulative gaussians (needed for weights calculation)

dots3DMP_fit_cgauss_OLD

%% and plot them
 
dots3DMP_plots_cgauss




%% nicer looking versions
% 
% dots3DMP_plots_cgauss_forTalk






