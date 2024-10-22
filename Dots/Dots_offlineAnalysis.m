% offline analysis wrapper for PLDAPS data, Dots protocols
% CF spawned it from dots3DMP_offlineAnalysis, begun 10-3-18

%*****************************
% skip this cell to analyze a data struct already in the workspace!
%*****************************

clear; close all

% these will specify which (previously saved) mat file to load
subject = 'hanzo';

folder = '/Users/chris/Documents/MATLAB/PLDAPS_data/';

% all Hanzo good data
% dateRange = 20190401:20201231; % all 'good' data

% Hanzo for R01
% dateRange = 20200728:20200930; % great PDW; misnamed (not really up to 9/30)
dateRange = 20200820:20200922; % nice RT, but PDW kinda flat
% dateRange = 20200901:20200922; % good for both, with TrNo<600 or even 800

% some indiv months
% dateRange = 20190501:20190531; % var dur example
% dateRange = 20201101:20201130; % RT example

% should be best!
% dateRange = 20200801:20201130;

% dateRange = 20210208:20210212; % last week

% maxTrialNum = 700; % set to ~600-800 to omit late trials

if sum(diff(dateRange)>1)==0
    file = [subject '_' num2str(dateRange(1)) '-' num2str(dateRange(end)) '.mat'];
elseif sum(diff(dateRange)>1)==1
    file = [subject '_' num2str(dateRange(1)) '-' num2str(dateRange(diff(dateRange)>1)) '+' num2str(dateRange(find(diff(dateRange)>1)+1)) '-' num2str(dateRange(end)) '.mat'];
else
    file = [subject '_' num2str(dateRange(1)) '---' num2str(dateRange(end)) '.mat'];
end

load([folder file], 'data');


%% more cleanup

% clean up RT variable for mixed RT/non-RT datasets
if isfield(data,'RT')
    data.RT(data.RT>2)=nan; % anything over 2s is almost certainly invalid,
                              % and this also gets rid of the infs
    data.RT(data.RT<0.01)=nan; % zeroes too
else
    data.RT = nan(size(data.choice));
end

% not doing anyhting with these for now:
try
    data = rmfield(data,{'forgivePDW','trainingSeqPDW','confidenceTime'});
catch
end

% remove invalid trials (fixation breaks, etc)
removethese = isnan(data.choice);
try removethese = removethese | data.oneTargPDW==1; catch; end 
fnames = fieldnames(data);
for F = 1:length(fnames)
    eval(['data.' fnames{F} '(removethese) = [];']);
end

% define signed coherence; positive = rightward (0 deg)
data.scoh = data.coherence;
data.scoh(data.direction==180) = -data.scoh(data.direction==180);

% random fixes
data.scoh(abs(data.scoh)<0.001)=0;
if max(data.choice)==2
    data.choice = data.choice-1; % ensure choice is 0 or 1 (sometimes coded as 1..2)
end

% quick look at blocks, for when some need to be excluded
blocks = unique(data.filename);
nTrialsByBlock = nan(length(blocks),1);
data.trialNum = zeros(size(data.filename));
for u = 1:length(blocks)
    nTrialsByBlock(u) = sum(ismember(data.filename,blocks(u)));
    % also number trials within a block, to look at early vs late trials
    data.trialNum(strcmp(blocks(u),data.filename)) = 1:sum(strcmp(blocks(u),data.filename))';
end

% we can be pretty sure blocks with <N trials (say, 10) are to be discarded
removethese = ismember(data.filename,blocks(nTrialsByBlock<10));

% remove early/late trials
if exist('maxTrialNum','var')
    removethese = removethese | data.trialNum>maxTrialNum;
    fnames = fieldnames(data);
    for F = 1:length(fnames)
        eval(['data.' fnames{F} '(removethese) = [];']);
    end
end

% TEMP: also remove non-PDW trials
removethese = isnan(data.PDW) | isnan(data.RT);
fnames = fieldnames(data);
for F = 1:length(fnames)
    eval(['data.' fnames{F} '(removethese) = [];']);
end

%% parse data
conftask = 2; % pdw
RTtask = 1;
RTCorrOnly = 0;
parsedData = Dots_parseData(data,conftask,RTtask,RTCorrOnly);
 


%% some basic stats

% t test on slopes (is sensitivity greater on high-bet trials?):
mu = abs(parsedData.B2(2)-parsedData.B3(2));
se = sqrt(parsedData.stats2.se(2)^2 + parsedData.stats3.se(2)^2);
t = mu/se;
df = sum(~isnan(data.PDW))-length(parsedData.B2); 
pval_ttest_slopes = 2*(1-tcdf(t,df)) % two-tailed

% alternative slope test: indicator version:
t = parsedData.B4(3)/parsedData.stats4.se(3);
df = sum(~isnan(data.PDW))-length(parsedData.B4); 
pval_ttest_slopes_alt = 2*(1-tcdf(t,df)) % two-tailed


% t test on pHigh, corr vs err (is confidence higher on correct trials?)
M = ~isnan(data.PDW) & data.correct==1;
pHighCorr_all = sum(M & data.PDW==1) / sum(M);
MM = ~isnan(data.PDW) & data.correct==0;
pHighErr_all = sum(MM & data.PDW==1) / sum(MM);

pHighSEcorr_all = sqrt( (pHighCorr_all.*(1-pHighCorr_all)) ./ sum(M) );
pHighSEerr_all = sqrt( (pHighErr_all.*(1-pHighErr_all)) ./ sum(MM) );

mu = abs(pHighCorr_all-pHighErr_all);
se = sqrt(pHighSEcorr_all^2 + pHighSEerr_all^2);
t = mu/se; 
df = sum(~isnan(data.PDW))-length(parsedData.B2); 
pval_ttest_pHighCorrErr = 2*(1-tcdf(t,df)) % two-tailed
% [should really use a two-sample Z test for proportions?]



%% plot

forTalk = 0;
cohs = unique(data.scoh);
Dots_plot(parsedData,cohs,conftask,RTtask,0,forTalk);

% %% if var dur, check conf/accuracy vs. dur
% if sum(isnan(data.RT))>0.8*length(data.RT) % arbitrary minimum proportion of non-RT trials; eventually make it a flag
%     Dots_TDA; % time-dependent accuracy function
% end



%% fit simplest DDM

% for RT data, can start here (ignores PDW, but a decent tutorial, based on
% Shadlen et al. 2006 and Palmer et al. 2005)
if sum(isnan(data.RT))<0.8*length(data.RT) % arbitrary minimum proportion of RT trials; eventually make it a flag
    % % params: initial guess
    %     % when stimval is coherence from -0.5..0.5, k is around 0.4
    %     % so to estimate k for a different X variable, say heading, 
    %     % simply scale it down by a factor of max(stimval)
    %     % (this is tricky though because it trades off w the bound)
    % k = max(data.scoh); % 'drift rate' or sensitivity term: a constant converting stimulus
    %                      % strength into units of momentary evidence
    % B = 20; % height of the bound, or threshold, for decision termination
    % Tnd = 225; %?
    

    % or just take vals from sim
    
    k = 0.3; % 'drift rate' or sensitivity term: a constant converting stimulus
             % strength into units of momentary evidence
    B = 25; % height of the bound, or threshold, for decision termination
    Tnd = 300;

    guess = [k B Tnd];
    [b,~] = Dots_fitDDM_1D_noConf(guess,data.scoh,data.choice,round(data.RT*1000));
end



%% fit DDM with confidence (with or without RT)

Dots_fitDDM_wrapper



%% testing out Gabe's code

% convert data struct to Gabe's format:

%   d = 
% 
%   1 x nStimulusConditions struct array with fields:
%     coherence   (1x1, signifies coherence for this stim cond.)
%     RTs         (nTrials x 1, Measured reaction time for each trial in this stim cond.)
%     choice      (nTrials x 1, Measured choice for each trial. 1 = right, 0 = left)
%     nTrials     (1x1, number of trials for this stimulus condition)
% 

u= unique(data.coherence);
for i=1:length(u)
  cohIdx{i}=find(data.coherence==u(i)); %get indeces into cell array
  data2(i).coherence = u(i); % one coherence value for this set of trials
  data2(i).RTs = data.RT(cohIdx{i}); % all of the trials for this coh
  data2(i).choice = data.choice(cohIdx{i}); % all choices for this coh
  data2(i).nTrials = length(cohIdx{i});
end




