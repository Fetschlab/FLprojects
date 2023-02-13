%% dots3DMP behavioral paper

% first demonstration of confidence and RT in self-motion paradigm
% key result is putative common mechanism giving rising to all 3 behavioral
% measures

clear; clc; close all
cd /Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/dataStructs
addpath(genpath('/Users/stevenjerjian/Desktop/FetschLab/Analysis/codes/'))

printFigs = 0;

%% Load in relevant data

subject = 'human'; conftask = 1; RTtask   = 1;      % human RT
% subject = 'human'; conftask = 1; RTtask   = 0;    % human non-RT (fix dur)
% subject = 'lucio'; conftask = 2; RTtask   = 1;    % monkey 1 (PDW + RT)
% subject = 'zarya'; conftask = 2; RTtask   = 1;    % monkey 2 (PDW + RT)

data = dots3DMP_loadBehaviorData(subject,conftask,RTtask);


if strcmp(subject, 'human')
    data.oneTargConf = false(size(data.heading)); % move to cleanUP

    subjs2keep = [1 2 3 4 5];
    subjs = unique(data.subj);
    
    fnames = fieldnames(data);
    removethese = ~ismember(data.subj,subjs(subjs2keep));
    for f=1:length(fnames), data.(fnames{f})(removethese) = []; end

    % remove zero heading trials, insufficient data
    removethese = data.heading==0;
    for f=1:length(fnames), data.(fnames{f})(removethese) = []; end

end



% if doing pooled analysis, then just go forward with data
% otherwise re-assign data to be from a given subject and proceed the same
% way (except maybe we want to use raw conf values?) (so re-assign confRaw
% to conf)
% data.confNorm = data.conf;
% data.conf = data.confRaw;

% to do

% - review each subjects individual curves
% - compare means to mean errors, particularly for cue conflict
% - visualization of accuracy, conf, RT relationships, as below
% - regression analyses
%       - effect of combined on accuracy
%       - conf to accuracy
%       - cue conflict to conf, 

% 2D Acc model fit
% alternative models - intuition + maths



%% basic data parsing and visualization

mods   = unique(data.modality); 
cohs   = unique(data.coherence); 
deltas = unique(data.delta);
hdgs   = unique(data.heading);

% deltas = [-3 3]; % for showing only non-zero deltas

% means per condition, logistic and gaussian fits
parsedData  = dots3DMP_parseData(data,mods,cohs,deltas,hdgs,conftask,RTtask); 
gfit        = dots3DMP_fit_cgauss(data,mods,cohs,deltas,conftask,RTtask); 

% plotting
fig_logistic = dots3DMP_plots(parsedData,mods,cohs,deltas,hdgs,conftask,RTtask); % logistic fits
fig_gaussFits = dots3DMP_plots_cgauss_byCoh(gfit,parsedData,mods,cohs,deltas,hdgs,conftask,RTtask); % eventually will just go straight to plotting fits of model?
    
%% 1. relationship between confidence and accuracy

% 1. conf vs P(right) for each modality
fig_ConfvsAcc = dots3DMP_plots_confchoice(parsedData,cohs,deltas,conftask);
 
%% 2. relationship between confidence and RT

% TODO review confidence distributions within each subject to determine
% threshold for high/lo (currently just the median normalized value)?

% for human data, 

% D = find(deltas==0); 
D = length(deltas)+1;

parsedData_byConf = dots3DMP_parseData_byConf(data,mods,cohs,deltas,hdgs,conftask,RTtask); 
gfit_byConf       = dots3DMP_fit_cgauss_byConf(data,mods,cohs,deltas,conftask,RTtask,D);
dots3DMP_plots_cgauss_byConf(gfit_byConf,parsedData_byConf,mods,2,deltas,hdgs,conftask,RTtask,D) % just 1 coh

%% 2b. RT quantiles
% 0 - plot errors/low bet only, 1 - plot correct/high bet only, 2 - plot correct/error or high/low bet separately, -1 - plot all trials

if RTtask
    dots3DMP_RTquantiles(data,conftask,2); 
%     confRT_distrs(data,mods,cohs,conftask,RTtask)
end

%% meta-sensitivity

% how good are subjects at reflecting on their own confidence
% across all
[d,meta_d] = dots3DMP_dprime_measures(data,mods,cohs,deltas); % NEED THE OUTPUTS!
%
% 
% 
% %% psychophysical cue weights
% % assume wvis always = 1-wves
% 
% wves = dots3DMP_cueWeights(gfit,cohs,deltas,conftask);
% 
% % bootstrapping for error bars, set rng for reproducibility
% rng(28);
% nboots = 20;
% [gfitBoot,wvesBoot] = dots3DMP_cgauss_bootstrap_func(data,gfit,mods,cohs,deltas,nboots,conftask,RTtask);
% dots3DMP_plotCueWeights(wves,wvesBoot,cohs,conftask) % plot the weights


%% individual days (NHP) or subjects (human)

switch subject
    case 'human'
        uF = data.subj;
    case {'lucio','zarya'}
        uF = data.subjDate;
end
        
uss = unique(uF);
fnames = fieldnames(data);

mods   = unique(data.modality); 
cohs   = unique(data.coherence); 
deltas = unique(data.delta);
hdgs   = unique(data.heading);

clear parsedData_bySess gfit_bySess
for s = 1:length(uss)
    temp = data;
    for f = 1:length(fnames)
        temp.(fnames{f})(~strcmp(uF,uss{s})) = [];
    end

    parsedData_bySess(s) = dots3DMP_parseData(temp,mods,cohs,deltas,hdgs,conftask,RTtask);
    gfit_bySess(s)       = dots3DMP_fit_cgauss(temp,mods,cohs,deltas,conftask,RTtask);
%     f1           = dots3DMP_plots_cgauss_byCoh(gfit_bySess(s),parsedData_bySess(s),mods,cohs,deltas,hdgs,conftask,RTtask); 
end


%% plot 1 subject/session

s=5;
dots3DMP_plots_cgauss_byCoh(gfit_bySess(s),parsedData(s),mods,cohs,deltas,hdgs,conftask,RTtask); 


%% 







%% 
% this should be another analysis function
% plot subject cue conflict shifts

clear shift error
% examine absolute size of shifts
for s = 1:length(subjs)
    for c = 1:length(cohs)
    shift(s,c,1) = (gfit(s).choice.mu(3,c,1)-gfit(s).choice.mu(3,c,2));
    shift(s,c,2) = (gfit(s).choice.mu(3,c,3)-gfit(s).choice.mu(3,c,2));
    
%     relshift(s,c) = (gfit(s).choice.mu(3,c,1)-gfit(s).choice.mu(3,c,3));
    
%     shift(s,c,1) = gfit(s).choice.mu(3,c,1);
%     shift(s,c,2) = gfit(s).choice.mu(3,c,3);
    
    error(s,c,1) = (gfit(s).choice.muSE(3,c,1));
    error(s,c,2) = (gfit(s).choice.muSE(3,c,3));
    end
end

shift2 = reshape(shift,size(shift,1),[]);
error2 = reshape(error,size(error,1),[]);
%%
figure;
alim = 3;
ssz = 70;
mkr = 'osd^p<>h';
subplot(211); axis square; 
for s = 1:length(subjs)
    scatter(abs(shift2(s,:)),error2(s,:),ssz,'k',mkr(s)); hold on
end
hr=refline(1,0); set(hr,'linew',1.5,'linestyle','--','color','k')
axis([0 alim 0 alim])
xlabel('|mu|')
ylabel('muSE')
tidyaxes; 
changeAxesFontSize(gca,14,14);

subplot(212); hold on; 
for s = 1:length(subjs)
    scatter(shift(s,1,1),error(s,1,1),ssz,'b',mkr(s));
    scatter(shift(s,1,2),error(s,1,2),ssz,'g',mkr(s));
    scatter(shift(s,2,1),error(s,2,1),ssz,'b',mkr(s),'filled');
    scatter(shift(s,2,2),error(s,2,2),ssz,'g',mkr(s),'filled');
end
xlabel('mu')
ylabel('muSE')
hr=refline(1,0); set(hr,'linew',1.5,'linestyle','--','color','k')
axis([-1 1 -1 1]*alim)
tidyaxes; 
changeAxesFontSize(gca,14,14);
prettyfig;
set(gcf,'position',[20 15 11 18]);





