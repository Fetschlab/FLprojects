%% dots3DMP human paper

% first, make sure we have cleaned datasets

% human non-RT
% human RT

% clean means no breakfix trials, only 2-target trials, only 'good'
% subjects, coherence and headings standardized, normalized conf
% see dots3DMP_cleanHumanData
% (keeping raw cohs, headings and conf values as well)

clear; clc; close all
cd /Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/dataStructs
addpath(genpath('/Users/stevenjerjian/Desktop/FetschLab/Analysis/codes/'))

%%

printFigs = 0;

subject = 'human';
conftask = 1;
RTtask   = 1; 


if ~RTtask
    load('human_20190625-20191231_nonRT_clean_Apr2022.mat');% human non-RT
    subjs2keep = [1 2 3 4 5];

else
    load('human_20200213-20220317_RT_clean_Apr2022.mat') % human RT
    subjs2keep = [1 2 4 5];
    
%     data.Confnorm = data.conf;
%     data.conf = data.confRaw;
    
end

subjs = unique(data.subj);
fnames = fieldnames(data);
removethese = ~ismember(data.subj,subjs(subjs2keep));
for f=1:length(fnames), data.(fnames{f})(removethese) = []; end

% remove zero heading trials, insufficient data
% fnames = fieldnames(data);
% for f=1:length(fnames), data.(fnames{f})(data.heading==0) = []; end


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



%% basic data parsing and visualization

% these will be useful in multiple places
mods   = unique(data.modality); 
cohs   = unique(data.coherence); 
deltas = unique(data.delta);
% deltas = [-3 3]; % for showing only non-zero deltas
hdgs   = unique(data.heading);

% means per condition, logistic and gaussian fits
parsedData  = dots3DMP_parseData(data,mods,cohs,deltas,hdgs,conftask,RTtask); 
gfit        = dots3DMP_fit_cgauss(data,mods,cohs,deltas,conftask,RTtask); 
% dots3DMP_plots(parsedData,mods,cohs,deltas,hdgs,conftask,RTtask)
f1=dots3DMP_plots_cgauss_byCoh(gfit,parsedData,mods,cohs,deltas,hdgs,conftask,RTtask); % eventually will just go straight to plotting fits of model?
if printFigs
    printfig(f1(1),'Behavior_AllSubj_Modalities','pdf');
    printfig(f1(2),'Behavior_AllSubj_CueConflict','pdf');
end
    
%% 1. relationship between confidence and accuracy

% 1. conf vs P(right) for each modality
f2=dots3DMP_plots_confchoice(parsedData,cohs,deltas,conftask);
if printFigs, printfig(f2,'ConfVsPRight','pdf'), end

%% 2. relationship between confidence and RT

% TODO review confidence distributions within each subject to determine
% threshold for high/lo (currently just the median normalized value)?

parsedData_byConf = dots3DMP_parseData_byConf(data,mods,cohs,deltas,hdgs,conftask,RTtask); 
gfit_byConf       = dots3DMP_fit_cgauss_byConf(data,mods,cohs,deltas,conftask,RTtask);
dots3DMP_plots_cgauss_byConf(gfit_byConf,parsedData_byConf,mods,cohs,deltas,hdgs,conftask,RTtask)

% 0 - plot errors/low bet only, 1 - plot correct/high bet only, 2 - plot correct/error or high/low bet separately, -1 - plot all trials
if RTtask
    dots3DMP_RTquantiles(data,conftask,1); 
%     confRT_distrs(data,mods,cohs,conftask,RTtask)
end
if printFigs
    printfig(fh(1),'RTquantiles_Conf','pdf');
    printfig(fh(2),'RTquantiles_Acc','pdf')
end


% Technical TODO
% make sure each figure is numbered and sized and saves, then make a
% wrapper that does analyses with each subject's data, by sub-sampling
% relevant trials from data



%% meta-sensitivity

% how good are subjects at reflecting on their own confidence
% across all
% dots3DMP_dprime_measures(data,mods,cohs,deltas) % NEED THE OUTPUTS!
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


%% individual subjects psychometric curve, and cue conflict shift/SE comparisons

subjs = unique(data.subj);
fnames = fieldnames(data);

mods   = unique(data.modality); 
cohs   = unique(data.coherence); 
deltas = unique(data.delta);
hdgs   = unique(data.heading);

clear parsedData gfit
for s = 1:length(subjs)
    temp = data;
    for f = 1:length(fnames)
        temp.(fnames{f})(~strcmp(data.subj,subjs{s})) = [];
    end
    
    
    parsedData(s)= dots3DMP_parseData(temp,mods,cohs,deltas,hdgs,conftask,RTtask);
    gfit(s)      = dots3DMP_fit_cgauss(temp,mods,cohs,deltas,conftask,RTtask);
%     f1           = dots3DMP_plots_cgauss_byCoh(gfit(s),parsedData(s),mods,cohs,deltas,hdgs,conftask,RTtask); 
%     set(f1(1),'Units','centimeters'); set(f2(1),'Units','centimeters');
%     printfig(f1(1),sprintf('Behavior_IndivSubjs_Modalities_%s',subjs{s}),'pdf');
%     printfig(f1(2),sprintf('Behavior_IndivSubjs_CueConflict_%s',subjs{s}),'pdf');
end

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





