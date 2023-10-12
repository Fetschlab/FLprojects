% class 'tuned' neurons during dots3DMPtuning paradigm
% e.g. so that we can reference CPs to preference during tuning task

clear; clc

%% Load data

subject   = 'lucio';
dateRange = 20220512:20230602;

dataPath = '/Users/stevenjerjian/Desktop/FetschLab/Analysis/data/lucio_neuro_datasets';
dataFileName = sprintf('%s_%d-%d_neuralData.mat',subject,dateRange(1),dateRange(end));
load(fullfile(dataPath,dataFileName));

%% 'clean up' - keep units recorded in both pars
% crude diagnostic criteria could be revised
parSelect  = {'dots3DMPtuning','dots3DMP'}; 
minRate    = 2;
minTrs     = [5 8];
dataStruct = dots3DMP_NeuralStruct_runCleanUp(dataStruct,parSelect,minRate,minTrs);

%% define conditions

mods = [1 2 3];
% cohs = [0.2 0.6];
% cohs = [1 2]; % use inds instead, real-val cohs were sometimes different
cohs = 1; % collapse across cohs later anyway
deltas = 0;

% use actual headings vector, and then remove/ignore NaNs
% should get this vector from matching behavioral dataset eventually
% or generate behavioral dataset from neural dataStruct
hdgsTuning = [-90 -60 -45 -25 -22.5 -12 0 12 22.5 25 45 60 90];
hdgsTuning([1 end]) = []; % drop +/- 90, weird platform motion

% generate the whole list of conditions we want to include
[hdg,modality,coh,~,~] = dots3DMP_create_trial_list(hdgsTuning,mods,cohs,deltas,1,0); 
condsTuning = [modality,coh,hdg];
condlabels  = {'modality','coherenceInd','heading'}; % must be consistent order with conds lists below

% collapse across cohs
condsTuning = [modality, hdg];
condlabels = {'modality','heading'};

%% create FR matrix

useFullDur = 0;
optsMS.smoothFR = 0;
optsMS.convKernel = fspecial('gaussian', [1 40], 5);
optsMS.keepSingleTrials = 1; % needed for tuning statistics 
optsMS.keepMU = 1;

if useFullDur
    % use full duration of the stimulus
    eventInfo.alignEvent  = {{'stimOn','stimOff'}}; % align to stimOn
    eventInfo.tStart = 0; % start at stimOn  ...+100ms for sensory delay?
    eventInfo.tEnd   = 0; % end at stimOff
else
    % use the middle of the stimulus
    eventInfo.alignEvent  = {{'fixation'},{'stimOn','stimOff',0.5}}; % align to middle of stimOn and stimOff
    eventInfo.tStart = [-0.3, -0.5]; % forward and back
    eventInfo.tEnd   = [0.1, 0.5];
end

eventInfo.otherEvents = {{'stimOn'},{'fixation'}};
eventInfo.otherEventnames = {{'stimOn'},{'FixAcq'}};
eventInfo.binSize = 0.2;

[allUnitsTuning, ~, nUnits] = dots3DMP_FRmatrix_fromDataStruct(dataStruct,'dots3DMPtuning',eventInfo,condsTuning,condlabels,optsMS);


%% statistics

% tests
% 1. compare trial FRs in two time periods (baseline and stim), across all headings (signrank)
% 2. compare trial FRs across headings, within time period (2a. bsln or 2b. stim)

% consider neuron as task-responsive if 1 is significant
% consider neuron as directionally tuned if 2 is significant
% consider neuron as task-mod and dir tuned if 1 & 2 are significant

% TO DO LIST
% validate/sanity-check against rh_plots
% check anova assumptions
% post-hocs?


% assumptions of ANOVA
% normality - each sample is drawn from normal dist
% equal var - each population has equal variance
% indep     - each sample per group is independent and observations within
% group are random

[unqModCoh,ia,ic] = unique(condsTuning(:,contains(condlabels,{'modality','coherence'})),'rows');

p_bslnVSstim     = nan(size(unqModCoh,1),nUnits);       % test 1
p_bslnStimXhdgs  = nan(size(unqModCoh,1),nUnits,2);     % test 2a,2b
p_bslnStimXhdgs_perm = nan(size(unqModCoh,1),nUnits,2); % test 2a,2b

p_bslnStimHdgs_liltest = p_bslnStimXhdgs;
p_bslnStimHdgs_vartest = p_bslnStimXhdgs;

prefDir = nan(size(unqModCoh,1),nUnits);
prefAmp = nan(size(unqModCoh,1),nUnits);


for u = 1:nUnits

    for uc = 1:size(unqModCoh,1)

        mean_stimFR = allUnitsTuning.data.FRmean{2}(uc==ic,u);

        bslnFR = cell2mat(allUnitsTuning.data.FRtrialmean{1}(uc==ic,u));
        stimFR = cell2mat(allUnitsTuning.data.FRtrialmean{2}(uc==ic,u));

        if isempty(bslnFR) || isempty(stimFR), continue, end

        % create headings vector to correspond to trial frs
        nTrs_per_hdg = allUnitsTuning.data.condntrs{2}(uc==ic,u);

        thisHdg = condsTuning(uc==ic,contains(condlabels,'heading'));
        thisHdg = repelem(thisHdg,nTrs_per_hdg);

        % compare firing rates to baseline, across all headings
        [p_bslnVSstim(uc,u),h] = signrank(bslnFR,stimFR);

%         % check for normally distributed samples
%         [~, p_bslnStimHdgs_liltest(uc,u,1)] = lillietest(bslnFR);
%         [~, p_bslnStimHdgs_liltest(uc,u,2)] = lillietest(stimFR);
% 
%         % check for equal variances
%         p_bslnStimHdgs_vartest(uc,u,1) = vartestn(bslnFR,thisHdg,'display','off');
%         p_bslnStimHdgs_vartest(uc,u,2) = vartestn(stimFR,thisHdg,'display','off');

        % compare across headings
        p_bslnStimXhdgs(uc,u,1) = anova1(bslnFR,thisHdg,'off');
        p_bslnStimXhdgs(uc,u,2) = anova1(stimFR,thisHdg,'off');

        % permutation ANOVA, gives almost exactly the same results
%         p_bslnStimXhdgs_perm(uc,u,1) = randanova1(bslnFR,thisHdg);
%         p_bslnStimXhdgs_perm(uc,u,2) = randanova1(stimFR,thisHdg);

        % heading preference direction and magnitude (R or L)

                                        % R>L + 1 means right pref will be 2, left pref will be 1, like behavior             
        prefDir(uc,u) = double( nanmean(mean_stimFR(hdgsTuning>0)) > nanmean(mean_stimFR(hdgsTuning<0)) ) + 1;

                        % (R-L) / (R+L)
        prefAmp(uc,u) = ( nanmean(mean_stimFR(hdgsTuning>0)) - nanmean(mean_stimFR(hdgsTuning<0) ) ./ nanmean(mean_stimFR(hdgsTuning~=0)));


    end
    
end


% impute pref dir in combined condition if not available
% if it prefers different directions in two modalities, assume prefDir in
% comb matches
% whichever unisensory preference is stronger

% within coh?

if any(contains(condlabels, 'coherence'))
    
end

[~, prefmod] = max(prefAmp);
nancomb = isnan(prefDir(3,:));

prefDir_imputed = prefDir;
prefDir_imputed(3,nancomb) = prefDir(sub2ind(size(prefDir),prefmod(nancomb),find(nancomb)));

prefDir = prefDir_imputed;

%% CP analysis (using entire interval mean FRs)

% make this a function with inputs e.g. "conditions_matrix"

% list of conditions for CPs

[hdg,modality,coh,~,ntr] = dots3DMP_create_trial_list(0,mods,1,deltas,1,0); 
condsTask = [modality,coh,hdg];
condsTask      = repmat(condsTask,4,1);                                 % for 4 possible choice-wager combinations
condsTask(:,4) = [ones(ntr*2,1); 2*ones(ntr*2,1)];                          % CHOICE
condsTask(:,5) = [zeros(ntr,1); ones(ntr,1); zeros(ntr,1); ones(ntr,1)];    % WAGER
condsTask(:,6) = zeros(size(condsTask,1),1);                                % remove 1-target wagers
condTasklabels = {'modality','coherenceInd','heading', 'choice','PDW','oneTargConf'};

condsTask(:,2)=[];
condTasklabels(2)=[];

clear opts
opts.keepMU = 1; 
opts.keepSingleTrials = 1; % need single trials for CPs
%opts.smoothFR = 1;
%optsMS.convKernel = fspecial('gaussian', [1 40], 5);

% Task
clear eventInfo
eventInfo.alignEvent  = {{'stimOn'},{'stimOn','saccOnset'},{'saccOnset'},{'saccOnset'}};
eventInfo.otherEvents = {{''},{''},{''},{''}};
eventInfo.otherEventnames = {{''},{''},{''},{''}};
eventInfo.tStart = [-0.3, 0, -0.3,    -0.1]; 
eventInfo.tEnd   = [   0, 0, -0.1,    +0.1];
eventInfo.binSize = 0.2;
% eventInfo.overlap = 0.75;

allUnitsTask = dots3DMP_FRmatrix_fromDataStruct(dataStruct,'dots3DMP',eventInfo,condsTask,condTasklabels,opts);

% now use prefDir from above, & FRtrialmean fields to calculate CP/WPs

[unqModCoh,ia,ic] = unique(condsTask(:,contains(condTasklabels,{'modality','coherence'})),'rows');

choiceCol = strcmpi(condTasklabels,'choice');
pdwCol    = strcmpi(condTasklabels,'pdw') | strcmpi(condTasklabels,'wager');

choiceP = nan(size(unqModCoh,1),nUnits,length(eventInfo.alignEvent),3);
wagerP  = nan(size(unqModCoh,1),nUnits,length(eventInfo.alignEvent),3);

for u = 1:nUnits

    for iae = 1:length(eventInfo.alignEvent)

        for uc = 1:size(unqModCoh,1)

            % skip if condition wasn't present
            if isnan(prefDir(uc,u)), continue, end

            fr_mean = allUnitsTask.data.FRtrialmean{iae}(:,u);

            if ~sum(uc==ic & ~cellfun(@isempty,fr_mean))
                continue
            end

            % choice only
            pref_trs = uc==ic & condsTask(:,choiceCol)==prefDir(uc,u);
            null_trs = uc==ic & condsTask(:,choiceCol)~=prefDir(uc,u);

            try
                choiceP(uc,u,iae,1) = rocN(cell2mat(fr_mean(pref_trs)), cell2mat(fr_mean(null_trs)), 100);
            end

            % wager only
            hi_trs = uc==ic & condsTask(:,pdwCol)==1;
            lo_trs = uc==ic & condsTask(:,pdwCol)==0;

            try
                wagerP(uc,u,iae,1) = rocN(cell2mat(fr_mean(hi_trs)), cell2mat(fr_mean(lo_trs)), 100);
            end

            % choice/wager, conditioned
            pref_hi = uc==ic & condsTask(:,choiceCol)==prefDir(uc,u) & condsTask(:,pdwCol)==1;
            null_hi = uc==ic & condsTask(:,choiceCol)~=prefDir(uc,u) & condsTask(:,pdwCol)==1;

            pref_lo = uc==ic & condsTask(:,choiceCol)==prefDir(uc,u) & condsTask(:,pdwCol)==0;
            null_lo = uc==ic & condsTask(:,choiceCol)~=prefDir(uc,u) & condsTask(:,pdwCol)==0;

            try
                % CPs by high and low bets separately
                choiceP(uc,u,iae,2) = rocN(cell2mat(fr_mean(pref_hi)), cell2mat(fr_mean(null_hi)), 100);
                choiceP(uc,u,iae,3) = rocN(cell2mat(fr_mean(pref_lo)), cell2mat(fr_mean(null_lo)), 100);
            end

            if sum(pref_hi) && sum(pref_lo) && sum(null_hi) && sum(null_lo)
                % WPs by pref and null dir separately
                wagerP(uc,u,iae,2) = rocN(cell2mat(fr_mean(pref_hi)), cell2mat(fr_mean(pref_lo)), 100);
                wagerP(uc,u,iae,3) = rocN(cell2mat(fr_mean(null_hi)), cell2mat(fr_mean(null_lo)), 100);
            end

        end
    end
end


%% time-resolved version

% Task
eventInfoTR.alignEvent  = {{'stimOn'},{'saccOnset'}};
eventInfoTR.otherEvents = {{'fixation'},{'postTargHold'}};
eventInfoTR.otherEventnames = {{'FixAcq'},{'WagerHold'}};
eventInfoTR.tStart = [-0.6 -0.3]; 
eventInfoTR.tEnd   = [+0.9 0.6];
eventInfoTR.binSize = 0.01;

opts.keepMU = 1; 
opts.keepSingleTrials = 1; % need single trials for CPs
opts.smoothFR = 1;
opts.convKernel = fspecial('gaussian', [1 40], 5);


allUnitsTaskTR = dots3DMP_FRmatrix_fromDataStruct(dataStruct,'dots3DMP',eventInfoTR,condsTask,condTasklabels,opts);

% now use prefDir from above, & FRtrial fields to calculate CP/WPs in bins
[unqModCoh,ia,ic] = unique(condsTask(:,contains(condlabels,{'modality','coherence'})),'rows');

choiceCol = strcmpi(condTasklabels,'choice');
pdwCol    = strcmpi(condTasklabels,'pdw') | strcmpi(condTasklabels,'wager');

choiceP_TR = cell(length(eventInfoTR.alignEvent),1);
wagerP_TR  = cell(length(eventInfoTR.alignEvent),1);


for iae = 1:length(eventInfoTR.alignEvent)

    nBins = length(allUnitsTaskTR.times.xvec{iae});

    choiceP_TR{iae} = nan(size(unqModCoh,1),nUnits,nBins,3);
    wagerP_TR{iae} = nan(size(unqModCoh,1),nUnits,nBins,3);

    for u = 1:nUnits


        for uc = 1:size(unqModCoh,1)

            % skip, didn't have this condition
            if isnan(prefDir(uc,u)), continue, end

            fr_binned = allUnitsTaskTR.data.FRtrial{iae}(:,u);

            % choice only
            pref_trs = uc==ic & condsTask(:,choiceCol)==prefDir(uc,u);
            null_trs = uc==ic & condsTask(:,choiceCol)~=prefDir(uc,u);

            % wager only
            hi_trs = uc==ic & condsTask(:,pdwCol)==1;
            lo_trs = uc==ic & condsTask(:,pdwCol)==0;

            % choice/wager, conditioned
            pref_hi = uc==ic & condsTask(:,choiceCol)==prefDir(uc,u) & condsTask(:,pdwCol)==1;
            null_hi = uc==ic & condsTask(:,choiceCol)~=prefDir(uc,u) & condsTask(:,pdwCol)==1;

            pref_lo = uc==ic & condsTask(:,choiceCol)==prefDir(uc,u) & condsTask(:,pdwCol)==0;
            null_lo = uc==ic & condsTask(:,choiceCol)~=prefDir(uc,u) & condsTask(:,pdwCol)==0;

            pref_un_frs = cell2mat(fr_binned(pref_trs));
            null_un_frs = cell2mat(fr_binned(null_trs));

            hi_un_frs = cell2mat(fr_binned(hi_trs));
            lo_un_frs = cell2mat(fr_binned(lo_trs));

            pref_hi_frs = cell2mat(fr_binned(pref_hi));
            null_hi_frs = cell2mat(fr_binned(null_hi));

            pref_lo_frs = cell2mat(fr_binned(pref_lo));
            null_lo_frs = cell2mat(fr_binned(null_lo));

  
            % loop over bins
            for t = 1:size(pref_un_frs,2)

                choiceP_TR{iae}(uc,u,t,1) = rocN(pref_un_frs(:,t),null_un_frs(:,t), 100);

                wagerP_TR{iae}(uc,u,t,1) = rocN(hi_un_frs(:,t),lo_un_frs(:,t), 100);


                % CPs by high and low bets separately
                choiceP_TR{iae}(uc,u,t,2) = rocN(pref_hi_frs(:,t),null_hi_frs(:,t), 100);
                choiceP_TR{iae}(uc,u,t,3) = rocN(pref_lo_frs(:,t),null_lo_frs(:,t), 100);

                if sum(pref_hi) && sum(pref_lo) && sum(null_hi) && sum(null_lo)
                    % WPs by pref and null dir separately
                    wagerP_TR{iae}(uc,u,t,2) = rocN(pref_hi_frs(:,t),pref_lo_frs(:,t), 100);
                    wagerP_TR{iae}(uc,u,t,3) = rocN(null_hi_frs(:,t),null_lo_frs(:,t), 100);
                end
            end

        end
    end
end


% save TuningStats_w_ChoiceWagerProbabilities choiceP wagerP p_bslnStimXhdgs p_bslnVSstim prefDir prefAmp -mat



%% PLOT EPOCH-CPs

% nansum(p_bslnStimXhdgs(:,:,2)<0.05,2)./sum(~isnan(p_bslnStimXhdgs(:,:,2)),2)
% sum(~isnan(p_bslnStimXhdgs(:,:,2)),2)


% sigTuned = any(p_bslnStimXhdgs(:,:,2)<0.05,1);
% sigTuned = all(p_bslnStimXhdgs(:,:,2)<0.05,1); 
% sigTuned = true(1, size(choiceP,2));

% sigTuned = p_bslnStimXhdgs(1,:,2)<0.05 & any(p_bslnStimXhdgs(2:3,:,2)>0.05); % ves-tuned but NOT vis-tuned
% sigTuned = p_bslnStimXhdgs(1,:,2)>0.05 & any(p_bslnStimXhdgs(2:3,:,2)<0.05); % vis-tuned but NOT ves-tuned
% sigTuned = p_bslnStimXhdgs(1,:,2)<0.05 & any(p_bslnStimXhdgs(2:3,:,2)<0.05); % ves- and vis-tuned

inds = [1 2 3];
sigTuned = any(p_bslnStimXhdgs(inds,:,2)<0.05);
% sigTuned = all(p_bslnStimXhdgs(2:3,:,2)>0.05); 

areaMST = strcmp(allUnitsTask.hdr.area, 'MSTd');

selUnits = sigTuned & areaMST;

condnames = {'Ves','Vis-L','Vis-H','Comb-L','Comb-H'};
condcols = 'kmrcb';
condnames = {'Ves','Vis','Comb'};
condcols = 'krb';
spcx = -0.1:0.05:0.1;

%%  plot High conditioned only, and high coh only
pt = 2;
figure('position',[500 200 600 400],'color','w'); hold on
for uc = [1 2 3]%1:size(choiceP,1)

    yy = choiceP(uc,selUnits,:,pt);
    errorbar((1:4)+spcx(uc),squeeze(nanmean(yy,2)),squeeze(nanstd(yy,[],2))./sqrt(sum(selUnits)),'color',condcols(uc),'marker','o','linewidth',2);
    plot([0.5 4.5],[0.5 0.5],'k--')
    set(gca,'xtick',1:4); title(pt_titles{pt})
    set(gca,'xticklabel',{'500ms pre-motion','Stim Period','-300 to -100ms pre-saccade','Saccade to PDW Hold'})
    text(0.8,0.6-0.01*(uc-1),condnames{uc},'color',condcols(uc),'fontsize',16,'fontweight','bold')
    ylabel('Choice Probability')
    changeAxesFontSize(gca,15,15);

end

%% plot all

pt_titles = {'Unconditioned','Choice-High','Choice-Low'};

figure('position',[500 200 1400 400],'color','w')

for pt=1:3
    subplot(1,3,pt); hold on; axis([0.5 4.5 0.4 0.6]);
    for uc = 1:size(choiceP,1)

        yy = choiceP(uc,selUnits,:,pt);
        errorbar((1:4)+spcx(uc),squeeze(nanmean(yy,2)),squeeze(nanstd(yy,[],2))./sqrt(sum(selUnits)),'color',condcols(uc),'marker','o','linewidth',1.5);
        plot([0.5 4.5],[0.5 0.5],'k--')
        set(gca,'xtick',1:4); title(pt_titles{pt})
        set(gca,'xticklabel',{'500ms pre-motion','Stim Period','-300 to -100ms pre-saccade','Saccade to PDW Hold'})
        if pt>1
            set(gca,'yticklabel',[])
        else
            text(0.8,0.6-0.01*(uc-1),condnames{uc},'color',condcols(uc),'fontsize',12,'fontweight','bold')
            ylabel('Choice Probability')
        end
        changeAxesFontSize(gca,15,15);

    end
end


pt_titles = {'Unconditioned','Wager-Pref','Wager-Null'};

figure('position',[500 800 1400 400],'color','w')

for pt=1:3
    subplot(1,3,pt); hold on; axis([0.5 4.5 0.4 0.6]);
    for uc = 1:size(wagerP,1)

        yy = wagerP(uc,selUnits,:,pt);
        errorbar((1:4)+spcx(uc),squeeze(nanmean(yy,2)),squeeze(nanstd(yy,[],2))./sqrt(sum(selUnits)),'color',condcols(uc),'marker','o','linewidth',1.5);
        plot([0.5 4.5],[0.5 0.5],'k--')
        set(gca,'xtick',1:4); title(pt_titles{pt})
        set(gca,'xticklabel',{'500ms pre-motion','Stim Period','-300 to -100ms pre-saccade','Saccade to PDW Hold'})
            
        if pt>1
            set(gca,'yticklabel',[])
        else
            text(0.8,0.6-0.01*(uc-1),condnames{uc},'color',condcols(uc),'fontsize',12,'fontweight','bold')
            ylabel('Wager Probability')
        end
        
        changeAxesFontSize(gca,15,15);
    end
end


% Figs

% Epoch_CPs/WPs_1 02/14/2023
% sigTuned defined as SUs with any condition 'tuned' for heading (62SUs)

% Epoch_CPs/WPs_all 02/14
% as above, but all SUs (109)


%% PLOT TIME-RESOLVED VERSION

sigTuned = any(p_bslnStimXhdgs(:,:,2)<0.05,1);
% sigTuned = true(size(choiceP,2),1);


% subplotInd = [1 3 4 5 6];
subplotInd = [1 2 3];


% relative subplot widths for alignments
xlens = cellfun(@range,allUnitsTaskTR.times.xvec);
sprc  = xlens./sum(xlens);

cols = 'kbg';

figure('position',[500 200 1000 800],'color','w')

for uc = 1:size(unqModCoh,1)

    axMain = subplot(3,2,subplotInd(uc));
    pos=get(axMain,'Position');

    for iae=1:length(choiceP_TR)
        
        t = allUnitsTaskTR.times.xvec{iae};

        if iae==1, p1 = pos(1);
        else, p1 = p1 + pos(3)*sprc(1:iae-1);
        end
        hh=subplot('position',[p1 pos(2)+0.05 pos(3)*sprc(iae)*0.95 pos(4)-0.05]); %axis off;
        hold on;

        %hh.XTick = [fliplr(-0.5:-0.5:t(1)) 0:0.5:t(end)];
        hh.XTickLabelRotation = 0;
        plot([0 0],[0 1],'k','linestyle','-','linewidth',2);

%         for ioe = 1:size(auMat.times.evTimes_byUnit{iae},2)
%             evMu = mean(evTimes_formatted(:,c,ioe,u)); % what are we averaging over here??
% 
%             if evMu<thisTmax && evMu > thisTmin
%                 plot([evMu evMu],[0 my],'k','linestyle','--','linewidth',2);
%                 if c==1
%                     text(evMu*0.99,my*1,otherEventnames{iae}{ioe},'fontsize',fsz,'horizo','right','verti','bottom','rotation',90);
%                 %             text(evTimes{iae}(c,ioe,u),my*1.02,otherEventnames{iae}{ioe},'fontsize',12,'horizo','center','verti','bottom');
%                 end
%             end
%         end

        xlabel(hh,sprintf('Time from %s [s]',eventInfoTR.alignEvent{iae}{1}),'fontsize',14)

        if iae==1
                %ylabel('Firing rate (sp/s)','fontsize',14)
        else
            hh.YAxis.Visible = 'off';
        end
        hh.XLim = t([1 end]);
        hh.YLim = [0.45 0.55];


        % Plot the actual means

        meanCP = squeeze(nanmean(choiceP_TR{iae}(uc,:,:,:),2));
        semCP  = squeeze(nanstd(choiceP_TR{iae}(uc,:,:,:),[],2))./sqrt(sum(sigTuned));


%         meanCP = squeeze(nanmean(wagerP_TR{iae}(uc,:,:,:),2));
%         semCP  = squeeze(nanstd(choiceP_TR{iae}(uc,:,:,:),[],2))./sqrt(sum(sigTuned));

        for r=1:size(meanCP,2)
            plot(t,meanCP(:,r),'color',cols(r),'linew',2);
        end

        changeAxesFontSize(gca,14,14);
        set(gca,'tickdir','out')

    end

end


%%  ======= conditioned PSTH averages =======

% list of conditions for conditioned PSTHs
% choice/wager
hdgs = [0];
[hdg,modality,coh,~,ntr] = dots3DMP_create_trial_list(hdgs,mods,1,deltas,1,0); % 0 heading
condsTask = [modality,coh,hdg];
condsTask      = repmat(condsTask,4,1);                                 % for 4 possible choice-wager combinations
condsTask(:,4) = [ones(ntr*2,1); 2*ones(ntr*2,1)];                          % CHOICE
condsTask(:,5) = [zeros(ntr,1); ones(ntr,1); zeros(ntr,1); ones(ntr,1)];    % WAGER
condsTask(:,6) = zeros(size(condsTask,1),1);                                % remove 1-target wagers
condTasklabels = {'modality','coherenceInd','heading','choice','PDW','oneTargConf'};
% condTasklabels = {'modality','coherenceInd','choice','PDW','oneTargConf'};

% ignore coh
condsTask(:,2) = [];
condTasklabels(2) = [];

% optsTR.collapse_conds = [0 1 1 0 0 0]; % collapse across hdgs
optsTR.collapse_conds = [0 1 0 0 0];


%%
% heading conditions
hdgs = [-12 -6 -3 -1.5 0 1.5 3 6 12];
[hdg,modality,coh,~,ntr] = dots3DMP_create_trial_list(hdgs,mods,1,deltas,1,0); % 0 heading
condsTask = [modality,coh,hdg];
condsTask      = repmat(condsTask,2,1);                                 
condsTask(:,end+1) = [zeros(ntr,1); ones(ntr,1)];                          % correct
% condsTask(:,end+1) = [ones(ntr,1); 2*ones(ntr,1)];                          % choice
condTasklabels = {'modality','coherenceInd','heading','correct'};

% ignore coh
condsTask(:,2) = [];
condTasklabels(2) = [];

optsTR.collapse_conds = [0 0 0];


%%
eventInfoTR.alignLabel      = { 'Baseline' , 'stimulus on',             'saccade'};  % labels, for plotting
eventInfoTR.alignEvent      = {{'fixation'}, {'motionOn'} ,             {'saccOnset'}};
eventInfoTR.otherEvents     = {{'targsOn'} , {'fixation','saccOnset'} , {'postTargHold'}};
eventInfoTR.otherEventnames = {{'targsOn'}  , {'FixAcq', 'RT'}  ,       {'WagerHold'}}; % labels, for plotting
eventInfoTR.tStart          = [-0.5             -0.7                    -0.3]; 
eventInfoTR.tEnd            = [0                1.2                     1.0];
eventInfoTR.binSize         = 0.001; 

optsTR.smoothFR = 1;
optsTR.convKernel = fspecial('gaussian', [1 150], 25);
% NOTE: conv with gaussian is really slow!

[allUnitsTaskTR, condsTask, ~] = dots3DMP_FRmatrix_fromDataStruct(dataStruct,'dots3DMP',eventInfoTR,condsTask,condTasklabels,optsTR);


%% processing of PSTHs

% zscore against baseline, and recode for pref/null direction preference

% bsln is 1st
% mean across conditions, std across conditions in 'baseline' interval

% mean across time first, then mean or sd across conditions
% could use FRmean, but then have to deal with dims. Future problem :)
mu = nanmean(nanmean(allUnitsTaskTR.data.PSTHs{1}, 1), 2); % mean for each unit independently
sd = nanstd(nanmean(allUnitsTaskTR.data.PSTHs{1}, 1), [], 2); % std

zscored_psths = cellfun(@(x) (x-mu)./sd, allUnitsTaskTR.data.PSTHs, 'UniformOutput', false);

% recode psths to be pref/null instead of right/left

condTasklabels_collapsed = condTasklabels(~optsTR.collapse_conds);
% condTasklabels_collapsed = condTasklabels;

[unqModCoh,ia,ic] = unique(condsTask(:,contains(condTasklabels_collapsed,{'modality','coherence'})),'rows');
choiceCol = strcmpi(condTasklabels_collapsed,'choice');
% choiceCol = strcmpi(condTasklabels_collapsed,'correct');


recoded_psths = cell(size(zscored_psths));
for iae = 1:length(recoded_psths)
    recoded_psths{iae} = nan(size(zscored_psths{iae}));

    for u = 1:nUnits

        for uc = 1:size(unqModCoh,1)

            tempFR = zscored_psths{iae}(uc==ic, :, u);
            newFR  = nan(size(tempFR));
            temp_conds = condsTask(uc==ic, :);

            if prefDir(uc,u)==1

                % for wagering plots
                newFR(temp_conds(:,choiceCol)==2, :) = tempFR(temp_conds(:,choiceCol)==1, :);
                newFR(temp_conds(:,choiceCol)==1, :) = tempFR(temp_conds(:,choiceCol)==2, :);


                % for heading conditions
%                 newFR(temp_conds(:,2)==0, :) = tempFR(temp_conds(:,2)==0, :);
%                 newFR(sign(temp_conds(:,2))==1, :) = tempFR(sign(temp_conds(:,2))==-1, :);
%                 newFR(sign(temp_conds(:,2))==-1, :) = tempFR(sign(temp_conds(:,2))==1, :);

            else
                newFR = tempFR;
            end

            recoded_psths{iae}(uc==ic, :, u) = newFR;

        end
    end

end


%% plot average conditioned PSTH for example groups of conditions

% TODO
% zscore across ALL headings?
% collapse cohs
% split cells by 'area', MST or PIVC
% add relevant 'otherevents' to plots
% somehow flip sign of all neurons so that prefDir is always rightward?

ymax = 3.5;
condnames = {'Ves','Vis','Comb'};

areaSel = strcmp(allUnitsTaskTR.hdr.area, 'MSTd');
% areaSel = strcmp(allUnitsTaskTR.hdr.area, 'PIVC');

% sigTuned = any(p_bslnStimXhdgs(:,:,2)<0.05,1); % any cond significant
% sigTuned = all(p_bslnStimXhdgs(:,:,2)<0.05,1); % all conds sig
sigTuned = true(size(areaSel)); % ignore tuning sig

% sigTuned = p_bslnStimXhdgs(1,:,2)<0.05 & p_bslnStimXhdgs(2,:,2)>0.05; % ves-tuned but NOT vis-tuned
% sigTuned = p_bslnStimXhdgs(1,:,2)>0.05 & p_bslnStimXhdgs(2,:,2)<0.05; % vis-tuned but NOT ves-tuned
% sigTuned = p_bslnStimXhdgs(1,:,2)<0.05 & p_bslnStimXhdgs(2,:,2)<0.05; % ves- and vis-tuned

inds = [1 2 3];
% sigTuned = any(p_bslnStimXhdgs(inds,:,2)<0.05);
% sigTuned = all(p_bslnStimXhdgs(2:3,:,2)>0.05); 

selUnits = areaSel & sigTuned;

% for choice/wager (blue green shades, like targets in task)
hdgcols = [0.65 0.80 0.90;
        0.10 0.45 0.70;
        0.70 0.90 0.55;
        0.20 0.60 0.20];

hdgcols = cbrewer('qual','Set2',4);

% for headings, spectrum
% N = length(hdgs);
% cols = cbrewer('div','RdBu',N*2);
% cols = cols([1:floor(N/2) end-floor(N/2):end],:);
% cols(hdgs==0,:) = [0 0 0];

% hdgcols = repmat(cols, size(condsTask,1)/size(cols,1), 1);

% subplotInd = [1 3 4 5 6];
subplotInd = [1 2 3];

ls = '-';

% relative subplot widths for alignments
al_inds = 2:3;
xlens = cellfun(@range,allUnitsTaskTR.times.xvec(al_inds));
sprc(al_inds)  = xlens./sum(xlens);

figure('position',[500 200 700 700],'color','w')

for uc = 1:size(unqModCoh,1)

    
    axMain = subplot(3,1,subplotInd(uc));
%     axMain = subplot(3,2,subplotInd(uc));
    pos=get(axMain,'Position');

    for iae=al_inds
        
        t = allUnitsTaskTR.times.xvec{iae};

        if iae==al_inds(1), p1 = pos(1);
        else, p1 = p1 + pos(3)*sprc(iae-1);
        end
        hh=subplot('position',[p1 pos(2)+0.05 pos(3)*sprc(iae)*0.95 pos(4)-0.05]); %axis off;
        hold on;

        %hh.XTick = [fliplr(-0.5:-0.5:t(1)) 0:0.5:t(end)];
        hh.XTickLabelRotation = 0;

        if uc==size(unqModCoh,1)
            xlabel(hh,{'Time from', sprintf('%s [s]',eventInfoTR.alignEvent{iae}{1})},'fontsize',14)
        end

        if iae==al_inds(1)
            if uc==2
                %ylabel('Firing rate (sp/s)','fontsize',14)
                ylabel('Normalized population activity','fontsize', 14)
            end
        else
            hh.YAxis.Visible = 'off';
        end
        hh.XLim = t([1 end]);
        hh.YLim = [-1.5 ymax];

        nTrs_conds = allUnitsTaskTR.data.condntrs{iae}(ic==uc, :);
        enough_trials = all(nTrs_conds>=3);

        oevs = allUnitsTaskTR.times.evTimes_bySession{iae};
        oevs = nanmean(nanmean(oevs(ic==uc, :, :), 3), 1);

        % Plot the actual means
        temp_psth = recoded_psths{iae}(ic==uc,:,selUnits & enough_trials);
%         temp_psth = zscored_psths{iae}(ic==uc,:,selUnits);

        temp_conds = condsTask(ic==uc,:);

        for cc = 1:size(temp_psth, 1)
            %if cc>9, ls = '-'; else, ls = ':'; end
            %if ismember(cc,[1:3 7:9]), continue, end
            %if cc<=9, continue, end
            plot(t, nanmean(temp_psth(cc,:,:), 3), 'color', hdgcols(cc,:),'linew',1.5, 'linestyle',ls)
        end
        ylim = get(hh, 'ylim');

        line([oevs; oevs],repmat(ylim', 1, length(oevs)),'color','k','linestyle','--')

        for ioe = 1:length(oevs)
            if oevs(ioe) > t(1) && oevs(ioe) < t(end)
                text(oevs(ioe),ylim(2),eventInfoTR.otherEventnames{iae}{ioe},'verti','bottom','horizo','center')
            end
        end

        plot([0 0],ylim,'k','linestyle','-','linewidth',2);
        if iae==al_inds(1)
            if uc==2
            lh=legend('null-lo','null-hi','pref-lo','pref-hi');
            set(lh,'box','off')
            end
        elseif iae==al_inds(end)
            ht=title(condnames{uc});
            ht.Position(1) = ht.Position(1)+0.4;
        end
    

        changeAxesFontSize(gca,14,14);
        set(gca,'tickdir','out')

    end

end

%%
fsz=12;
hdgVec = hdgs;
figure('position',[1000 100 600 600],'color','w')
if max(hdgVec)==12
    sc = 3; len = 3;
elseif max(hdgVec)==90
    sc=1; len = 3;
end

hdgXY = len .* [sind(hdgVec*sc); cosd(hdgVec*sc)];
textVec = hdgXY .* [1.05; 1.1];
startPoint = [0 0];
% axis([-6 6 -6 6]);
axis([-1 1 -1 1]*4); axis square
hold on
for h=1:length(hdgVec)
    plot(startPoint(1)+[0; hdgXY(1,h)],startPoint(2)+[0; hdgXY(2,h)],'linew',2,'color',hdgcols(h,:));

    [th,r] = cart2pol(startPoint(1)+textVec(1,:),startPoint(2)+textVec(2,:));
    th = rad2deg(th);
    if hdgVec(h)>0, ha = 'left'; ra = th(h); va = 'middle';
    elseif hdgVec(h)<0, ha = 'right'; ra = th(h)-180; va = 'middle';
    else , ha = 'center'; ra = 0; va = 'bottom';
    end
    if hdgVec(h)==0 || abs(hdgVec(h))>2
        text(startPoint(1)+textVec(1,h),startPoint(2)+textVec(2,h),num2str(hdgVec(h)),'fontsize',fsz,'horizo',ha,'verti',va,'rotation',ra,'color',hdgcols(h,:),'fontweight','bold')
    end
end
set(gca,'Visible','Off');
