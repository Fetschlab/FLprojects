% class 'tuned' neurons during dots3DMPtuning paradigm

% e.g. so that we can reference CPs to preference during tuning task

clear;clc;close all
addpath(genpath('/Users/stevenjerjian/Desktop/FetschLab/Analysis/codes/'))

%% Load in the data

subject   = 'lucio';
area      = 'MST';
dateRange = 20220512:20230131;

dataPath = '/Users/stevenjerjian/Desktop/FetschLab/Analysis/data/';
dataFileName = sprintf('%s_%d-%d_neuralData_%s.mat',subject,dateRange(1),dateRange(end),area);
load(fullfile(dataPath,dataFileName));

% 20220512-20230131, remove PIVC and single electrode sessions
removethese = [4 5 8:24 37:41];
dataStruct(removethese) = [];

% % 'clean up' - keep units recorded in both pars
parSelect  = {'dots3DMPtuning','dots3DMP'}; 
minRate    = 5;
minTrs     = [5 10];
dataStruct = dots3DMP_NeuralStruct_runCleanUp(dataStruct,parSelect,minRate,minTrs);

%% define conditions

mods = [1 2 3];
% cohs = [0.2 0.6];
cohs = [1 2]; % use inds instead, real-val cohs were sometimes different
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

%% create FR matrix

useFullDur = 0;
optsMS.smoothFR    = 0;
optsMS.keepSingleTrials = 1;
optsMS.keepMU = 0;

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
eventInfo.binSize = 0.05;

[allUnitsTuning, nUnits] = dots3DMP_FRmatrix_fromDataStruct(dataStruct,'dots3DMPtuning',eventInfo,condsTuning,condlabels,optsMS);



%% statistics

% tests
% 1. compare trial FRs in two time periods (baseline and stim), across all headings (signrank)
% 2. compare trial FRs across headings, within time period (2a. bsln or 2b. stim)

% consider neuron as task-responsive if 1 is significant
% consider neuron as directionally tuned if 2 is significant
% consider neuron as task-mod and dir tuned if 1 & 2 are significant

% TO DO LIST
% validate against rh_plots
% check anova assumptions
% post-hocs?
% switch statistics to R, analysis to Python, viz to Python?

% assumptions of ANOVA
% normality - each sample is drawn from normal dist
% equal var - each population has equal variance
% indep     - each sample per group is independent and observations within
% group are random

[unqModCoh,ia,ic] = unique(condsTuning(:,~contains(condlabels,'heading')),'rows');

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


%% CP analysis (using entire interval mean FRs)

% list of conditions for CPs
[hdg,modality,coh,~,ntr] = dots3DMP_create_trial_list(0,mods,cohs,deltas,1,0); 
condsTask = [modality,coh];
condsTask      = repmat(condsTask,4,1);                                 % for 4 possible choice-wager combinations
condsTask(:,3) = [ones(ntr*2,1); 2*ones(ntr*2,1)];                          % CHOICE
condsTask(:,4) = [zeros(ntr,1); ones(ntr,1); zeros(ntr,1); ones(ntr,1)];    % WAGER
condsTask(:,5) = zeros(size(condsTask,1),1);                                % remove 1-target wagers
condTasklabels = {'modality','coherenceInd','choice','PDW','oneTargConf'};

opts.keepMU = 0; 
opts.keepSingleTrials = 1; % need single trials for CPs

% Task
eventInfo.alignEvent  = {{'stimOn'},{'stimOn','saccOnset'},{'saccOnset'},{'saccOnset','postTargHold'}};
eventInfo.otherEvents = {{''},{''},{''},{''}};
eventInfo.otherEventnames = {{''},{''},{''},{''}};
eventInfo.tStart = [-0.3, 0, -0.3,    0]; 
eventInfo.tEnd   = [   0, 0, -0.1,    0];
eventInfo.binSize = 0.4;
eventInfo.overlap = 0.75;

allUnitsTask = dots3DMP_FRmatrix_fromDataStruct(dataStruct,'dots3DMP',eventInfo,condsTask,condTasklabels,opts);

% now use prefDir from above, & FRtrialmean fields to calculate CP/WPs

[unqModCoh,ia,ic] = unique(condsTask(:,1:2),'rows');

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

            % choice only
            pref_trs = uc==ic & condsTask(:,choiceCol)==prefDir(uc,u);
            null_trs = uc==ic & condsTask(:,choiceCol)~=prefDir(uc,u);

            choiceP(uc,u,iae,1) = rocN(cell2mat(fr_mean(pref_trs)), cell2mat(fr_mean(null_trs)), 100);

            % wager only
            hi_trs = uc==ic & condsTask(:,pdwCol)==1;
            lo_trs = uc==ic & condsTask(:,pdwCol)==0;

            wagerP(uc,u,iae,1) = rocN(cell2mat(fr_mean(hi_trs)), cell2mat(fr_mean(lo_trs)), 100);

            % choice/wager, conditioned
            pref_hi = uc==ic & condsTask(:,choiceCol)==prefDir(uc,u) & condsTask(:,pdwCol)==1;
            null_hi = uc==ic & condsTask(:,choiceCol)~=prefDir(uc,u) & condsTask(:,pdwCol)==1;

            pref_lo = uc==ic & condsTask(:,choiceCol)==prefDir(uc,u) & condsTask(:,pdwCol)==0;
            null_lo = uc==ic & condsTask(:,choiceCol)~=prefDir(uc,u) & condsTask(:,pdwCol)==0;

            % CPs by high and low bets separately
            choiceP(uc,u,iae,2) = rocN(cell2mat(fr_mean(pref_hi)), cell2mat(fr_mean(null_hi)), 100);
            choiceP(uc,u,iae,3) = rocN(cell2mat(fr_mean(pref_lo)), cell2mat(fr_mean(null_lo)), 100);

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
eventInfoTR.binSize = 0.3;
eventInfoTR.overlap = 0.5; % 50% overlap

allUnitsTaskTR = dots3DMP_FRmatrix_fromDataStruct(dataStruct,'dots3DMP',eventInfoTR,condsTask,condTasklabels,opts);

% now use prefDir from above, & FRtrial fields to calculate CP/WPs in bins

[unqModCoh,ia,ic] = unique(condsTask(:,1:2),'rows');

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
% sigTuned = true(size(choiceP,2),1);

% sigTuned = p_bslnStimXhdgs(1,:,2)<0.05 & any(p_bslnStimXhdgs(2:3,:,2)>0.05); % ves-tuned but NOT vis-tuned
% sigTuned = p_bslnStimXhdgs(1,:,2)>0.05 & any(p_bslnStimXhdgs(2:3,:,2)<0.05); % vis-tuned but NOT ves-tuned
sigTuned = p_bslnStimXhdgs(1,:,2)<0.05 & any(p_bslnStimXhdgs(2:3,:,2)<0.05); % ves- and vis-tuned

% sigTuned = any(p_bslnStimXhdgs([1 3 5],:,2)<0.05);
% sigTuned = all(p_bslnStimXhdgs(2:3,:,2)>0.05); 

condnames = {'Ves','Vis-L','Vis-H','Comb-L','Comb-H'};
condcols = 'kmrcb';
spcx = -0.1:0.05:0.1;

pt_titles = {'Unconditioned','Choice-High','Choice-Low'};

% plot High conditioned only, and high coh only
pt = 2;
figure('position',[500 200 600 400],'color','w'); hold on
for uc = [1 3 5]%1:size(choiceP,1)

    yy = choiceP(uc,sigTuned,:,pt);
    errorbar((1:4)+spcx(uc),squeeze(nanmean(yy,2)),squeeze(nanstd(yy,[],2))./sqrt(sum(sigTuned)),'color',condcols(uc),'marker','o','linewidth',2);
    plot([0.5 4.5],[0.5 0.5],'k--')
    set(gca,'xtick',1:4); title(pt_titles{pt})
    set(gca,'xticklabel',{'500ms pre-motion','Stim Period','-300 to -100ms pre-saccade','Saccade to PDW Hold'})
    text(0.8,0.6-0.01*(uc-1),condnames{uc},'color',condcols(uc),'fontsize',16,'fontweight','bold')
    ylabel('Choice Probability')
    changeAxesFontSize(gca,15,15);

end

%%

figure('position',[500 200 1400 400],'color','w')

for pt=1:3
    subplot(1,3,pt); hold on; axis([0.5 4.5 0.4 0.6]);
    for uc = 1:size(choiceP,1)

        yy = choiceP(uc,sigTuned,:,pt);
        errorbar((1:4)+spcx(uc),squeeze(nanmean(yy,2)),squeeze(nanstd(yy,[],2))./sqrt(sum(sigTuned)),'color',condcols(uc),'marker','o','linewidth',1.5);
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

        yy = wagerP(uc,sigTuned,:,pt);
        errorbar((1:4)+spcx(uc),squeeze(nanmean(yy,2)),squeeze(nanstd(yy,[],2))./sqrt(sum(sigTuned)),'color',condcols(uc),'marker','o','linewidth',1.5);
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


subplotInd = [1 3 4 5 6];


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