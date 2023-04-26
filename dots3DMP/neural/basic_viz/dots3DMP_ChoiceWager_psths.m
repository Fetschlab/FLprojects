% choice-wager conditioned PSTHs

clear;clc;close all
addpath(genpath('/Users/stevenjerjian/Desktop/FetschLab/Analysis/codes/'))

% Load in the data

subject   = 'lucio';
dateRange = 20220512:20230131;

dataPath = '/Users/stevenjerjian/Desktop/FetschLab/Analysis/data/';
dataFileName = sprintf('%s_%d-%d_neuralData.mat',subject,dateRange(1),dateRange(end));
load(fullfile(dataPath,dataFileName));

% 20220512-20230131, remove 20220520
dataStruct(4) = [];


load(fullfile(dataPath,dataFileName));

% inputs to dots3DMP_NeuralStruct_runCleanUp
parSelect  = {'dots3DMPtuning','dots3DMP'}; 
minRate    = 5;
minTrs     = [5 10];
dataStruct = dots3DMP_NeuralStruct_runCleanUp(dataStruct,parSelect,minRate,minTrs);

%% get PSTHs for each neuron for each choice and wager outcome

% need to find a way to code option to include e.g. multiple headings, but
% not separate by them for the conditions

% ALSO NEED simple way to plot only selected mods/cohs without having to
% rerun this part

optsTR.smoothFR    = 1;
optsTR.convKernel  = fspecial('average', [1 20]); 

mods = [1 2 3];
cohs = [1 2];

par = 'dots3DMP';
% hdgs = [-12 -6 -3 -1.5 0 1.5 3 6 12]; 
hdgs = [0]; % just use zero heading!
% or eventually, a few small headings, but need to subtract mean and
% collapse

[hdg,modality,coh,delta,nconds] = dots3DMP_create_trial_list(hdgs,mods,cohs,0,1,0);
condsTask = [modality,coh,hdg]; 
condsTask(:,end+1) = 0; % 2-target confidence only!
nCol = size(condsTask,2);

% choice and wager
condsTask = repmat(condsTask,4,1);
condsTask(1:nconds*2,nCol+1)     = 1; 
condsTask(nconds*2+1:end,nCol+1) = 2;
condsTask([1:nconds nconds*2+1:nconds*2+nconds],nCol+2)   = 1;
condsTask([nconds+1:nconds*2 nconds*2+nconds+1:end],nCol+2) = 0; % this is not really necessary...
ltxt = {'L-hi','L-lo','R-hi','R-lo'};

% choice only
% condsTask = repmat(condsTask,2,1);
% condsTask(1:nconds,nCol+1) = 1;
% condsTask(nconds+1:end,nCol+1) = 2;
% ltxt = {'L','R'};

condlabels  = {'modality','coherenceInd','heading','oneTargConf','choice','PDW'}; 
% condlabels  = {'modality','coherenceInd','heading','choice'}; 


% time-resolved psths, separate alignments to stimOn and saccOnset
TaskeventInfoTR.alignEvent      = {{'stimOn'},{'saccOnset'}};
TaskeventInfoTR.otherEvents     = {{'fixation','saccOnset'},{'postTargHold'}};
TaskeventInfoTR.otherEventnames = {{'FixAcq','RT'},{'Wager'}};
TaskeventInfoTR.tStart          = [-1,-0.5];
TaskeventInfoTR.tEnd            = [1.0,1];
TaskeventInfoTR.binSize = 0.02;

allUnitsChoiceWager = dots3DMP_FRmatrix_fromDataStruct(dataStruct,par,TaskeventInfoTR,condsTask,condlabels,optsTR);

% nconds = length(unique([modality,coh],'rows'));
%%


% in dots3DMP task

% 14, 15, 20, 28, 36, 42, 44, 46, 47, 58, 60, 61, 75, 90, 92, 105, 106, 125, 151, 153, 157, 158, 159, 161,
% 163, 166,184, 187, 207, 217, 226, 233, 236, 242, 247, 257, 260, 262, 263, 264, 265, 272, 275, 278, 282, 284

% 263
fsz = 20;

u = 197;
auMat = allUnitsChoiceWager;
muFRs = auMat.data.PSTHs;

alignEvent = TaskeventInfoTR.alignEvent;
otherEventnames = TaskeventInfoTR.otherEventnames;

cols = cbrewer('qual','Paired',6);
% cols = cols(end:-1:end-3,:);
cols = cols([6 5 2 1],:);

% cols = [1 0 0; 0 0 1];

% figure('position',[100 100 1000 800],'color','w')
figure('position',[100 100 900 500],'color','w')

if length(cohs)==2
    ucond = [1 cohs(1); 2 cohs(1); 2 cohs(2); 3 cohs(1); 3 cohs(2)];
    titles = {'Ves';'Vis (Low Coh)';'Vis (High Coh)';'Comb (Low Coh)';'Comb (High Coh)';'All'};
    subplotInd = [1 3 4 5 6];
else
    ucond = [1 cohs(1); 2 cohs(1); 3 cohs(1)];
    titles = {'Ves';'Vis';'Comb'};
    subplotInd = [1 2 3];
end

[uc,ia,ic] = unique(condsTask(:,1:3),'rows');


clear tempFR
for iae=1:length(muFRs)
    tempFR{iae} = muFRs{iae}(:,:,u);
% 
%     for c = 1:length(uc)
%         theseFR = tempFR{iae}(ic==c,:);
%         tempFR{iae}(ic==c,:) = theseFR - nanmean(theseFR,1); 
%     end
% 
%     dmFRs{iae}(:,:,u) = tempFR{iae};
end

mxFR = cellfun(@(x) max(x(:)), tempFR);
my   = max(mxFR)*1.2;

xlens = cellfun(@range,auMat.times.xvec);
sprc  = xlens./sum(xlens);

for c=1:size(ucond,1)
%     axMain = gca;
    axMain = subplot(3,2,subplotInd(c));
    pos=get(axMain,'Position');

    for iae=1:length(muFRs)
        thisTmax = auMat.times.xvec{iae}(end);
        thisTmin = auMat.times.xvec{iae}(1);

        auFR_formatted  = reshape(muFRs{iae},nconds,[],size(muFRs{iae},2),size(muFRs{iae},3));
        conds_formatted = reshape(condsTask,nconds,[],size(condsTask,2));
        evTimes_formatted = reshape(auMat.times.evTimes_byUnit{iae},nconds,[],size(auMat.times.evTimes_byUnit{iae},2),size(auMat.times.evTimes_byUnit{iae},3));

        if iae==1, p1 = pos(1);
        else p1 = p1 + pos(3)*sprc(1:iae-1);
        end
        hh=subplot('position',[p1 pos(2)+0.05 pos(3)*sprc(iae)*0.95 pos(4)-0.05]); %axis off;
        hold on;

        temp = squeeze(auFR_formatted(c,:,:,u));
        for h=1:size(temp,1)
            plot(auMat.times.xvec{iae},temp(h,:),'linew',2,'color',cols(h,:))
        end

        hh.XTick = [fliplr(-0.5:-0.5:thisTmin) 0:0.5:thisTmax];
        hh.XTickLabelRotation = 0;
        plot([0 0],[0 my],'k','linestyle','-','linewidth',2);

        for ioe = 1:size(auMat.times.evTimes_byUnit{iae},2)
            evMu = mean(evTimes_formatted(c,:,ioe,u)); % average over conditions

            if evMu<thisTmax && evMu > thisTmin
                plot([evMu evMu],[0 my],'k','linestyle','--','linewidth',2);
                text(evMu*0.99,my*1,otherEventnames{iae}{ioe},'fontsize',fsz,'horizo','right','verti','bottom','rotation',90);
                %             text(evTimes{iae}(c,ioe,u),my*1.02,otherEventnames{iae}{ioe},'fontsize',12,'horizo','center','verti','bottom');
            end
        end

        if ucond(c,1)==max(ucond(:,1))
            xlabel(hh,sprintf('Time from %s [s]',alignEvent{iae}{1}),'fontsize',fsz)
        end
        if iae==length(muFRs)
            ht=title(titles{c});
            ht.Position(1)=thisTmax;
            ht.HorizontalAlignment = 'right';
%             lh=legend(ltxt,'fontsize',fsz);
%             lx.Box = 'off';

            
%             inset = axes('position',[0.6 0.6 0.25 0.2]);
%             axis(inset,[-1 1 -1 1]); axis square; hold on;
%             scatter(inset,[-1 -1 1 1],[1 -1 1 -1]*0.6,300,cols,'filled')
%             scatter(inset,0,0,150,'k','filled')
%             inset.YAxis.Visible = 'off';
%             inset.XAxis.Visible = 'off';
%             inset.Title.String = 'PDW Targets';
%             inset.Title.fontsize = fsz;
        end

        if iae==1
%             if ucond(c,2)==cohs(1) && (ucond(c,1)==2 || size(ucond,1)==1)
                ylabel(hh,'Firing rate (sp/s)','fontsize',fsz)
%             end
        else
            hh.YAxis.Visible = 'off';
        end
        hh.XLim = auMat.times.xvec{iae}([1 end]);
        hh.YLim = [0 my];

%         if c==1 && iae==length(muFRs)
%             text(thisTmax,my*0.9,sprintf('%d-%d, %s',auMat.hdr.unitDate(u),auMat.hdr.unitID(u),dataStruct(1).data.(par).units.cluster_labels{auMat.hdr.unitType(u)}),...
%                 'horizo','right','fontsize',12);
%         end
        changeAxesFontSize(hh,fsz,fsz);
        set(hh,'tickdir','out')
        
    end

    %         hst = suptitle(sprintf('%d-%d, unit %d (%s)',unitDate(u),unitSet(u),unitID(u),dataStruct(1).data.(par).units.cluster_labels{unitType(u)}));
    %         hst.FontSize = 14;

end

%%
% inset = axes('position',[0.8 0.7 0.1 0.15]);
figure('color','w')
axis([-1 1 -1 1]); axis square; hold on;
scatter([-1 -1 1 1],[1 -1 1 -1]*0.6,500,cols,'filled')
scatter(0,0,200,'k','filled')
set(get(gca,'YAxis'),'Visible','off');
set(get(gca,'XAxis'),'Visible','off');
title('PDW targets','fontsize',24)
