
clear;clc;close all
addpath(genpath('/Users/stevenjerjian/Desktop/FetschLab/Analysis/codes/'))

% Load in the data

subject   = 'lucio';
% dateRange = 20220615:20221003;
dateRange = 20220805:20221003;

% dateRange = 20221003;

dataPath = '/Users/stevenjerjian/Desktop/FetschLab/Analysis/data/';
% dataFileName = sprintf('%s_%d-%d_neuralData.mat',subject,dateRange(1),dateRange(end));
dataFileName = sprintf('%s_%d-%d_neuralData_ks25.mat',subject,dateRange(1),dateRange(end));

% dataFileName = 'lucio_20220923_tempneuralData.mat'; 

load(fullfile(dataPath,dataFileName));

dataStruct = dataStruct([2 3 4 7 9 10 11 12 16 17]);

% inputs to dots3DMP_NeuralStruct_runCleanUp
parSelect  = {'dots3DMPtuning','dots3DMP'}; 
minRate    = 5;
minTrs     = [5 10];
dataStruct = dots3DMP_NeuralStruct_runCleanUp(dataStruct,parSelect,minRate,minTrs);

%% pre-set some things

% for tuning functions, unused for now
% tuning_vonMises = @(b,dir) b(1) * exp(b(2)*cosd(dir-b(3))) / (2*pi*besseli(0,b(2))) + b(4);
% tuning_vonMises_err = @(b,dir,FR) nansum((tuning_vonMises(b,dir)-FR).^2); % sum squared error

% these should match across tuning and task
mods = [1 2 3];
% cohs = [0.2 0.6];
cohs = [1 2]; % use inds instead
deltas = 0;

% some sessions for tuning only recorded one
% coherence (typically 'high' i.e. >0.6), but this will be labelled as
% cohInd 1, and in current code set-up, averaging across conditions *across
% neurons* would therefore sometimes average low cohs (in 2 coh session)
% with high cohs from 1-coh sessions...NEED TO FIX THIS before population
% analyses

%% create PSTHs and calculate average activity for units during tuning paradigm
% useFullDur --> stimOn-stimOfyf
% otherwise, middle 1s of stimulus

par     = 'dots3DMPtuning';
% hdgsTuning = [-90 -45 -22.5 -12 0 12 22.5 45 90]; % some early sessions had +/-12 instead
hdgsTuning = [-90 -45 -22.5 0 22.5 45 90];

% generate the whole list of conditions we want to include
[hdg,modality,coh,~,~] = dots3DMP_create_trial_list(hdgsTuning,mods,cohs,deltas,1,0); 
condsTuning = [modality,coh,hdg];
condlabels  = {'modality','coherenceInd','heading'}; % must be consistent order with conds lists below

useFullDur = 0;

optsTR.calcTuning  = 0;
optsTR.smoothFR    = 1;
optsTR.convKernel  = fspecial('average', [1 20]); 

optsMS.calcTuning  = 1;
optsMS.smoothFR    = 0;

if useFullDur
    % use full duration of the stimulus
    eventInfo.alignEvent  = {{'stimOn','stimOff'}}; % align to stimOn
    eventInfo.tStart = 0; % start at stimOn  ...+100ms for sensory delay?
    eventInfo.tEnd   = 0; % end at stimOff
else
    % use the middle 1.5s of the stimulus
    eventInfo.alignEvent  = {{'stimOn','stimOff',0.5}}; % align to middle of stimOn and stimOff
    eventInfo.tStart = -0.75; % 0.75s forward and back
    eventInfo.tEnd   = 0.75;
end

% same for both
eventInfo.otherEvents = {{'fixation','stimOff'}};
eventInfo.otherEventnames = {{'FixAcq','stimOff'}};
eventInfo.binSize = 0.2;

% time-resolved psths, separate alignments to stimOn and saccOnset
TuningeventInfoTR.alignEvent      = {{'stimOn'}};
TuningeventInfoTR.otherEvents     = {{'fixation','stimOff'}};
TuningeventInfoTR.otherEventnames = {{'FixAcq','stimOff'}};
TuningeventInfoTR.tStart          = [-1.5];
TuningeventInfoTR.tEnd            = [3];
TuningeventInfoTR.binSize = 0.01;

allUnitsTuning_MeanStimResp = dots3DMP_FRmatrix_fromDataStruct(dataStruct,par,eventInfo,condsTuning,condlabels,optsMS);
allUnitsTuning_timeResolved = dots3DMP_FRmatrix_fromDataStruct(dataStruct,par,TuningeventInfoTR,condsTuning,condlabels,optsTR);

%% === repeat for Task

par = 'dots3DMP';
hdgsTask = [-12 -6 -3 -1.5 0 1.5 3 6 12];

[hdg,modality,coh,delta,~] = dots3DMP_create_trial_list(hdgsTask,mods,cohs,deltas,1,0);
condsTask = [modality,coh,hdg,ones(size(hdg))];
condlabels  = {'modality','coherenceInd','heading','correct'}; % use only correct trials

% time-resolved psths, separate alignments to stimOn and saccOnset
TaskeventInfoTR.alignEvent      = {{'stimOn'},{'saccOnset'}};
TaskeventInfoTR.otherEvents     = {{'fixation','saccOnset'},{'postTargHold'}};
TaskeventInfoTR.otherEventnames = {{'FixAcq','RT'},{'Wager'}};
TaskeventInfoTR.tStart          = [-1,-0.5];
TaskeventInfoTR.tEnd            = [1.0,1.0];
TaskeventInfoTR.binSize = 0.02;

% mean response during stim duration (i.e. stimOn to saccOnset)
eventInfoMS.alignEvent      = {{'stimOn','saccOnset'}};
eventInfoMS.otherEvents     = {{'fixation','saccOnset','postTargHold'}};
eventInfoMS.otherEventnames = {{'FixAcq','RT','Wager'}};
eventInfoMS.tStart = 0; % start at stimOn + 100ms for sensory delays?
eventInfoMS.tEnd   = 0; % end at stimOff
eventInfoMS.binSize = 0.1;
 
allUnitsTask_MeanStimResp = dots3DMP_FRmatrix_fromDataStruct(dataStruct,par,eventInfoMS,condsTask,condlabels,optsMS);
allUnitsTask_timeResolved = dots3DMP_FRmatrix_fromDataStruct(dataStruct,par,TaskeventInfoTR,condsTask,condlabels,optsTR);


%% ====== plot mean responses during stimulus

% select tuning or Task

par = 'dots3DMPtuning';
% par = 'dots3DMP';

switch par
    case 'dots3DMPtuning'
        auMat   = allUnitsTuning_MeanStimResp;
        condMat = condsTuning;
        hdgVec  = hdgsTuning;

    case 'dots3DMP'
        auMat   = allUnitsTask_MeanStimResp;
        condMat = condsTask;
        hdgVec  = hdgsTask;
end

iae = 1;
auFR_formatted  = reshape(auMat.data.muFRs{iae},length(hdgVec),[],size(auMat.data.muFRs{iae},2));
conds_formatted = reshape(condMat,length(hdgVec),[],size(condMat,2));

condcols = 'kmcrb';
numUnits = size(auFR_formatted,3);

upp = 80; % units per 'page'
[p,q] = numSubplots(upp);    

for u=1:numUnits
   
    iu=u;

    if upp > 16
        fignum = ceil(u/upp); h_fig=figure(fignum);
        set(h_fig,'color','w','position',[100 100 1600 1000]);

        snum = mod(u,upp)+upp*(mod(u,upp)==0);       
    else
        fignum = 1;
        figure(fignum); set(gcf,'color','w','position',[100 100 1400 800]);
        snum = u;
   end

   ax=subplot(p(1),p(2),snum); hold on

   for c = 1:size(conds_formatted,2)
       plot(hdgVec,auFR_formatted(:,c,iu),'color',condcols(conds_formatted(1,c,1)+2*(find(conds_formatted(1,c,2)==cohs)-1)),'linew',1.5,'marker','o')
   end
   set(gca,'XTick',hdgVec);
%    set(gca,'XTickLabel',[]);
   ylim([0 max(10,max(max(auFR_formatted(:,:,iu))))])
   title(sprintf('%d, %d-%d\nunit %d (%s)',iu,auMat.hdr.unitDate(iu),auMat.hdr.unitSet(iu),...
       auMat.hdr.unitID(iu),dataStruct(1).data.(par).units.cluster_labels{auMat.hdr.unitType(iu)}))

end

%% =====  plot mean response during stimulus of one unit

% single tuning curve

iae = 1;
auFR_formatted  = reshape(auMat.data.muFRs{iae},length(hdgVec),[],size(auMat.data.muFRs{iae},2));
seFR_formatted  = reshape(auMat.data.seFRs{iae},length(hdgVec),[],size(auMat.data.seFRs{iae},2));
conds_formatted = reshape(condMat,length(hdgVec),[],size(condMat,2));

condcols = 'kmcrb';
% condcols = 'krb';
legtxt = {'Ves','Vis (Low Coh)','Vis (High Coh)','Comb (Low Coh)','Comb (High Coh)'};

figure;
set(gcf,'color','w','position',[100 100 400 400]);   hold on;

iu = 352;
for c = 1:size(conds_formatted,2)
    errorbar(hdgVec,auFR_formatted(:,c,iu),seFR_formatted(:,c,iu),'color',condcols(conds_formatted(1,c,1)+2*(find(conds_formatted(1,c,2)==cohs)-1)),'linew',2,'marker','o')
end
set(gca,'XTick',hdgVec);
title(sprintf('%d, %d-%d\nunit %d (%s)',iu,auMat.hdr.unitDate(iu),auMat.hdr.unitSet(iu),...
       auMat.hdr.unitID(iu),dataStruct(1).data.(par).units.cluster_labels{auMat.hdr.unitType(iu)}));
changeAxesFontSize(gca,15,15);
hl=legend(legtxt); 
set(hl,'box','off','location','northeast')
% hl.Position(2) = 0.5;
set(gca,'tickdir','out')
ylim([0 1.25 * max(max(auFR_formatted(:,:,iu)))])
ylabel('Firing rate [spikes per sec]','fontsize',14)
xlabel(sprintf('Heading angle [%s]',char(176)))
           


%% ====== plot time resolved PSTH

par = 'dots3DMPtuning';
% par = 'dots3DMP';

switch par
    case 'dots3DMPtuning'
        auMat   = allUnitsTuning_timeResolved;
        condMat = condsTuning;
        hdgVec  = hdgsTuning;

        alignEvent = TuningeventInfoTR.alignEvent;
        otherEventnames = TuningeventInfoTR.otherEventnames;

    case 'dots3DMP'
        auMat   = allUnitsTask_timeResolved;
        condMat = condsTask;
        hdgVec  = hdgsTask;

        alignEvent = TaskeventInfoTR.alignEvent;
        otherEventnames = TaskeventInfoTR.otherEventnames;
end


N = length(hdgVec);
hdgcols = cbrewer('div','RdBu',N*2);
hdgcols = hdgcols([1:floor(N/2) end-floor(N/2):end],:);
hdgcols(hdgVec==0,:) = [0 0 0];

figure('position',[100 100 1000 800],'color','w')

cohs = [1 2];
if length(cohs)==2
    ucond = [1 cohs(1); 2 cohs(1); 2 cohs(2); 3 cohs(1); 3 cohs(2)];
    titles = {'Ves';'Vis (Low Coh)';'Vis (High Coh)';'Comb (Low Coh)';'Comb (High Coh)';'All'};
    subplotInd = [1 3 4 5 6];
else
    ucond = [1 cohs(1); 2 cohs(1); 3 cohs(1)];
    titles = {'Ves';'Vis';'Comb'};
    subplotInd = [1 2 3];
end
% 
% ucond = [1 cohs(1); 2 cohs(1); 2 cohs(2)];
% titles = {'Ves';'Vis (Low Coh)';'Vis (High Coh)'};
% subplotInd = [1 2 3];

% ucond = [1 cohs(1)];
% titles = {'Ves'};

% mcols = {'Greys','Reds','Reds','Blues','Blues','Purples'};

muFRs = auMat.data.PSTHs;

% u = 352;
% u = 292;

u = 58;

clear tempFR
for iae=1:length(muFRs)
    tempFR{iae} = muFRs{iae}(:,:,u);
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

        auFR_formatted  = reshape(muFRs{iae},length(hdgVec),[],size(muFRs{iae},2),size(muFRs{iae},3));
        conds_formatted = reshape(condMat,length(hdgVec),[],size(condMat,2));
        evTimes_formatted = reshape(auMat.times.evTimes{iae},N,[],size(auMat.times.evTimes{iae},2),size(auMat.times.evTimes{iae},3));

        if iae==1, p1 = pos(1);
        else p1 = p1 + pos(3)*sprc(1:iae-1);
        end
        hh=subplot('position',[p1 pos(2) pos(3)*sprc(iae)*0.95 pos(4)]); %axis off;
        hold on;

        temp = squeeze(auFR_formatted(:,c,:,u));
        for h=1:length(hdgVec)
            plot(auMat.times.xvec{iae},temp(h,:),'linew',2,'color',hdgcols(h,:))
        end

        plot([0 0],[0 my],'k','linestyle','-');

        for ioe = 1:size(auMat.times.evTimes{iae},2)
            evMu = mean(evTimes_formatted(:,c,ioe,u)); % what are we averaging over here??

            if evMu<thisTmax && evMu > thisTmin
                plot([evMu evMu],[0 my],'k','linestyle','--');
                text(evMu*0.99,my*1,otherEventnames{iae}{ioe},'fontsize',12,'horizo','right','verti','bottom','rotation',90);
                %             text(evTimes{iae}(c,ioe,u),my*1.02,otherEventnames{iae}{ioe},'fontsize',12,'horizo','center','verti','bottom');
            end
        end

        if ucond(c,1)==max(ucond(:,1))
            xlabel(hh,sprintf('Time from %s [s]',alignEvent{iae}{1}),'fontsize',14)
        end
        if iae==length(muFRs)
            ht=title(titles{c});
            ht.Position(1)=thisTmax;
            ht.HorizontalAlignment = 'right';
        end

        if iae==1
            if ucond(c,2)==cohs(1) && (ucond(c,1)==2 || size(ucond,1)==1)
                ylabel('Firing rate [spikes per sec]','fontsize',14)
            end
        else
            hh.YAxis.Visible = 'off';
        end
        hh.XLim = auMat.times.xvec{iae}([1 end]);


%         if c==1 && iae==length(muFRs)
%             text(thisTmax,my*0.9,sprintf('%d-%d, %s',auMat.hdr.unitDate(u),auMat.hdr.unitID(u),dataStruct(1).data.(par).units.cluster_labels{auMat.hdr.unitType(u)}),...
%                 'horizo','right','fontsize',12);
%         end
        changeAxesFontSize(gca,15,15);
        set(gca,'tickdir','out')

    end

    %         hst = suptitle(sprintf('%d-%d, unit %d (%s)',unitDate(u),unitSet(u),unitID(u),dataStruct(1).data.(par).units.cluster_labels{unitType(u)}));
    %         hst.FontSize = 14;

end

% plot schematic of heading vectors
% spin this out into a neat function at some point

if size(ucond,1)==3
    axes('position',[0.65 0.75 0.17 0.17])
else
    axes('position',[0.55 0.6 0.3 0.3])
end
if max(hdgVec)==12
    sc = 4; len = 4;
elseif max(hdgVec)==90
    sc=1; len = 3;
end

hdgXY = len .* [sind(hdgVec*sc); cosd(hdgVec*sc)];
textVec = hdgXY .* [1.05; 1.1];
startPoint = [0 0];
% axis([-6 6 -6 6]);
axis([-1 1 -1 1]*5); %axis square
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
        text(startPoint(1)+textVec(1,h),startPoint(2)+textVec(2,h),num2str(hdgVec(h)),'fontsize',11,'horizo',ha,'verti',va,'rotation',ra,'color',hdgcols(h,:),'fontweight','bold')
    end
end
set(gca,'Visible','Off');



