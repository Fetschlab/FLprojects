% SJ 08/2021

subject = 'zarya';
paradigm = 'dots3DMP';
dateRange = 20210813:20210825;

dateStr = num2str(dateRange(1));
for d = 2:length(dateRange)
    dateStr = [dateStr newline num2str(dateRange(d))];
end

localDir = ['/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/' subject '/nexonar/'];
remoteDir = ['/var/services/homes/fetschlab/data/' subject '/' subject '_nexonar/'];
% remoteDir = ['/var/services/homes/fetschlab/data/'];

subject = 'test';
useVPN = 0;
overwriteLocalFiles = 0; % set to 1 to always use the server copy
getDataFromServer % now also includes pdsCleanup to reduce file size and complexity

%%
localDir = ['/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/' subject '/'];
remoteDir = ['/var/services/homes/fetschlab/data/' subject '/' ];

% localDir = ['/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/'];
% remoteDir = ['/var/services/homes/fetschlab/data/'];

getDataFromServer % now also includes pdsCleanup to reduce file size and complexity


%%
clear all

subject = 'zarya';
localDirNex = ['/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/' subject '/nexonar/'];
localDirPDS = ['/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/' subject '/'];

filename = 'zarya20210817dots3DMP1538';
load(fullfile(localDirNex,[filename '_nexonar.mat']))
load(fullfile(localDirPDS,[filename '.mat']))

%%
data = createDataStructure_oneFile(PDS,subject);

PDS.conditions(end) = [];
PDS.data(end) = [];
% clean-up Nex data, UDP issue 08-2021
[nexPDS,nexClean] = dots3DMP_nexonarCleanUp(nex,PDS);

data.nexonar = nexPDS';

%%

% load('nexClean_SJtest_20210825.mat');
mods = [1 2 3];
% hdgs = [-12 -6 -3 -1.5 0 1.5 3 6 12];
hdgs = unique(data.heading);

hdgCol = flipud(cbrewer('div','RdBu',length(hdgs)));

% set small hdg to zero, change color of zero heading to mark it easily
data.heading(abs(data.heading)<0.01) = 0;
hdgs(abs(hdgs)<0.01) = 0;
hdgCol(hdgs==0,:) = [0.5 0.5 0.5];

maxlen = max(cellfun(@(x) size(x,1), nex.nexdata));

% nexDat = nex.nexdata;
nexDat = data.nexonar;
nexMat = nan(maxlen,size(nexDat{1},2),length(nexDat));

for t=1:length(nexDat)
    nexMat(1:size(nexDat{t},1),:,t) = nexDat{t};
end

nexMat = permute(nexMat,[1 3 2]);
mu = nexMat(1,:,:);
nexMat = bsxfun(@minus,nexMat,mu);
% nexMat(:,:,3) = -nexMat(:,:,3);
%%
nT = 200; %


figure(1);
clear hh
for m=1:length(mods)
    subplot(2,3,m); hold on;
%     axis([0 5000 -40 40]); % X pos vs time
%     axis([-40 40 0 180])
    for h=1:length(hdgs)
%         trs = nex.behavior.goodtrial & nex.conditions.heading==hdgs(h) & nex.conditions.modality==mods(m) & nex.conditions.delta==0;
        trs = ~isnan(data.choice) & data.heading==hdgs(h) & data.modality==mods(m) & data.delta==0;
%         plot(nexMat(:,trs,1),nexMat(:,trs,3),'color',hdgCol(h,:),'marker','o','linestyle','none');
%         plot(nexMat(:,trs,1),nexMat(:,trs,3),'color',hdgCol(h,:),'linewidth',0.5);
        plot(nexMat(1:nT,trs,3),'color',hdgCol(h,:),'linewidth',0.5);
%         plot(nexMat(:,trs,3),nexMat(:,trs,4),'color',hdgCol(h,:),'linewidth',0.5);
        
%         hold on;
%         plot(squeeze(mean(nexMat(:,trs,3),2)),squeeze(mean(nexMat(:,trs,4),2)),'color',hdgCol(h,:),'linewidth',1.5);
 
    end
end
%%
modlabels = {'Ves','Vis','Comb'};
figure('position',[100 100 800 800],'color','w');
for m=1:length(mods)
    for h=1:length(hdgs)
        trs = ~isnan(data.choice) & data.heading==hdgs(h) & data.modality==mods(m) & data.delta==0;
        subplot(3,3,m); hold on; axis([-40 40 0 180]); offsetAxes; title(sprintf('%s,x-y',modlabels{m}));
        plot(nexMat(1:nT,trs,3),nexMat(1:nT,trs,4),'color',hdgCol(h,:),'linewidth',0.5);
        subplot(3,3,m+3); hold on; axis([-40 40 0 5]); offsetAxes; title(sprintf('%s,x-z',modlabels{m}));
        plot(nexMat(1:nT,trs,3),nexMat(1:nT,trs,5),'color',hdgCol(h,:),'linewidth',0.5);
        subplot(3,3,m+6); hold on; axis([0 180 -2 5]); offsetAxes; title(sprintf('%s,y-z',modlabels{m}));
        plot(nexMat(1:nT,trs,4),nexMat(1:nT,trs,5),'color',hdgCol(h,:),'linewidth',0.5);
    end
end

figure('position',[100 100 800 800],'color','w');
for m=1:length(mods)
    for h=1:length(hdgs)
        trs = ~isnan(data.choice) & data.heading==hdgs(h) & data.modality==mods(m) & data.delta==0;
        subplot(3,3,m); hold on; axis([0 5000 -40 40]); offsetAxes; title(sprintf('%s,T-x',modlabels{m}));
        plot(nexMat(:,trs,1),nexMat(:,trs,3),'color',hdgCol(h,:),'linewidth',0.5);
        subplot(3,3,m+3); hold on; axis([0 5000 0 160]); offsetAxes; title(sprintf('%s,T-y',modlabels{m}));
        plot(nexMat(:,trs,1),nexMat(:,trs,4),'color',hdgCol(h,:),'linewidth',0.5);
        subplot(3,3,m+6); hold on; axis([0 5000 -2 5]); offsetAxes; title(sprintf('%s,T-z',modlabels{m}));
        plot(nexMat(:,trs,1),nexMat(:,trs,5),'color',hdgCol(h,:),'linewidth',0.5);
    end
end

%% 
% get predicted trajectories, should be in MP cmds...

ampl = 160;
finalX = 160.*sind(hdgs);
finalY = 160.*cos(hdgs); 
finalZ = zeros(size(hdgs));





%% check timestamps vs delivery stamps

tr = randsample(length(nex.nexdata),50);
temp=squeeze(nexMat(:,tr,1:2));
figure('position',[200 200 300 500])
subplot(211)
plot(squeeze(temp(:,:,1)-temp(1,:,1)),squeeze(temp(:,:,2)-temp(1,:,1)),'r.','linew',2)
hold on
refline(1,1)
ylabel('Delivery Stamp')
subplot(212)
plot(squeeze(temp(:,:,1)-temp(1,:,1)),nexMat(:,tr,3)-nexMat(1,tr,3),'k.','linew',2)
xlabel('Timestamp')
ylabel('x')
ylim([-50 50])



%% choice-conditioned motion trajectories - subtract the mean trajectory for each condition


for m=[1 3]%1:length(mods)
    figure(101); subplot(3,1,m); hold on;
    for h=1:length(hdgs)
        trs = ~isnan(data.choice) & data.heading==hdgs(h) & data.modality==mods(m) & data.delta==0;
            
        choices = data.choice(trs);
        tempMat = nexMat(:,trs,:);
        [tempMat_mu,meanTraj] = demean(tempMat,2);
        
        L_trs = choices==1;
        R_trs = choices==2; 
        
        plot(mean(tempMat_mu(1:100,L_trs,3),2),mean(tempMat(1:100,L_trs,4),2),'r','linew',2)
        plot(mean(tempMat_mu(1:100,R_trs,3),2),mean(tempMat(1:100,R_trs,4),2),'b','linew',2)
    
        
    end
end








%% kluge cleanup some old data for Arielle 
% SJ 07-19-2021

lens= cellfun(@(x) size(x,1),nex.nexdata);

wrong_iTrials = find(nex.pldaps.trialSeed==7);
nex.pldaps.iTrial(wrong_iTrials) = nex.pldaps.iTrial(wrong_iTrials-1)+1;

split_trials = find(diff(nex.pldaps.iTrial)==0);

% trialSeedSplit(:,1) = nex.pldaps.trialSeed(split_trials);
% trialSeedSplit(:,2) = nex.pldaps.trialSeed(split_trials+1);

discardInds = split_trials+1;
%%%%
for i=1:length(split_trials)
    nex.nexdata{split_trials(i)} = [nex.nexdata{split_trials(i)}; nex.nexdata{split_trials(i)+1}];
    nex.nexdata(split_trials(i)+1) = [];
end

fnames = fieldnames(nex.conditions);
for f=1:length(fnames)
    nex.conditions.(fnames{f})(split_trials+1) = [];
end

fnames = fieldnames(nex.behavior);
for f=1:length(fnames)
    nex.behavior.(fnames{f})(split_trials+1) = [];
end


fnames = fieldnames(nex.pldaps);
for f=1:length(fnames)
    nex.pldaps.(fnames{f})(split_trials+1) = [];
end

%%
clearvars -except nex hdr


