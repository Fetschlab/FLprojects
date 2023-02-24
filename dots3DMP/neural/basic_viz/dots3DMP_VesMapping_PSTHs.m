% ves mapping


clear;clc;close all
addpath(genpath('/Users/stevenjerjian/Desktop/FetschLab/Analysis/codes/'))

% Load in the data

subject   = 'lucio';
dateRange = 20220615:20220923;

dataPath = '/Users/stevenjerjian/Desktop/FetschLab/Analysis/data/';
dataFileName = sprintf('%s_%d-%d_neuralData.mat',subject,dateRange(1),dateRange(end));
load(fullfile(dataPath,dataFileName));

%%

par     = 'VesMapping';

hdgTheta = [90 -90 0 180 0 0]';
hdgPhi   = [0 0 0 0 90 -90]';

condsVesMapping = [hdgTheta,hdgPhi];
condlabels = {'headingTheta','headingPhi'};

eventInfo.alignEvent  = {{'stimOn','stimOff'}}; % align to stimOn
eventInfo.tStart = -0.5; 
eventInfo.tEnd   = 0.5; % go to 0.5s after stimOff

% same for both
eventInfo.otherEvents = {{'stimOff'}};
eventInfo.otherEventnames = {{'stimOff'}};
eventInfo.binSize = 0.01;

optsTR.calcTuning  = 0;
optsTR.smoothFR    = 1;
optsTR.normalizeFR = 0;
optsTR.convKernel  = fspecial('average', [1 20]); 

allUnitsVesMapping = dots3DMP_FRmatrix_fromDataStruct(dataStruct,par,eventInfo,condsVesMapping,condlabels,optsTR);


%% plot PSTH of example unit


% 6 conditions only

u = 69;
iae=1;

tempFR = allUnitsVesMapping.data.PSTHs{iae}(:,:,u);
xt     = allUnitsVesMapping.times.xvec{iae};

figure('position',[100 100 600 400],'color','w'); box off;
hold on
plot(xt,tempFR','linew',2);
ylabel('Firing rate [spikes per sec]','fontsize',14)
changeAxesFontSize(gca,15,15);
set(gca,'tickdir','out')

for c=1:size(condsVesMapping,1)
    lhtxt{c} = sprintf('%s=%d, %s=%d',char(952),hdgTheta(c),char(966),hdgPhi(c));
end

my = max(tempFR(:))*1.25;
axis([xt([1 end]) 0 my]);
plot([0 0],[0 my],'k','linestyle','-');

for ioe = 1:size(allUnitsVesMapping.times.evTimes{iae},2)
    evMu = mean(allUnitsVesMapping.times.evTimes{iae}(:,ioe,u));

  
    plot([evMu evMu],[0 my],'k','linestyle','--');
    text(evMu*0.99,my*1,eventInfo.otherEventnames{iae}{ioe},'fontsize',12,'horizo','right','verti','bottom','rotation',90);
end

legend(lhtxt)

text(1,my*0.9,sprintf('%d-%d, %s',allUnitsVesMapping.hdr.unitDate(u),allUnitsVesMapping.hdr.unitID(u),dataStruct(end).data.(par).units.cluster_labels{allUnitsVesMapping.hdr.unitType(u)}),...
                'horizo','left','fontsize',12);