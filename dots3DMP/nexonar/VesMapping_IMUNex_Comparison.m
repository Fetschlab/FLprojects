
% 2022/03/15 SJ recorded simultaneous simple IMU and Nexonar for Vestibular
% Mapping (sinusoid) paradigm.

% Script processes this data and visualizes it.

tdatafolder = '/Users/stevenjerjian/Desktop/FetschLab/Analysis/data/20220315'; % trellis data folder
pdatafolder = '/Users/stevenjerjian/Desktop/FetschLab/Analysis/data/nexonar_test'; % pldaps and nexonar data

ns2file      = 'test20220315dots3DMP0002.ns2';
nsEventsfile = fullfile(tdatafolder,'test20220315dots3DMP0002_RippleEvents.mat');
PDSfile      = fullfile(pdatafolder,'test20220315VesMapping1250.mat');
NEXfile      = fullfile(pdatafolder,'test20220315VesMapping1250_nexonar.mat');


%% 

load(PDSfile);
load(NEXfile);

nex.behavior = rmfield(nex.behavior,{'choice','correct','RT','conf'});
[nexPDS,nexClean,exitflag] = dots3DMP_nexonarCleanUp(nex,PDS);

% freq/ampl and theta/phi come in pairs
freqampl = [nexClean.conditions.headingFreq; nexClean.conditions.headingAmpl]';
unique_freqampls = unique(freqampl,'rows');

thetaphi = [nexClean.conditions.headingTheta; nexClean.conditions.headingPhi]';
unique_thetaphis = unique(thetaphi,'rows');

nfa = size(unique_freqampls,1);
ntp = size(unique_thetaphis,1);

nexMat = nan(2000,5,length(nexPDS));
for t=1:length(nexPDS)
    nexMat(1:size(nexPDS{t},1),:,t) = nexPDS{t};
end
nexMat(:,1,:,:) = nexMat(:,1,:,:) - nexMat(1,1,:,:);


nexMean = nan(2000,5,nfa,ntp);
nexErr = nan(2000,5,nfa,ntp);

clear n
for fa=1:nfa
    for tp = 1:ntp
        
        J = nexClean.conditions.headingFreq==unique_freqampls(fa,1) & nexClean.conditions.headingAmpl==unique_freqampls(fa,2) & ...
            nexClean.conditions.headingTheta==unique_thetaphis(tp,1) & nexClean.conditions.headingPhi==unique_thetaphis(tp,2);
        n(fa,tp) = sum(J); % confirm 3 of each condition...
        
        nexMean(:,:,fa,tp) = nanmean(nexMat(:,:,J),3);
        nexErr(:,:,fa,tp)  = nanstd(nexMat(:,:,J),[],3)/sqrt(n(fa,tp));
        
    end
end

%%
Fs = 1024;
figure;
for t=1:ntp
    for f=1:nfa
    subplot(6,3,f+(t-1)*3)
    plot(nexMean(:,1,f,t)/Fs,nexMean(:,3:5,f,t)-nexMean(1,3:5,f,t));
    tt = sprintf('%d, %d \t %.2f, %d',unique_thetaphis(t,1),unique_thetaphis(t,2),unique_freqampls(f,1),unique_freqampls(f,2));
    title(tt)
    xlim([-1 18])
    %legend('x','y','z')
    end
end

% looks reasonable for the most part (in terms of x,y,z)
% minor displacement along non-chosen axes
% some issues at lowest freq?
% amplitudes aren't as much as they should be?

% what's nexonar's sample rate? % 1kHz?


%% IMU 

movavtime = 0.1;
% read in analog data
dchs = 10243:10245;

completeFilePath = fullfile(tdatafolder,ns2file);
for c=1:length(dchs)
    
    [nsData(c)] = getdata_NS(completeFilePath, 'Analog 1k', dchs(c));
end

% create data segments aligned to stimOn?
load(nsEventsfile);

stimOns = nsEvents.Events.stimOn; % in seconds
win     = [-1 18]; % seconds around alignment event
Fs = nsData(1).analogInfo.SampleRate;

samps = round((win(1) * Fs):((win(2) * Fs)));
times = samps/Fs;

stimOnWins = bsxfun(@plus,samps,round(stimOns * Fs)');

trialData = nan(size(stimOnWins,1),size(stimOnWins,2),3);
for c = 1:3
    trialData(:,:,c) = nsData(c).analogData(stimOnWins);
end
trialData = movmean(trialData,movavtime*Fs,2);
trialData = trialData - mean(trialData(:,1:20,:),2);

%

figure;
for tp=1:ntp
    for fa=1:nfa
    subplot(6,3,fa+(tp-1)*3)
    
        J = nexClean.conditions.headingFreq==unique_freqampls(fa,1) & nexClean.conditions.headingAmpl==unique_freqampls(fa,2) & ...
            nexClean.conditions.headingTheta==unique_thetaphis(tp,1) & nexClean.conditions.headingPhi==unique_thetaphis(tp,2);
        n(fa,tp) = sum(J); % confirm 3 of each condition...
        
        temp = squeeze(nanmean(trialData(J,:,:)));
%         temp = cumsum(temp,1);
        
        plot(times,temp);
        tt = sprintf('%d, %d \t %.2f, %d',unique_thetaphis(tp,1),unique_thetaphis(tp,2),unique_freqampls(fa,1),unique_freqampls(fa,2));
        title(tt)
        axis([times([1 end]) -10 10])
        %legend('x','y','z')

    end
end







