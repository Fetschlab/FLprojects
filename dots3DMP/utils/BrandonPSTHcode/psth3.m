%Code written by Brandon Nanfito

%% Load Data

clear;
jobID = 3;
numPresentations = 1;
channels = 19; %Channel before remapping of probe
experiment = 1; %Not really experiment file but for LIP that part of the code might just work fine
offlineSorting = 1;
nasOn = 0; %Are we looking through the NAS system
dateOfInterest = '2021-02-23_15-11-30_exp';
chanOfInterest = '_Ch19_3_multi';
fileOfInterest = [dateOfInterest chanOfInterest '.psort'];

if nasOn == 0
    filename = ['/Volumes/jk3/' dateOfInterest '/' fileOfInterest]; %2021-01-06_14-37-57_map\2021-01-06_14-37-57_map_Ch1_Sorting.psort'; %2020-12-01_15-18-07_map_Ch20_Sorting.psort';
    foldername = ['/Volumes/jk3/' dateOfInterest '/'];
else
    filename = ['Z:\fetschlab\data\hanzo_neuro\' dateOfInterest '\' fileOfInterest];
    foldername = ['Z:\fetschlab\data\hanzo_neuro\' dateOfInterest '\'];
end
pSorterResults = Psort_read_psort(filename); %Remember this is a long Continous file with Time Stamps
[openEvents, ~, openCont, ~] = neuralExtractorMapping(jobID, numPresentations, channels, experiment, offlineSorting, 0, 1, [], foldername);

firstSpikeTimes = openCont(pSorterResults.topLevel_data.ss_index == 1); %spike Times from offline sorter    
moreThanOne = unique(pSorterResults.topLevel_data.ss_index);
runNumber = 1;
if length(moreThanOne) > 1 %This might not be useful for when there is multiple depths (Would make things very complicated)
    secondSpikeTimes = openCont(pSorterResults.topLevel_data.cs_index == 1);
    runNumber = 2;
end

%% Settings
%"zs" is a matrix containing the zeros of each trial, where each row
%corresponds to a particular window, and each column is a trial. Set the
%width of the windows of interest by changing the values of "windowArms_s"
%(in seconds). 

%Define your windows of interest
windowArms_s(1,1) = .1; zs(1,:) = openEvents.motionStart; windowArms_s(1,2) = .5;
windowArms_s(2,1) = .5; zs(2,:) = openEvents.motionEnd  ; windowArms_s(2,2) = .5;
windowArms_s(3,1) = .5; zs(3,:) = openEvents.endAcquire ; windowArms_s(3,2) = .1;

windows(:,:,1) = [zs(1,:)-windowArms_s(1,1);zs(1,:)+windowArms_s(1,2)];
windows(:,:,2) = [zs(2,:)-windowArms_s(2,1);zs(2,:)+windowArms_s(2,2)];
windows(:,:,3) = [zs(3,:)-windowArms_s(3,1);zs(3,:)+windowArms_s(3,2)];

winLength_ms(1) = round((max(windows(2,:,1)-windows(1,:,1)))*1000);
winLength_ms(2) = round((max(windows(2,:,2)-windows(1,:,2)))*1000);
winLength_ms(3) = round((max(windows(2,:,3)-windows(1,:,3)))*1000);

bins{1} = (0:winLength_ms(1))-(windowArms_s(1,1)*1000);
bins{2} = (0:winLength_ms(2))-(windowArms_s(2,1)*1000);
bins{3} = (0:winLength_ms(3))-(windowArms_s(3,1)*1000);

condition = 4;
cond_colors = {'g' 'b' 'r' 'k' 'c' 'm'};
win_zeros = {'motion start' 'motion end' 'end acquire'};

%% Extract trials for conditions
%Here you will need to create a vector filled with the indices of all the
%trials of interest for each condition you want to plot, and assign them to
%a corresponding cell in the cell array "cond"

highBetInd = find(openEvents.pdw == 1);
numHigh = length(highBetInd);
lowBetInd = find(openEvents.pdw == 0);
numLow = length(lowBetInd);

direction = openEvents.direction == 180;
score = direction - openEvents.choice;
correctInd = find(score == 0 | score == -2);
incorrectInd = find(score == -1);

leftchoiceInd = find(openEvents.choice == 1);
rightchoiceInd = find(openEvents.choice == 2);

if condition == 1
    cond{1} = 1:length(openEvents.trialNum);
    
    cond_name = 'all'; 
end

if condition == 2
    co = unique(openEvents.coherence);
    for i = 1:length(co)

        cond{i} = find(openEvents.coherence == co(i));

    end
    
    cond_name = 'coh';
end

if condition == 3
    highCorInd = highBetInd(ismember(highBetInd,correctInd));
    cond{1} = highCorInd;
    lowCorInd = lowBetInd(ismember(lowBetInd,correctInd));%lowBetInd(find(ismember(lowBetInd,correctInd)));
    cond{2} = lowCorInd;
    highIncorInd = highBetInd(ismember(highBetInd,incorrectInd));
    cond{3} = highIncorInd;
    lowIncorInd = lowBetInd(ismember(lowBetInd,incorrectInd));
    cond{4} = lowIncorInd;
    
    cond_name = 'score';
end

if condition == 4
    highrightInd = highBetInd(ismember(highBetInd,rightchoiceInd));
    cond{1} = highrightInd;
    highleftInd = highBetInd(ismember(highBetInd,leftchoiceInd));
    cond{2} = highleftInd;
    lowrightInd = lowBetInd(ismember(lowBetInd,rightchoiceInd));
    cond{3} = lowrightInd;
    lowleftInd = lowBetInd(ismember(lowBetInd,leftchoiceInd));
    cond{4} = lowleftInd;
    
    cond_name = 'choice';
end

%% Extract event aligned spike times
%For each condition, loop through all the trials in that condition and
%extract the event alligned spiketimes, for each window. Store the event
%aligned spike times in "spikeData" such that each cell contains a cell
%array for every condition, and for each condition's cell array, the rows
%correspond to a particular window and the columns correspond to a 
%particular trial. 

for c = 1:length(cond)
    
    ti = cond{c};
    spikeData{c} = cell(length(windows(1,1,:)),length(ti));
    
    for w = 1:length(windows(1,1,:))
        
        for t = 1:length(ti)
            
            curT = ti(t);
            spikes = round((firstSpikeTimes(find(firstSpikeTimes>windows(1,curT,w),1,'first'):find(firstSpikeTimes<windows(2,curT,w),1,'last'),1)-zs(w,curT))*1000);
            spikeData{c}(w,t) = mat2cell(spikes,length(spikes),1);
            
        end
        
    end
    
end


%% Make the plots

%c = 1; 
%w = 3;

%These are just the x labels for the plots sdf plots, hard coded but I'm
%sure theres a better way to do it
binXLabel{1} = {'-100' '0' '100' '200' '300' '400' '500'};
binXLabel{2} = {'-500' '-400' '-300' '-200' '-100' '0' '100' '200' '300' '400' '500'};
binXLabel{3} = {'-500' '-400' '-300' '-200' '-100' '0' '100'};

%Looping through all three windows can be slow, commenting this loop out
%and specifying a window to look at may be more efficient for quickly
%looking at data
for w = 1:length(windows(1,1,:))
f = figure('Units','normalized','Position',[0 0 1 1]);
trial_count = 0;

    %Loop through each condition
    for c = 1:length(cond)

    %load the event aligned spike times for the trials of interest
    %(belonging to a specfic condition) into sptimes
    sptimes = spikeData{c}(w,:);

    %Make Raster Plot
    ax = subplot(4,1,1); hold on
    for iTrial = 1:length(sptimes) %loop through all trials

        spks = cell2mat(sptimes(iTrial)); %extract spike times
        xspikes = transpose(repmat(spks,1,3)); %assign spike times to x values
        yspikes = nan(size(xspikes));

        if ~isempty(yspikes)
            yspikes(1,:) = iTrial-1;
            yspikes(2,:) = iTrial;
        end

        yspikes = yspikes+trial_count; %this should make the y value increase by one for each trial

        plot(xspikes,yspikes,'Color',cond_colors{c},'LineWidth',1)

    end
    ax.YLabel.String = 'trial';
    ax.XLabel.String = ['time (ms) relative to ' win_zeros{w}];

    % PSTH

    all = [];
    for iTrial = 1:length(sptimes) %concatenate the spike times from all trials within a condition
        all = [all; sptimes{iTrial}];
    end

    ax = subplot(4,1,2); hold on
    nbins = length(bins{w});
    h = histogram(all,nbins); %make the histogram

    slength = length(bins{w});
    bdur = slength/nbins;
    nobins = 1000/bdur;
    h2 = histogram('BinEdges',h.BinEdges,'BinCounts',(h.Values/length(sptimes))*nobins); %convert the histogram into Hz
    delete(h);
    h2.FaceAlpha = .5;
    h2.FaceColor = cond_colors{c};

    mVal = max(h2.Values) + round(max(h2.Values)*.1);
    ax.XLim = [bins{w}(1) bins{w}(end)];
    ax.YLim = [0 mVal];
    ax.XLabel.String = ['time (ms) relative to ' win_zeros{w}];

    ax.YLabel.String = 'Firing rate (Hz)';
    if c == length(cond) %only need to add the legend once
        if condition == 1 %legend needs to change depending on condition
           legend('all')
        elseif condition == 2
           legend('c = 0.0000','c = 0.0320','c = 0.0640','c = 0.1280','c = 0.2560','c = 0.5120')
        elseif condition == 3
           legend('high bet; correct','low bet; correct','high bet; incorrect','low bet; incorrect')
        elseif condition == 4 
           legend('high right','high left','low right','low left')
        end
    end

    % Spike Density Function

    %tstep = 1;
    sigma = 5;
    time = bins{w};%tstep-100:tstep:500;

    for iTrial = 1:length(sptimes)

        spks = [];
        gauss = [];
        spks = sptimes{iTrial}';

        if isempty(spks)
            out = zeros(1,length(time));
        else

            for iSpk = 1:length(spks) %for each spike in this trial, create a gaussian curve that models it

                mu = spks(iSpk);
                p1 = -.5 *((time-mu)/sigma).^2;
                p2 = (sigma*sqrt(2*pi));
                gauss(iSpk,:) = (exp(p1)./p2)*1000; %this should be one gaussian curve centered at the current spike

            end

            sdf(iTrial,:) = sum(gauss,1); %sum all the curves to create a sdf for this trial, with each trial being a row in sdf (specific to current condition)
            sdfAllCond(iTrial+trial_count,:) = sum(gauss,1); %use this for the total sdf plot (subplot 3) since sdf is updated with each condition
        end
    end

    if c == length(cond) %only plot this when you have the sdf for all trials (once youve cycled through all the conditions)
        ax = subplot(4,1,3);
        imagesc(sdfAllCond)
        ax.XLabel.String = ['time (ms) relative to ' win_zeros{w}];
        ax.YLabel.String = 'Trials';
        ax.XLim = [1 length(bins{w})-1];
        ax.XTick = [0:100:length(bins{w})-1];
        ax.XTickLabel = binXLabel{w};
    end

    % Average Response
    ax = subplot(4,1,4); hold on
    sErrorMean(c,:) = SEM(sdf);
    signal = mean(sdf); timePoints = bins{w}+find(bins{w} == 0); %average the sdf curve for each trial in this condition
    yplus = signal+sErrorMean(c,:);
    yminus = signal-sErrorMean(c,:);
    fill([timePoints fliplr(timePoints)],[yplus fliplr(yminus)],cond_colors{c},'edgecolor','none'); %plot the SEM around the mean sdf
    alpha(0.2);
    plot(mean(sdf),'Color',cond_colors{c},'LineWidth',1.5)
    mVal_store(c) = max(mean(sdf)) + max(mean(sdf))*.1; %max(mean(sdf)) + round(max(mean(sdf))*.1);

    if c == length(cond)
        mVal = max(mVal_store);
        ax.YLim  = [0 mVal];
    end

    ax.XLim = [0 length(bins{w})-1];
    ax.XTick = [0:100:length(bins{w})-1];
    ax.XTickLabel = binXLabel{w};
    ax.YLabel.String = 'Firing rate (Hz)';
    ax.XLabel.String = ['time (ms) relative to ' win_zeros{w}];

    trial_count = trial_count + length(sptimes);
    clear sdf;clear spks;clear gauss; clear h; clear h2;

    end

    f.PaperUnits = 'inches';
    f.PaperSize = [22 17];
%     f.PaperPosition = [.05 .2 11.8 8];
    saveas(f,[dateOfInterest chanOfInterest '_' cond_name '_' num2str(w) '.pdf'])
    clear sdfAllCond; clear sErrorMean; 
    
end

% c2 = zeros(1,length(bins{1,3}));
% for b = 1:length(bins{1,3})
%     c = 0;
%     x = bins{1,3}(b);
% for con = 1:length(cond)
%     for i = 1:length(spikeData{1,con}(3,:))
%         if ismember(x,spikeData{1,con}{3,i})
%             c = c+1;
%         end
%     end
% end
%     c2(1,b) = c;
% end
% 
% f2 = figure;
% plot(c2)
% f2.PaperUnits = 'inches';
% f2.PaperSize = [22 17];
% %     f.PaperPosition = [.05 .2 11.8 8];
% saveas(f2,[dateOfInterest chanOfInterest '_spikecounts_' num2str(w) '.pdf'])
% clear sdfAllCond; clear sErrorMean; 




