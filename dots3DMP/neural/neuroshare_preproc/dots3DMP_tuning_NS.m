% generate and plot tuning curves

% this is semi-online, i.e. run the short dots3DMP_tuning paradigm, then
% run this on saved data.

% run quick script to generate tuning curves and ISI plots for each channel
% (for any sorted or unsorted spikes)

% one subplot per channel, showing tuning curve for three modalities 

% SJ 05-2022 this will only work in current form if trellis recording only
% contains tuning paradigm (and in any case, would need to stop recording
% to run this semi-online, before starting task).

%% ENTER FILE INFO HERE

subject = 'lucio';
paradigm = 'dots3DMP';

date = 20220307;
trellis_filenum = '0001';


%%
filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(date) '\'];
filename = sprintf('%s%d%s%s_RippleEvents.mat',subject,date,paradigm,trellis_filenum);
load([filepath filename])


% THIS CODE SHOULD WORK FOR DOTS3DMP, or DOTS3DMP TUNING, just need to
% specify the appropriate headings used

% hdgvals = [-12 -6 -3 -1.5 0 1.5 3 6 12];
hdgvals = [-90 -45 -22.5 -12 0 12 22.5 45 90];
% hdgvals = [-45 -22.5 0 22.5 45];

showSorted = 1;

%% extract relevant variables for tuning curve plotting

% newer version 01/2022
motionStart = nsEvents.Events.stimOn;
motionEnd   = nsEvents.Events.stimOff;
hdg = nsEvents.Events.headingInd;
coh = nsEvents.Events.coherenceInd;
mod = nsEvents.Events.modality;
% motionStart = nsEvents.EventTimes.stimOn;
% motionEnd   = nsEvents.EventTimes.stimOff;
% hdg = nsEvents.trialInfo.headingInd;
% coh = nsEvents.trialInfo.coherenceInd;
% mod = nsEvents.trialInfo.modality;

hdgs = unique(hdg);
cohs = unique(coh);
mods = unique(mod);

ntrs = length(hdg);
nchs = length(nsEvents.spkData.chs);
% sorted_units = unique(nsEvents.spkData.unit_id); % nope, deal with this
% later
 
frRate      = nan(ntrs,nchs);
spkCount    = nan(ntrs,nchs);

spkCountSorted = [];
frRateSorted = [];
chanSorted = [];
unitSorted  = [];
sortedIndex = 0;

% prc2use = [25 75]; % center, width e,g, [25 75] use the middle 50% of the motion window 
prc2use = [0 100]; % use entire stim duration

for ch=1:nchs
    spktimes = nsEvents.spkData.data{ch}.spikeEventTime;
    sortedIndex = size(spkCountSorted,2);
    
    for t=1:ntrs
        
        if ~isnan(motionStart(t)) && ~isnan(motionEnd(t))
            
            % calculate middle segment of time
            temp = prctile(motionStart(t):1e-3:motionEnd(t),prc2use);
            
            spkCount(t,ch) = sum(spktimes > temp(1) & spktimes < temp(2));
            frRate(t,ch) = spkCount(t,ch) / (motionEnd(t)-motionStart(t));
            
            if showSorted
                sorted_ids = unique(nsEvents.spkData.data{ch}.unit_id);
                for u=sorted_ids
                   J = nsEvents.spkData.data{ch}.unit_id == u;
                   spkCountSorted(t,sortedIndex+1) =  sum(spktimes(J) > temp(1) & spktimes(J) < temp(2));
                   frRateSorted(t,sortedIndex+1)   =  spkCountSorted(t,end) / (motionEnd(t)-motionStart(t));
                   chanSorted(sortedIndex+1) = ch;
                   unitSorted(sortedIndex+1) = u;
                end
            end
        else
            if showSorted
                spkCountSorted(t,sortedIndex+1) =  nan;
                frRateSorted(t,sortedIndex+1)   =  nan;
                chanSorted(sortedIndex+1) = ch;
                unitSorted(sortedIndex+1) = u;
            end

        end
    end
%     disp(numSorted)
end

% spike count/firing rates averaged across trials for different conditions

frRate_byCond   = nan(length(mods),length(cohs),length(hdgs),nchs);
spkCount_byCond = nan(length(mods),length(cohs),length(hdgs),nchs);

frRateSorted_byCond   = nan(length(mods),length(cohs),length(hdgs),size(frRateSorted,2));
spkCountSorted_byCond = nan(length(mods),length(cohs),length(hdgs),size(spkCountSorted,2));

for m=1:length(mods)
    for c = 1:length(cohs)
        for h=1:length(hdgvals)
            J = mod==mods(m) & coh==cohs(c) & hdg==h; 
            N(m,c,h) = sum(J);
            
            if sum(J)==0
                continue
            end
            
            frRate_byCond(m,c,h,:) = nanmean(frRate(J,:),1);    
            spkCount_byCond(m,c,h,:) = nanmean(spkCount(J,:),1);
            
            frRateSorted_byCond(m,c,h,:) = nanmean(frRateSorted(J,:),1);
            spkCountSorted_byCond(m,c,h,:) = nanmean(spkCountSorted(J,:),1);
            
 
        end
    end
end



%% plotting

% one subplot per channel, fr vs hdg, each coh and mod plotted separately
% usual conventions

% black is ves
% magenta is low coh, red is high coh (vis)
% cyan is low coh, blue is high coh (comb)

% second plot is for online sorted units

if length(cohs)==2
    cols{1} = {'k','m','c'};
    cols{2} = {'k','r','b'};
elseif length(cohs)==1
    cols{1} = {'k','r','b'};
else
    error('not yet designed to handle >2 cohs')
end

% [p,q]=numSubplots(nchs);
p = [1,1];

figure('name','Unsorted MUA','color','w','position',[100 100 200*p(2) 200*p(1)]);
for ch=1:nchs
    subplot(p(1),p(2),ch); hold on

    for m=1:length(mods)
        for c=1:length(cohs)
            if mods(m)==1 && c>1, continue, end
            
            % why is spkCount producing modality-dependent results? it's
            % because of differences in motionEnd - motionStart
            % hmm, also some heading-specific ones, particularly for
            % combined (smaller headings have longer motion Time??
            % ah, could be because motionEnd in PLDAPS is checking
            % shownFrames AND currCmd...which could be different.
            
            % CHECK WHAT THE DEFINITION OF MOTION START AND END ARE FOR
            % EACH MODALITY (what is being sent to Ripple)
            
%             plot(hdgvals,squeeze(spkCount_byCond(m,c,:,ch)),'marker','.','markersize',10,'color',cols{c}{mods(m)});
            plot(hdgvals,squeeze(frRate_byCond(m,c,:,ch)),'marker','.','markersize',10,'color',cols{c}{mods(m)});
        end
    end
    xlim(hdgvals([1 end]));
    ylim([min(frRate_byCond(:)) max(frRate_byCond(:))].*[0.8 1.2])
    title(sprintf('ch %d',ch)) 
    set(gca,'xtick',hdgvals,'fontsize',8);
    xtickangle(45)
    if ch>1,set(gca,'xticklabels',''); end
end



if showSorted
    nSorted = size(frRateSorted,2);
    
    [p,q]=numSubplots(nSorted);
   
    figure('name','Sorted Units','color','w','position',[300 100 200*p(2) 200*p(1)]);
    for ch=1:nSorted
        subplot(p(1),p(2),ch); hold on
        
        for m=1:length(mods)
            for c=1:length(cohs)
                if mods(m)==1 && c>1, continue, end
                
                % why is spkCount producing modality-dependent results? it's
                % because of differences in motionEnd - motionStart
                % hmm, also some heading-specific ones, particularly for
                % combined (smaller headings have longer motion Time??
                % ah, could be because motionEnd in PLDAPS is checking
                % shownFrames AND currCmd...which could be different.
                
                % CHECK WHAT THE DEFINITION OF MOTION START AND END ARE FOR
                % EACH MODALITY (what is being sent to Ripple)
                
                %             plot(hdgvals,squeeze(spkCount_byCond(m,c,:,ch)),'marker','.','markersize',10,'color',cols{c}{mods(m)});
                plot(hdgvals,squeeze(frRateSorted_byCond(m,c,:,ch)),'marker','.','markersize',10,'color',cols{c}{mods(m)});
            end
        end
        xlim(hdgvals([1 end]));
        ylim([min(frRateSorted_byCond(:)) max(frRateSorted_byCond(:))].*[0.8 1.2])
        title(sprintf('ch %d, unit %d',chanSorted(ch),unitSorted(ch)))
        set(gca,'xtick',hdgvals,'fontsize',8);
        xtickangle(45)
        if ch>1,set(gca,'xticklabels',''); end
    end
end










