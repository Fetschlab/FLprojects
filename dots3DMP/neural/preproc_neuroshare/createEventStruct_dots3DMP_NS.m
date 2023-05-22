% function nsEvents = createEventStruct_dots3DMP_NS(subject,date,paradigm,filenum,pldaps_filetime,filepath)
function nsEvents = createEventStruct_dots3DMP_NS(filename,par,pldaps_filename,filepath)

% SJ 05-2021, last updated 02-2022
% 
% extract the timing of key task events from .nev files
% timings are in 'Trellis' time, i.e. synced with recording of neural signals
% creates data structure, with trial info, outcomes, key event times,
% similar to data struct created by PLDAPS_preprocessing


% NOTES/TODO
% 07/2021 verify that timings here are consistent with 'actual' timings of events
% in PLDAPS/Datapixx i.e. no delays
    % 07-2021 manually checked a few numbers...seems to be <1-2ms difference 
    % between events decoded from Trellis recordings and PLDAPS data structs
    % 11-2021 checked timing again, differences ~1ms at most
% 07-2021 improved the checks for 'corrupt' data (see kluges below)
% 02-2022 added logic for dealing with dots3DMP_tuning paradigm - skips
            % some fields related to behavior in dots3DMP task paradigm
% 02-2022 since heading, coherence and delta are stored as the index in the
% condition list, rather than the raw value, we should pull these in from
% PLDAPS data structures at some point.

% packet drop issue - switch to TCP
% 06-2021 /done
% seems to resolve it, now only packet drops occur while configuration is
% being changed, and should just avoid this while recording anyway!

% PRE-REQUISITES: NeuroShare toolbox and getdata_NS function

% addpath(genpath('C:\Program Files (x86)\Ripple\Trellis\Tools'));
% addpath(genpath('C:\Users\fetschlab\Documents\MATLAB\dots3DMP')); 
                    
% e.g.
% subject  = 'lucio';
% date     = 20210517;      
% filenum  = 1;             % filename in Trellis (auto-increment numbers)
% pldaps_filetime = 1340;   % matching PLDAPS filenames, need to be
% manually entered as a list of matching length
% filepath = 'C:\Users\fetschlab\Trellis\dataFiles\';

% these are inputs directly now, parsing and re-formulating was unnecessary % 02/2022
% pldaps_filename = sprintf('%s%s%s%d',subject,date,paradigm,pldaps_filetime);
% filename = sprintf('%s%s%s%.04d.nev',subject,date,paradigm,filenum);
completeFilePath = fullfile(filepath,filename);

hdr.pldaps_filename = pldaps_filename;
hdr.Ripple_filename = filename;
hdr.Ripple_filepath = completeFilePath;

filedate = filename(6:13); % SJ 03/2022

dataType    = 'Digital';
dataChannel = 'parallel';  

% grab the raw NS data
[nsEvents] = getdata_NS(completeFilePath, dataType, dataChannel);

if strcmp(hdr.Ripple_filename,'lucio20220314dots3DMP0002.nev')
    % forgot to stop the Trellis recording again! doh
    nsEvents.data(1722:end) = [];
    nsEvents.time(1722:end) = [];
    nsEvents.sz(1722:end) = [];
elseif strcmp(hdr.Ripple_filename,'lucio20220401dots3DMP0001.nev')
    nsEvents.data(1204:end) = [];
    nsEvents.time(1204:end) = [];
    nsEvents.sz(1204:end) = [];
end

tr_start = strfind(nsEvents.data>=2^9,[1 1 1 1 1 1 1]); % find clocktime for trial delimiter
ntrs = length(tr_start);

% switch par
%     case 'VesMapping'
%         hdr.infolabels = {};
%     case 'RFMapping'
%         hdr.infolabels = {'coherenceInd','numDirs'};
%     case 'dots3DMP'
%         hdr.infolabels = {'headingInd','modality','coherenceInd','deltaInd','choice','correct','PDW'};
%     case 'dots3DMPtuning'
%         hdr.infolabels = {'headingInd','modality','coherenceInd','deltaInd'}; % no behavioral task here
% end

% allow for multiple paradigms in same file
hdr.infolabels = {'headingInd','modality','coherenceInd','deltaInd','choice','correct','PDW','numDirs'};

%%

% pre-allocate

% same for all paradigms
nsEvents.Events.trStart  = nan(size(tr_start));
% nsEvents.Events.trEnd    = nan(size(tr_start));
nsEvents.Events.fpOn     = nan(size(tr_start));
nsEvents.Events.fixation = nan(size(tr_start));
nsEvents.Events.reward   = nan(size(tr_start));

if strcmp(par(1:8),'dots3DMP')
    nsEvents.Events.stimOn  = nan(size(tr_start));
    nsEvents.Events.stimOff = nan(size(tr_start));
elseif strcmp(par,'RFMapping')
    nsEvents.Events.stimOn  = nan(6,length(tr_start));
    nsEvents.Events.stimOff = nan(6,length(tr_start));
end

if strcmp(par,'dots3DMP')
    nsEvents.Events.targsOn   = nan(size(tr_start));
    nsEvents.Events.saccOnset = nan(size(tr_start));
    nsEvents.Events.targHold  = nan(size(tr_start));
    nsEvents.Events.postTargHold = nan(size(tr_start));
end

for i=1:length(hdr.infolabels)
    nsEvents.Events.(hdr.infolabels{i}) = nan(size(tr_start));
end

%%
for t=1:ntrs
    
    if t==ntrs, t_en = length(nsEvents.data);
    else,       t_en = tr_start(t+1)-1;
    end
    
    % extract one trial's data
 
    tr_events = nsEvents.data(tr_start(t):t_en);
    tr_times  = nsEvents.time(tr_start(t):t_en);

    % tr_info contains stimulus information (hdgInd, mod, cohInd,
    % deltaInd), and behavior (choice, correct, PDW)
    tr_info = tr_events(tr_events>=2^8 & tr_events<2^9) - 2^8;

    if strcmp(par,'RFMapping') && ~isempty(tr_info)
        numDirs = tr_info(end);
        if length(tr_info)==numDirs+2
            nsEvents.Events.dotOrder(:,t) = tr_info(1:numDirs);
        else
            nsEvents.Events.dotOrder(:,t) = nan(numDirs,1);
        end
        tr_info(1:end-2) = [];
    end

    nsEvents.pldaps.iTrial(t) = tr_events(7)-(2^9+2^8);
    if str2double(filedate)==20210813 % messed up because we are ignoring bit 5 for the nexonar!, changing this back to only ignore when Dout is already high
        nsEvents.pldaps.iTrial(t)=t;
    end
    
    nsEvents.pldaps.unique_trial_number{t} = tr_events(1:6)-(2^9); % clocktime
    
    % when Nexonar is on, bit vals will be offset. account for this and
    % re-assign to tr_events
    %
    % NEXONAR ACQUISITION USES BIT 5 on the Grapevine Ripple end, correct as of 02/2022
    % ****changes to FLNexonar.cs must be reflected here!!!****
    % 
    % SJ 02-2022, RFmapping sends higher bits corresponding to dots
    % presentation order within the trial so this would be affected by
    % nexonar...ofc will probably never use nexonar in this paradigm, but would need to consider if we do something similar with VesMapping

    tr_eventsN = tr_events(tr_events<2^6);
    tr_eventsN = tr_eventsN - (tr_eventsN>=2^5)*2^5;
    tr_events(tr_events<2^6) = tr_eventsN;

    % a few checks...
    % try to be sure that each tr_events only contains ONE trial's events,
    % i.e. there should only be two 1 strobes for start/end of PLDAPS trial
    % PLDAPS trial numbers should also be incrementing by 1
    % number of info events should be consistent
   
    % 07-14-2021 this throws an error at 256 (2^8) trials because the lowWord gets screwed by the modulo...obviously. 
    % 07-16-2021 - see end of function for a posthoc fix for early tests, but should be ok from now on - changed strobe code to shift highWord to factor this in
    % if it's not fixed i.e. if pldaps.iTrial is not incrementing by 1 at any point, this will throw control to keyboard
    if ( t > 1 && (nsEvents.pldaps.iTrial(t)-nsEvents.pldaps.iTrial(t-1) ~= 1))  && str2double(filedate)>=20210716
        warning('PLDAPS iTrials not incrementing correctly, trial %d iTrial = %d, trial %d iTrial = %d',t-1,nsEvents.pldaps.iTrial(t-1),t,nsEvents.pldaps.iTrial(t));
%         keyboard
    end
    
    % packet drop issue was resolved by switching to TCP, see main notes
%     if sum(tr_events==1) ~= 2 || sum(tr_events==5) ~=2
%         warning('packet drop likely');
%     end
    
    % only reason a mismatch should occur is if trial was a breakfix and
    % behavioral outcomes were not recorded, in that case skip this trial
    % SJ 03-2022 if PLDAPS quit happens mid-trial and last trial is
    % aborted, could also happen, in which case skip
    if length(tr_info)~=length(hdr.infolabels)
        if t==ntrs % last trial
            disp('PLDAPS quit mid-trial on the last trial, skipping it')
            continue
        elseif sum(tr_events==10)==0 % brfix
            warning('Mismatch in trial info for trial %d but wasn''t a breakfix...need to check',t)
            continue
        end
    end
    
    %
    
    % store trialInfo
    for i=1:length(tr_info)
        nsEvents.Events.(hdr.infolabels{i})(t)  = tr_info(i);
    end
       
    % store EventTimes
    try nsEvents.Events.trStart(t)   	= tr_times(find(tr_events==1,1)); catch; end
%     if sum(tr_events==1)==2
%         try nsEvents.Events.trEnd(t)        = tr_times(find(tr_events==1,1,'last'));  catch; end
%     end
    
    try nsEvents.Events.fpOn(t)         = tr_times(tr_events==2);           catch; end
    try nsEvents.Events.fixation(t)     = tr_times(tr_events==3);           catch; end
    try nsEvents.Events.targsOn(t)      = tr_times(tr_events==4);           catch; end

    % STIMULUS ONSET/OFFSET
    if strcmp(par(1:8),'dots3DMP') || strcmp(par,'VesMapping')
        try nsEvents.Events.stimOn(t)       = tr_times(find(tr_events==5,1));   catch; end

        % if there are two motion events, the last one is stimOff (should be the same as saccOnset in RT task)
        if sum(tr_events==5)==2
            try nsEvents.Events.stimOff(t)      = tr_times(find(tr_events==5,1,'last'));  catch; end
        end

    elseif strcmp(par,'RFMapping') % should be numDirs*2 motion events
        try 
            stimOnOff = tr_times(tr_events==5); 
            if length(stimOnOff)==numDirs*2
                nsEvents.Events.stimOn(1:numDirs,t) = stimOnOff(1:2:end);
                nsEvents.Events.stimOff(1:numDirs,t) = stimOnOff(2:2:end);
            end 
        catch
            % should include a warning here? 
            % currently if numDirs doesn't match with the number of
            % stimOnOffs then they will all be left as NaN. this is probably fine,
            % the only way this happens is breakfix or expt ended
            % mid-trial, either way we don't care about these trials
        end
    end
    
    try nsEvents.Events.saccOnset(t)    = tr_times(tr_events==6);  catch; end % RT
    try nsEvents.Events.targHold(t)    	= tr_times(tr_events==7);  catch; end
    try nsEvents.Events.postTargHold(t) = tr_times(tr_events==8);  catch; end % PDW
    try nsEvents.Events.reward(t)       = tr_times(tr_events==9);  catch; end
    try nsEvents.Events.breakfix(t)     = tr_times(tr_events==10);  catch; end

    try nsEvents.Events.nexStart(t)       = tr_times(find(tr_events==12,1));   catch; end
    try nsEvents.Events.nexEnd(t)       = tr_times(find(tr_events==12,1,'last'));   catch; end

    % 01-2022 maybe should add in postTargsaccOnset (confRT)?
end

% store header info
nsEvents.hdr = hdr;

% 01/2022, clear these raw data fields, they are obsolete
nsEvents = rmfield(nsEvents,{'sz','data','time'});

% kluge for an old bug with trial number digIn
if ~all(diff(nsEvents.pldaps.iTrial)==1) && str2double(filedate)<20210716
   addFactor = repmat((0:ceil(ntrs/2^8)-1) * 2^8, 2^8-1, 1); % should be [0, 256, ...]
   addFactor = addFactor(:); 
   addFactor = addFactor(1:ntrs);
   nsEvents.pldaps.iTrial = nsEvents.pldaps.iTrial + addFactor;
end


