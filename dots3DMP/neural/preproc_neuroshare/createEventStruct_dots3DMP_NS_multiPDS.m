function nsEvents = createEventStruct_dots3DMP_NS_multiPDS(filename,paradigms,pldaps_filenames,filepath)
%
% spawned and adapted from createEventStruct_dots3DMP_NS, 04-2022
% 
% extract the timing of key task events from .nev files
% timings are in 'Trellis' time, i.e. synced with recording of neural signals
% creates data structure, with trial info, outcomes, key event times,
% this is now agnostic to paradigm and number of PDS files, so multiple PDS files can be recorded
% within the same trellis file!

% INPUTS:
% filename - trellis filename
% paradigms - list of paradigms associated with this filename
% pldaps_filenames - pldaps filenames matching paradigms above
% filepath 

% PRE-REQUISITES: NeuroShare toolbox and getdata_NS function

completeFilePath = fullfile(filepath,filename);

hdr.pldaps_filenames = pldaps_filenames;
hdr.Ripple_filename = filename;
hdr.Ripple_filepath = completeFilePath;

filedate = filename(6:13); % SJ 03/2022

dataType    = 'Digital';
dataChannel = 'parallel';  

% parNums sent from PLDAPS are 1-4 as follows
parOptions = {'dots3DMP','dots3DMPtuning','RFMapping','VesMapping'};

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
tr_end   = [tr_start(2:end)-1 length(nsEvents.data)];
ntrs     = length(tr_start);

% these can be easily outside the loop... 05-2022
% just need the indices relative to tr_start in strobe sequence, and
% subtract the appropriate highwords used for strobing
clockinds = bsxfun(@plus,tr_start',0:5);
unique_trial_number = nsEvents.data(clockinds)-2^9;
% unique_trial_number = mat2cell(unique_trial_number,ones(size(unique_trial_number,1),1));

iTrial = nsEvents.data(tr_start+6)-(2^9+2^8);

if str2double(filedate)>=20220504
    parNum = nsEvents.data(tr_start+8)-(2^9+2^8); % par number sent directly via strobe
else
    parNum = nan(1,ntrs);
end

% set paradigm number and trials to keep based on pldaps times 
% this way, discard any trials we don't want based on pldaps times in info
% (e.g. if we aborted a recording, or don't want to include a certain
% pldaps file, we will only keep trials pertaining to pldaps files listed
% in pldaps_filetimes

blockNum = nan(1,ntrs);
keeptrs = false(ntrs,1);

for p=1:length(pldaps_filenames)
    hhmm = pldaps_filenames{p}(end-3:end);
    % first trials should be at first hhmm after file start time, and when
    % pldaps iTrial is reset ([Inf iTrial], rather than iTrial==1, allows for pldaps to start
    % before Trellis recording, although this should be avoided if possible!)
    
    try
    firstTrial(p) = find( unique_trial_number(:,4)>=str2double(hhmm(1:2)) & unique_trial_number(:,5)>=str2double(hhmm(3:4)) & diff([Inf iTrial]')<=0,1);
%     firstTrial(p) = find( unique_trial_number(:,4)>=str2double(hhmm(1:2)) & unique_trial_number(:,5)>=str2double(hhmm(3:4)) & iTrial==1,1);
    catch
        disp('couldn''t find the relevant trials...check that pldaps_filenames are inputted correctly!')
        keyboard
    end
    try
        lastTrial(p) = firstTrial(p) - 1 + find(diff(iTrial(firstTrial(p):end))<0,1); 
    catch
        lastTrial(p) = ntrs;
    end
    keeptrs(firstTrial(p):lastTrial(p)) = 1;

    % set parNum based on par list, if it wasn't sent explicitly from PLDAPS
    if str2double(filedate)<=20220504% || ~strcmp(paradigms(p),'dots3DMP')
        parNum(firstTrial(p):lastTrial(p)) = find(strcmpi(parOptions,paradigms(p)));
    end
    blockNum(firstTrial(p):lastTrial(p)) = p;
end

tr_start = tr_start(keeptrs);
tr_end = tr_end(keeptrs);

ntrs = length(tr_start);
nsEvents.pldaps.unique_trial_number = unique_trial_number(keeptrs,:);
% nsEvents.pldaps.unique_trial_number = unique_trial_number(keeptrs); % ifcell
nsEvents.pldaps.iTrial = iTrial(keeptrs);
nsEvents.pldaps.parNum = parNum(keeptrs);
nsEvents.pldaps.blockNum = blockNum(keeptrs);

nsEvents.pldaps.parName = parOptions(nsEvents.pldaps.parNum);

% covers all paradigms
hdr.infolabels = {'headingInd','modality','coherenceInd','deltaInd','choice','correct','PDW','numDirs','nr','apertureDiam','targetR','targetTheta'};

%%

% pre-allocate

% same for all paradigms
nsEvents.Events.trStart  = nan(size(tr_start));
% nsEvents.Events.trEnd    = nan(size(tr_start));
nsEvents.Events.fpOn     = nan(size(tr_start));
nsEvents.Events.fixation = nan(size(tr_start));
nsEvents.Events.reward   = nan(size(tr_start));

nsEvents.Events.stimOn  = nan(size(tr_start));
nsEvents.Events.stimOff = nan(size(tr_start));
nsEvents.Events.motionOn = nan(size(tr_start));

nsEvents.Events.targsOn   = nan(size(tr_start));
nsEvents.Events.saccOnset = nan(size(tr_start));
nsEvents.Events.targHold     = nan(size(tr_start));
nsEvents.Events.postTargHold = nan(size(tr_start));
nsEvents.Events.reward    = nan(size(tr_start));
nsEvents.Events.breakfix  = nan(size(tr_start));
nsEvents.Events.nexStart  = nan(size(tr_start));
nsEvents.Events.nexEnd    = nan(size(tr_start));
nsEvents.Events.return    = nan(size(tr_start));


for i=1:length(hdr.infolabels)
    nsEvents.Events.(hdr.infolabels{i}) = nan(size(tr_start));
end

%%
for t=1:ntrs
    
    % extract one trial's data
    tr_events = nsEvents.data(tr_start(t):tr_end(t));
    tr_times  = nsEvents.time(tr_start(t):tr_end(t));

    % tr_info contains stimulus conditions and behavior outcomes
    tr_info = tr_events(tr_events>=2^8 & tr_events<2^9) - 2^8;

    % if we add 'nr' as a tr_info element this will need to change, based
    % on date of recording
    if strcmp(nsEvents.pldaps.parName{t},'RFMapping')
        
        % before ... wasn't sending dot coh and numDirs, just the
        % individual dot presentation order
        if str2double(filedate)<=20220301 
            numDirs = length(tr_info);
            try nsEvents.Events.dotOrder{t} = tr_info(1:numDirs); end
        else
            numDirs = tr_info(end);
%             nsEvents.Events.dotOrder{t} = nan(numDirs,1);


            % SJ revised, 08-05-2022, brfixes were causing an issue as not all
            % the dotOrder indices were being sent!
            % could fix this on the PLDAPS end too...
            if str2double(filedate)>=20220610
                numDotsShown = length(tr_info)-4;
            elseif length(tr_info)==2 
                numDotsShown = numDirs;
            else
                numDotsShown = length(tr_info)-2; % is this right? when did I add nr and apertureDiam
            end

%             if numDotsShown~=numDirs && ~any(tr_events==10)
%                 disp('numDotsShown does not match numDirs, but no brfix\n...')
%                 keyboard
%             end
            
            if length(tr_info)>2
                nsEvents.Events.dotOrder{t} = tr_info(1:numDotsShown);
                tr_info(1:numDotsShown) = [];
            end

%             if str2double(filedate)>=20220610
%                 if length(tr_info)==numDirs+4
%                     nsEvents.Events.dotOrder{t} = tr_info(1:numDirs);
%                     tr_info(1:end-4) = []; 
%                 elseif length(tr_info)==numDirs+2
%                     nsEvents.Events.dotOrder{t} = tr_info(1:numDirs);
%                     tr_info(1:end-2) = [];
%                 else 
%                     nsEvents.Events.dotOrder{t} = NaN;
%                     tr_info = tr_info(end-3:end);
%                 end
%             end
        end
    else
        nsEvents.Events.dotOrder{t}=NaN;
    end

    % assigned outside of the loop now SJ 04-2022
%     nsEvents.pldaps.unique_trial_number(t,:) = tr_events(1:6)-(2^9); % clocktime
%     nsEvents.pldaps.iTrial(t) = tr_events(7)-(2^9+2^8);

    if str2double(filedate)==20210813,nsEvents.pldaps.iTrial(t)=t; end % messed up because we are ignoring bit 5 for the nexonar!, changing this back to only ignore when Dout is already high 

    % when Nexonar is on, bit vals during trial will be offset by constant value. account for this and re-assign to tr_events
    %
    % NEXONAR ACQUISITION USES BIT 5 on the Grapevine Ripple end, correct as of 02/2022
    % ****changes to FLNexonar.cs must be reflected here!!!****
    % 
    % SJ 02-2022, RFmapping sends higher bits corresponding to dots presentation order within the trial so this would be affected by
    % nexonar...ofc will probably never use nexonar in this paradigm, but would need to consider if we do something similar with VesMapping

    nb = 5; % nexbit 5
    tr_eventsN = tr_events(tr_events<2^(nb+1));
    tr_eventsN = tr_eventsN - (tr_eventsN>=2^nb)*2^nb;
    tr_events(tr_events<2^(nb+1)) = tr_eventsN;

    % store trialInfo, paradigm-specific
    switch nsEvents.pldaps.parName{t}
        case 'dots3DMP'        
            infolabelInds = 1:7;
            if length(tr_info)==4 || ismember(10,tr_info)
                tr_info(end+1:end+3) = NaN; % put NaNs in place for brfix trials
            end
        case 'dots3DMPtuning',  infolabelInds = 1:4;
        case 'RFMapping'
%             if str2double(filedate)>=20220803
%                 infolabelInds = [3 11 12 10 9 8];
            if str2double(filedate)>=20220610
                infolabelInds = [3 10 9 8];
            elseif str2double(filedate)>20220301
                infolabelInds = [3 8];
            else 
                infolabelInds = [];
            end
        case 'VesMapping',      infolabelInds = [9];   
    end

    try
        for i=1:length(infolabelInds)
            try
                nsEvents.Events.(hdr.infolabels{infolabelInds(i)})(t)  = tr_info(i);
            catch
                nsEvents.Events.(hdr.infolabels{infolabelInds(i)})(t) = NaN;
            end
        end
    catch
        keyboard
    end

    % store EventTimes
    try nsEvents.Events.trStart(t)   	= tr_times(find(tr_events==1,1)); catch; end
%     if sum(tr_events==1)==2
%         try nsEvents.Events.trEnd(t)        = tr_times(find(tr_events==1,1,'last'));  catch; end
%     end
    if sum(tr_events==1)>2
        keyboard
        warning('more than 2 trial start/end events in trial %d...\n',t)
    end

    try nsEvents.Events.fpOn(t)         = tr_times(tr_events==2);           catch; end 
    try nsEvents.Events.fixation(t)     = tr_times(tr_events==3);           catch; end
    try nsEvents.Events.targsOn(t)      = tr_times(tr_events==4);           catch; end

    % STIMULUS ONSET/OFFSET
    % for RFMapping, this will be the very first stimOn and the final stimOff
    % stimOn_all and stimOff_all will include the individual starts and
    % ends of each dot direction presentation
    try nsEvents.Events.stimOn(t)   = tr_times(find(tr_events==5,1));   catch; end
    try nsEvents.Events.stimOff(t)  = tr_times(find(tr_events==5,1,'last'));  catch; end

    % SJ 2022-11-29
    % for dots3DMP RT task, RT is calculated from the time when motion
    % reaches 1/100th of peak acceleration (which occurs some time after
    % the 'stimulus'/motion period begins. So to match the RTs returned by
    % PLDAPS here, we want to use stimOff - 'motionOn'
    % motionOn is a new field added as of 0 
    if sum(tr_events==5)==3 %&& sum(tr_events==6)==1
        ev5 = find(tr_events==5);
        try nsEvents.Events.motionOn(t) = tr_times(ev5(2));   catch; end
    end

    switch nsEvents.pldaps.parName{t}
        case 'RFMapping' % should be numDirs*2 motion events
            try
                stimOnOff = tr_times(tr_events==5);
                nsEvents.Events.stimOn_all{t} = stimOnOff(1:2:end);
                nsEvents.Events.stimOff_all{t} = stimOnOff(2:2:end);
            catch
                fprintf('trial %d (RFMapping) did not have expected number of stim on-offs\n',t)
                % currently if numDirs doesn't match with the number of stimOnOffs, then they will all be left as NaN. this is probably fine,
                % the only way this happens is breakfix or expt was ended mid-trial, either way we would exclude these trials
            end

        case 'VesMapping'
            % ~start of each cycle SJ 2022-09-27
            % this should help to fold cycles later for IFR calcs
            nsEvents.Events.stimOn_all{t} = tr_times(tr_events==5);
        otherwise
            nsEvents.Events.stimOn_all{t} = nsEvents.Events.stimOn(t);
            nsEvents.Events.stimOff_all{t} = nsEvents.Events.stimOff(t);
    end
    
    try nsEvents.Events.saccOnset(t)    = tr_times(tr_events==6);  catch; end % RT
    try nsEvents.Events.targHold(t)    	= tr_times(tr_events==7);  catch; end
    try nsEvents.Events.postTargHold(t) = tr_times(tr_events==8);  catch; end % PDW
    try nsEvents.Events.reward(t)       = tr_times(tr_events==9);  catch; end
    try nsEvents.Events.breakfix(t)     = tr_times(tr_events==10);  catch; end

    if sum(tr_events==12)==2
        try nsEvents.Events.nexStart(t)     = tr_times(find(tr_events==12,1));   catch; end
        try nsEvents.Events.nexEnd(t)       = tr_times(find(tr_events==12,1,'last'));   catch; end
    end 

    try nsEvents.Events.return(t)     = tr_times(tr_events==13);  catch; end % 08-30-2022

end

% changed from reward, because reward time will be NaN on error trials, we want 'good trials' just to denote non-brfix (i.e. completed) trials
% in any case, this likely gets overwritten in the offline processing (see nsEventsConditions.m)
nsEvents.Events.goodtrial = isnan(nsEvents.Events.breakfix);

% store header info
nsEvents.hdr = hdr;

% 01/2022, clear these raw data fields, they are obsolete now
nsEvents = rmfield(nsEvents,{'sz','data','time'});
