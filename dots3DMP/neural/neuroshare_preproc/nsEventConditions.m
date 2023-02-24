function [nsEvents] = nsEventConditions(nsEvents,PDS,par)
% [nsEvents] = NSEVENTCONDITIONS(nsEvents,PDS,PAR)
%
% pulls in useful trial info from relevant PLDAPS file into nsEvent struct
% 
% since events are sent to Ripple via binary DigIn, conditions except
% modality are sent as indices in condition list. Here we pull the actual
% info from PLDAPS, so that eventual neural dataStruct has all the task
% information needed for analysis
%
% SJ updates
% 05-2022 fixed issues with mismatched trials
% 06-2022 added 'goodtrial' from PDS, since this isn't in all Ripple files
% 08-2022 added import of 1-targ behavior from PLDAPS, also not sent to Trellis

%% check that unique trial numbers match completely, if not, select NS entries to match
% note that if multiple PDS files correspond to one Trellis file, this
% will return just the trials in nsEvents which correspond to the selected
% PDS file.
% this is important in the wrapper func, to make sure we are not
% overwriting nsEvents, or at least concatenating the output

for t=1:length(PDS.data)
    PDSutn(t,:) = PDS.data{t}.unique_number;
end
NSutn = nsEvents.pldaps.unique_trial_number;
if iscell(NSutn)
    NSutn = cell2mat(NSutn');
end

matchTrials = ismember(NSutn,PDSutn,'rows');
% fprintf('%d trials found in NSevents for this paradigm\n',sum(matchTrials))
% if sum(matchTrials)<length(matchTrials)
%     fprintf('%d trials did not match between PDS and NS, removing these...\n',sum(~matchTrials))
% end

% 06-2022
% breakfix vector was 1 entry short in one case, which messes up the
% matchTrials stuff. let's just remove it here, and just add in the goodtrial logical
% later, which serves the same purpose (this might be good to fix at some point)
try nsEvents.Events = rmfield(nsEvents.Events,'breakfix'); catch; end

% 02-2023 actually don't do this, in case we want to concatenate events from all sessions together

% remove fields only relevant to RFmapping if par is not RFmapping
% if ~strcmpi(par,'RFmapping')
%     try nsEvents.Events = rmfield(nsEvents.Events,'dotOrder'); catch; end
%     try nsEvents.Events = rmfield(nsEvents.Events,'stimOn_all'); catch; end
%     try nsEvents.Events = rmfield(nsEvents.Events,'stimOff_all'); catch; end
% end
% remove fields which are all nans, irrelevent for chosen paradigm
% (presumably)
% fnames = fieldnames(nsEvents.Events);
% for f = 1:length(fnames)
%     if ~iscell(nsEvents.Events.(fnames{f})) && all(isnan(nsEvents.Events.(fnames{f})))
%         try nsEvents.Events = rmfield(nsEvents.Events,fnames{f}); catch; end
%     end
% end

% seems like VesMapping events didn't have stimOff_all, so if VesMapping is
% done at the end the vector was the wrong length, so assign empties
% 02-22-2023
if nsEvents.pldaps.parNum(end)==4
    nsEvents.Events.stimOff_all(end+1:length(nsEvents.Events.stimOn_all)) = {[]};
end

% now replace all empties with stimOff if there's only one
empty_inds = find(cellfun(@isempty,nsEvents.Events.stimOff_all));
for c = 1:length(empty_inds)
    nsEvents.Events.stimOff_all(empty_inds(c)) = {nsEvents.Events.stimOff(empty_inds(c))};
end

%%
fnamesC = fieldnames(nsEvents.Events);
for f = 1:length(fnamesC)
    nsEvents.Events.(fnamesC{f}) = nsEvents.Events.(fnamesC{f})(matchTrials);
end
fnamesC = fieldnames(nsEvents.pldaps);
for f = 1:length(fnamesC)
    if strcmp(fnamesC{f},'unique_trial_number') && ~iscell(nsEvents.pldaps.unique_trial_number) % SJ 05-2022 some early files used a cell instead of mat, can re-create info files to fix this
        nsEvents.pldaps.unique_trial_number = nsEvents.pldaps.unique_trial_number(matchTrials,:);
    else
        nsEvents.pldaps.(fnamesC{f}) = nsEvents.pldaps.(fnamesC{f})(matchTrials);
    end
end

PDSdata = PDS.data;
PDSconditions = PDS.conditions;

NSutn=NSutn(matchTrials,:);

% in case PDS started before Trellis
matchTrials2 = ismember(PDSutn,NSutn,'rows');
PDSutn=PDSutn(matchTrials2,:);

if ~isequal(PDSutn,NSutn)
    disp('unique tr nums still do not match between PDS and NS...something more serious is wrong!\n'); keyboard
end

PDSconditions = PDSconditions(matchTrials2);
PDSdata = PDSdata(matchTrials2);

%% 



% refactor trial data from individual cells into matrix for easier insertion into nsEvents
% should be fixed across trials within a paradigm
fnamesC = fieldnames(PDS.conditions{1}.stimulus);

clear tempConds
for t = 1:length(PDSconditions)
    
    for f=1:length(fnamesC)
        tempConds.(fnamesC{f})(:,t) = PDSconditions{t}.stimulus.(fnamesC{f});
        try
            tempBehav.RT(:,t) = PDSdata{t}.behavior.RT;
        catch
            % fine
        end
    end
    
    % datapixx time is what is sent to Ripple
    PDStrstart(t) = PDSdata{t}.datapixx.unique_trial_time(1);
    
    nsEvents.Events.goodtrial(t) = PDSdata{t}.behavior.goodtrial;
    
    % SJ added 08-23-2022
    if strcmp(nsEvents.pldaps.parName{t},'dots3DMP')
        nsEvents.Events.oneTargChoice(t) = PDSdata{t}.behavior.oneTargChoice;
        nsEvents.Events.oneTargConf(t)   = PDSdata{t}.behavior.oneTargConf;
    else
        nsEvents.Events.oneTargChoice(t) = NaN;
        nsEvents.Events.oneTargConf(t) = NaN;
    end
end


% just add in the fields from PDS to nsEvents, without any further checks (except
% for the ones done above)
% TO DO
% maybe just do a sanity check on some, e.g. modality (since this is the
% same for both, or hdgOrder/dotOrder for RFmapping)

fnamesC = fieldnames(PDSconditions{1}.stimulus);
for f=1:length(fnamesC)
    if any(strcmp(fnamesC{f},{'headingTheta','hdgOrder','headingPhi'}))
        nsEvents.Events.(fnamesC{f}) = mat2cell(tempConds.(fnamesC{f})',ones(sum(matchTrials),1))';
    else
        nsEvents.Events.(fnamesC{f}) = tempConds.(fnamesC{f});
    end
end

try
    nsEvents.Events.RT = tempBehav.RT;
catch
end

% SJ 09-30-2022 should do the same with behavior - sanity check that chocie
% and PDW are transferred correctly, and also to take the RT from here!



%-----

% this is the old version, it has more checks and tests, but won't work for
% RFmapping, so needs adjustment..
%{
fnames = fieldnames(PDSconditions{1}.stimulus);
for f=1:length(fnames)
    
    % indices of the stimulus condition in nsEvents and temp should match,
    % but need to loop over blocks because list of indices may be different
    % for different blocks, even of the same paradigm
    
    % this doesn't work for RFmapping, oh and it assumes that all condition
    % are presented at least once!
    nsInds = nan(size(nsEvents.pldaps.blockNum));
    for b=1:length(nsEvents.hdr.pldaps_filenames)
        I = nsEvents.pldaps.blockNum==b;
        [~,~,nsInds(I)] = unique(temp.(fnames{f})(I));
    end
    
    if isfield(nsEvents.Events, ([fnames{f} 'Ind']))
        
        % this will also replace NaNs in Ripple Events with condition from PDS
        % (NaN in Ripple Events on brfix trials, but PDS condition is
        % specified from the beginning of the trial)
        F = nsEvents.Events.([fnames{f} 'Ind']);
        F(isnan(F)) = nsInds(isnan(F));
        
        if isequal(F,nsInds)
            nsEvents.Events.(fnames{f}) = temp.(fnames{f});
        else
            fprintf('%s indices do not match between PDS and NS\n',fnames{f})
            keyboard
        end
        
    elseif strcmp(fnames{f},'modality')
        
        F = nsEvents.Events.modality;
        nsEvents.Events.modality(isnan(F)) = temp.modality(isnan(F));
        
        if ~isequal(nsEvents.Events.modality,temp.modality)
            fprintf('modality does not match between PDS and NS\n')
            keyboard
        end
        
    end
end
%}


% PDStrstart = PDStrstart-PDStrstart(1)+nsEvents.Events.trStart(1);

% check relative timing of PDS trial start and Ripple trial start
% could use 'linearinterp' if there is any sense that shift is different
% across trials, but already seems super small (<1ms)

% fo = fit(PDStrstart',nsEvents.Events.trStart','poly1');
% fprintf('PDS/dpixx trStart and nsEvents trstart offset is ~%.3fms',fo.p2*1e3)

% would need to go through and shift ALL the timed events
%timingEvents = {'trStart',''};
% trStartsrefit = fo(nsEvents.Events.trStart');

    
   


    

