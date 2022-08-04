%% create SessionData

% create a struct with neural Data for each session, see
% dots3DMP_NeuralPreProcessing for explanation of dataCell structure

% fields in 'events' which are time-based
tEvs   = {'trStart','fpOn','fixation','reward','stimOn','stimOff','saccOnset','targsOn','targHold','postTargHold','reward','breakfix'};

% dataCell = {};
dataCell = struct();

% sess = length(dataCell); % eventually allow for this code to append to existing dataCell if desired, instead of always starting from blank
sess = 0;

sflds = {'subject','date','pen','gridxy','probe_type','probe_ID'};
for n = 1:length(currentFolderList)
%     disp(currentFolderList{n})
    if isempty(strfind(currentFolderList{n},'20')) || contains(currentFolderList{n},'Impedance'); continue; end
    
%     try
    clear info
    load(fullfile(localDir,[subject currentFolderList{n} 'dots3DMP_info.mat']));
    fprintf('Adding data from %s, %d of %d\n',currentFolderList{n},n,length(currentFolderList))
    
    % we want 1 entry in dataCell for each unique 'recording session' -
    % each set
    
    [unique_sets,~,ic] = unique(info.rec_group);
    
    
    % SJ TO DO add some fprintf statements to report progress and what is
    % being added
    
    for u=1:length(unique_sets)
        
        % load neural data, if from kilosort (one file per set, either from the beginning, or after concatenation)
        if ~contains(info.probe_type{1},'Single')
            remoteDirSpikes = sprintf('/var/services/homes/fetschlab/data/%s/%s_neuro/%d/%s%d_%d/',subject,subject,info.date,subject,info.date,unique_sets(u));
            mountDir = sprintf('/Volumes/homes/fetschlab/data/%s/%s_neuro/%d/%s%d_%d/',subject,subject,info.date,subject,info.date,unique_sets(u));
            try 
                sp = loadKSdir(mountDir);
                [spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] = templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);
                [~,max_site] = max(max(abs(sp.temps),[],2),[],3);
                chs = sp.ycoords(max_site);
                depths = templateDepths;
            catch
                fprintf('Could not load kilosort sp struct for %s, set %d',info.date,unique_sets(u));
                continue
            end
        end
        
        sess = sess+1;
        
        % assign fields relevant to given recording 'set'/session
%         for f=1:length(sflds)
%             dataCell{sess}.info.(sflds{f}) = info.(sflds{f});
%         end
%         dataCell{sess}.info.set = unique_sets(u);
        
        dataCell(sess).date = info.date;
        dataCell(sess).info = info;
        dataCell(sess).set = unique_sets(u);
        
        % loop over paradigms (this is slightly different to how the
        % nsEvents are created!)
        for par=1:length(paradigms)
            
            theseFiles =  find((ic'==unique_sets(u)) & ismember(lower(info.par),lower(paradigms{par})) & (~isnan(info.pldaps_filetimes)));
            
            if isempty(theseFiles), continue, end
            
            %fprintf('%d files for %s\n',length(theseFiles),paradigms{par});
            
            % concatenate all the data and condition fields of PDS files for given paradigm which are marked as part of the same recording group
            clear allPDS
            st = 1;
            for pf=1:length(theseFiles)
                clear PDS
                if strcmp(info.par{theseFiles(pf)},'RFMapping')
                    info.par{theseFiles(pf)} = 'RFmapping';
                end
                PDSfilenames{pf} =  [info.subject num2str(info.date) info.par{theseFiles(pf)} num2str(info.pldaps_filetimes(theseFiles(pf))) '.mat'];
                try load(fullfile(PDSdir,PDSfilenames{pf}),'PDS');
                catch, fprintf('PDS file not found %s\n',PDSfilenames{pf}); return;
                end
                en = st-1+length(PDS.data);
                allPDS.data(st:en)       = PDS.data;
                allPDS.conditions(st:en) = PDS.conditions;
                st = en+1;
            end
            
            % for each trellis file within a given set+paradigm,
            % concatenate and store the events, and compute the necessary
            % shift for the spiketimes (if there were multiple trellis
            % files concatenated for kilosort, we will need to subtract the
            % shifted timestamps to line the spikes up with events again)
            [unique_trellis_files,~,ii] = unique(info.trellis_filenums(theseFiles));
            
            currPos = 0;
            
            if ~contains(info.probe_type{1},'Single')
                thisParSpikes  = false(size(sp.st));
                shiftSpikeTime = zeros(size(sp.st));
            else
                sp.st = [];
                sp.clu = [];
                sp.waves = [];
                shiftSpikeTime = [];
            end
            
            % now loop over each trellis file within a particular paradigm
            for utf=1:length(unique_trellis_files)
                NSfilename  = sprintf('%s%ddots3DMP%04d_RippleEvents.mat',info.subject,info.date,unique_trellis_files(utf));
                load(fullfile(localDir,NSfilename));
                
                % pull in relevant condition data from PLDAPS and sub-select trials from this paradigm
                
                % I think I messed up here and (permanently) deleted the wrong RFmapping
                % PLDAPS file (should have kept 1344 not 1338 as it seems
                % that is the one matched to 0008 trellis file)
                if strcmp(NSfilename,'lucio20220719dots3DMP0008_RippleEvents.mat')
                    
                    continue
                    
                    % don't have targetR and targetTheta though...adding
                    % this to Ripple events from now on
%                     thisParEvents = nsEvents;
%                     hdgThetas = [90 -90 0 180 0 0]; hdgPhis = [0 0 0 0 90 -90];
%                     for t = 1:length(thisParEvents.Events.dotOrder)
%                         if all(~isnan(thisParEvents.Events.dotOrder{t}))
%                             thisParEvents.Events.headingTheta{t} = hdgThetas(thisParEvents.Events.dotOrder{t});
%                             thisParEvents.Events.headingPhi{t} = hdgPhis(thisParEvents.Events.dotOrder{t});
%                         else
%                             thisParEvents.Events.headingTheta{t} = NaN;
%                             thisParEvents.Events.headingPhi{t} = NaN;
%                         end
%                     end
%                     thisParEvents.Events.hdgOrder = thisParEvents.Events.dotOrder;
                else
                    [thisParEvents] = nsEventConditions(nsEvents,allPDS,lower(paradigms{par}));
                end
                
                timeStampsShifted = thisParEvents.analogInfo.timeStampsShifted ./ double(thisParEvents.analogInfo.Fs);

                % do some concatenation in pldaps and events fields in case
                % the same par+block is split over multiple files (unlikely)
                nTr    = length(thisParEvents.Events.trStart);
                fnames = fieldnames(thisParEvents.Events);
                for f=1:length(fnames)
                    
                    if ismember(fnames{f},tEvs)
                        %                     dataCell{sess}.data.(paradigms{par}).events.(fnames{f})(currPos+1:currPos+nTr) = thisParEvents.Events.(fnames{f}) + timeStampsShifted(1);
                        dataCell(sess).data.(paradigms{par}).events.(fnames{f})(currPos+1:currPos+nTr) = thisParEvents.Events.(fnames{f})  + timeStampsShifted(1);
                    else
                        %                     dataCell{sess}.data.(paradigms{par}).events.(fnames{f})(currPos+1:currPos+nTr) = thisParEvents.Events.(fnames{f});
                        dataCell(sess).data.(paradigms{par}).events.(fnames{f})(currPos+1:currPos+nTr) = thisParEvents.Events.(fnames{f});
                    end
                    
                end
                fnames = fieldnames(thisParEvents.pldaps);
                for f=1:length(fnames)
                    if strcmp(fnames{f},'unique_trial_number') && ~iscell(nsEvents.pldaps.unique_trial_number)
%                         dataCell{sess}.data.(paradigms{par}).pldaps.(fnames{f})(currPos+1:currPos+nTr,:) = thisParEvents.pldaps.(fnames{f});
                        dataCell(sess).data.(paradigms{par}).pldaps.unique_trial_number(currPos+1:currPos+nTr,:) = thisParEvents.pldaps.(fnames{f});
                    else
%                         dataCell{sess}.data.(paradigms{par}).pldaps.(fnames{f})(currPos+1:currPos+nTr) = thisParEvents.pldaps.(fnames{f});
                        dataCell(sess).data.(paradigms{par}).pldaps.(fnames{f})(currPos+1:currPos+nTr) = thisParEvents.pldaps.(fnames{f});
                    end
                end
                dataCell(sess).data.(paradigms{par}).pldaps.blockNum2(currPos+1:currPos+nTr,:) = utf;
                currPos = currPos+nTr;

                % deal with single electrode recording neural data
              
                if contains(info.probe_type{1},'Single')
                    remoteDirSpikes = sprintf('/var/services/homes/fetschlab/data/%s/%s_neuro/%d/%s%ddots3DMP%04d/',subject,subject,info.date,subject,info.date,unique_trellis_files(utf));
                    mountDir = sprintf('/Volumes/homes/fetschlab/data/%s/%s_neuro/%d/%s%ddots3DMP%04d/',subject,subject,info.date,subject,info.date,unique_trellis_files(utf));
                    
                    cmd = ['ssh ' IPadd ' ls ' remoteDirSpikes];
                    [~,remoteFileList] = system(cmd);
                    
                    newlines = strfind(remoteFileList,newline);
                    remoteFiles = cell(length(newlines),1);
                    for nl = 1:length(newlines)
                        if nl==1, remoteFiles{nl} = remoteFileList(1:newlines(1)-1);
                        else, remoteFiles{nl} = remoteFileList(newlines(nl-1)+1:newlines(nl)-1);
                        end
                    end
                    
                    try
                        load([mountDir remoteFiles{contains(remoteFiles,'waveforms')}]);
                    catch
                        if unique_trellis_files(utf)==6, keyboard, end
                        fprintf('Could not find or load waveforms file %d, id%04d\n',info.date,unique_trellis_files(utf))
                        continue
                    end
                        
                    % rename/refactor some vars to match kilosort style
                    % mksort allows up to 4 units per ch, treat these
                    % clusters as the 'cluster ids'
                    
                    % concatenate spike times across recordings of the same
                    % paradigm within a set (we will shift the times below)
                    sp.st   = [sp.st; (waveforms.spikeTimes')/1000]; % mksort waveform times are in ms it seems
                    sp.clu  = [sp.clu; waveforms.units'];
                    
                    sp.waves = [sp.waves waveforms.waves];
                    
                    % for mksort, assume that clustering is consistent
                    % throughout!
                    if utf==1 
                        sp.cids = unique(sp.clu(sp.clu>0));
                        
                        % manual ratings as SU/MU/noise (4 or 3 --> 2 for SU, 2 or 1 --> 1 for MU, 0 as noise clusters --> 3)
                        sp.cgs  = ceil(waveforms.ratings.ratings(sp.cids) / 2);
                        sp.cgs(sp.cgs==0) = 3;
                    end
                    thisParSpikes  = true(size(sp.st));
                    shiftSpikeTime = [shiftSpikeTime; timeStampsShifted(1)*ones(size(waveforms.spikeTimes'))];
        
                else
                    % pick out spikes from the sp.st vector which can be linked
                    % to this paradigm's timeframe (with a good buffer on either
                    % side, e.g. 20secs), and what time to shift the spike
                    % times by (if any) so that they align with events again
                    % this shift is only necessary if multiple Trellis
                    % recordings were made for same location - these will
                    % have been concatenated for kilosort sorting, but
                    % events will still be separate
                    
                    timeLims = timeStampsShifted(1) + thisParEvents.Events.trStart([1 end]) + [-1 1]*20;
                    
                    thisFileSpikes    = (sp.st >= timeLims(1) & sp.st < timeLims(2));
                    thisParSpikes     = thisParSpikes | thisFileSpikes;
                    
                    % set the shifted time for each spike associated with a particular
                    % file to reconcile with nsEvents. NOTE: This is
                    % file- rather than paradigm specific, as there may be
                    % more than one file for a given paradigm
                    % also note that this will miss spikes outside the
                    % timeLims above, but we don't really care about them
                    % anyway because we select for thisParSpikes below
                    shiftSpikeTime(thisFileSpikes) = timeStampsShifted(1);
                    
                end
            end
            
            if isempty(sp.st), continue, end
            
            % shift the spike times now
            % add shifted time because we have shifted the nsEvents too
            % when storing them above.
            spikeTimes = sp.st + shiftSpikeTime;
          
            if keepMU, inds = sp.cgs<3;
            else,      inds = sp.cgs==2;
            end
            
            cids = sp.cids(inds);
            cgs  = sp.cgs(inds); 
            
%             dataCell{sess}.data.(paradigms{par}).cluster_id = cids;
%             dataCell{sess}.data.(paradigms{par}).cluster_type = cgs;
            
            dataCell(sess).data.(paradigms{par}).cluster_id = cids;
            dataCell(sess).data.(paradigms{par}).cluster_type = cgs;

            % SJ 06/13/2022 
            % store channel and depth information for each cluster (should do this if we
            % want to reconstruct putative recording locations)
%             dataCell(sess).data.(paradigms{par}).chs = chs;
%             dataCell(sess).data.(paradigms{par}).depths = depths;

            
            % now, assign spiketimes for each unit to given cell
            for unit=1:sum(inds)
                
                theseSpikes = sp.clu==cids(unit) & thisParSpikes;   
%                 dataCell{sess}.data.(paradigms{par}).spiketimes{unit} = spikeTimes(theseSpikes);
                dataCell(sess).data.(paradigms{par}).spiketimes{unit} = spikeTimes(theseSpikes);
                dataCell(sess).data.(paradigms{par}).waveforms{unit} = sp.waves(:,theseSpikes);
            end
            
        end
    end
%     catch
%         fprintf('Something not working, possibly file not found for %s...',currentFolderList{n})
%         continue
%     end
        
end

file = [subject '_' num2str(dateRange(1)) '-' num2str(dateRange(end)) '_neuralData.mat'];
disp('saving...');
save([localDir(1:length(localDir)-length(subject)-7) file], 'dataCell');