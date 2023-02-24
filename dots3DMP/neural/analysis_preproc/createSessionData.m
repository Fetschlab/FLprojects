%% create SessionData

% create a struct with neural Data for each session, see
% dots3DMP_NeuralPreProcessing for explanation of dataStruct structure

% IMPORTANT: 
% For spike sorting purposes, successive recordings with the same putative
% units (i.e. at the same location) are concatenated.
% The timestamps for concatenated recordings are shifted according to the
% length of the overall data so that the range of events and spikes is
% matched for a given recording (nsEvents.analogData.timeStamps and .timeStampsShifted).
% i.e. if a recording is the first in the set, it's timestamps should start from around 0, whereas 
% % if it is later in the set they will start from some other time >0. 
% The spiketimes will be matched accordingly, so the range of spiketimes should approximately match the range of (and be aligned to) task events. 


% fields in 'events' which contain event times, this will be important later
tEvs   = {'trStart','fpOn','fixation','reward','stimOn','stimOff','saccOnset',...
    'targsOn','targHold','postTargHold','reward','breakfix','nexStart','nexEnd','return'};

% sess = length(dataStruct); % eventually allow for this code to append to
% existing dataStruct if desired, instead of always starting from blank?

dataStruct = struct();
sess = 0;

sflds = {'subject','date','pen','gridxy','probe_type','probe_ID'};
for n = 1:length(currentFolderList)
%     disp(currentFolderList{n})
    if isempty(strfind(currentFolderList{n},'20')) || contains(currentFolderList{n},'Impedance'); continue; end
    
    clear info
    load(fullfile(localDir,[subject currentFolderList{n} 'dots3DMP_info.mat']));
    fprintf('Adding data from %s, %d of %d\n',currentFolderList{n},n,length(currentFolderList))
    
    % we want 1 row in dataStruct for each unique 'recording set'
    [unique_sets,~,ic] = unique(info.rec_group);
    
    for u=1:length(unique_sets)

        sess = sess+1; % increment the row in dataStruct

        remoteDirSpikes = sprintf('/var/services/homes/fetschlab/data/%s/%s_neuro/%d/%s%d_%d/',subject,subject,info.date,subject,info.date,unique_sets(u));
        mountDir = sprintf('/Volumes/homes/fetschlab/data/%s/%s_neuro/%d/%s%d_%d/',subject,subject,info.date,subject,info.date,unique_sets(u));

        if contains(info.probe_type{1},'Single')
            mountDir = [mountDir 'phy_WC/']; end

        try
            disp(mountDir)
            sp = loadKSdir(mountDir);
        catch
            warning('dots3DMP:createSessionData:loadKSdir','Could not load kilosort sp struct for %d, set %d\n\n',info.date,unique_sets(u));
        end
        try
            unitInfo = getUnitInfo(mountDir, keepMU);
        catch
            warning('dots3DMP:createSessionData:getUnitInfo','Could not get cluster info for this ks file..file has probably not been manually curated\n')
        end

        dataStruct(sess).date = info.date;
        dataStruct(sess).info = info;
        dataStruct(sess).set = unique_sets(u);
        
        % loop over paradigms 
        % (NOTE that the logic here differs from how nsEvents is initially created on the experiment rig)
        for par=1:length(paradigms)
            
            theseFiles =  find((ic'==unique_sets(u)) & ismember(lower(info.par),lower(paradigms{par})) & (~isnan(info.pldaps_filetimes)));
            
            if isempty(theseFiles), continue, end
                        
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
                catch, fprintf('PDS file not found..are you connected to the NAS?%s\n',PDSfilenames{pf}); 
                    keyboard
                    return;
                end
                en = st-1+length(PDS.data);
                allPDS.data(st:en)       = PDS.data;
                allPDS.conditions(st:en) = PDS.conditions;
                st = en+1;
            end
            
            % for each trellis file within a given set+paradigm, concatenate and store the events, and compute the necessary
            % shift for the spiketimes as well

            [unique_trellis_files,~,ii] = unique(info.trellis_filenums(theseFiles));
            
            currPos = 0;
            
            %if ~contains(info.probe_type{1},'Single')
                thisParSpikes  = false(size(sp.st));
%                 shiftSpikeTime = zeros(size(sp.st));
            %else
            %    sp.st = [];
            %    sp.clu = [];
            %    shiftSpikeTime = [];
            %end
            
            % now loop over each trellis file within a particular paradigm
            for utf=1:length(unique_trellis_files)
                NSfilename  = sprintf('%s%ddots3DMP%04d_RippleEvents.mat',info.subject,info.date,unique_trellis_files(utf));
                
                % messed up with PDS files on this one, oops
                if strcmp(NSfilename, 'lucio20220719dots3DMP0008_RippleEvents.mat'), continue, end

                try
                    load(fullfile(localDir,NSfilename));
                catch
                    fprintf('Could not load %s, skipping...\n\n', NSfilename)
                    continue
                end
                
                fprintf('adding data from %s (%s)\n\n', NSfilename, paradigms{par})

   
                % pull in relevant condition data from PLDAPS and sub-select trials from this paradigm

                [thisParEvents]   = nsEventConditions(nsEvents,allPDS,lower(paradigms{par})); % % SJ added 08-22-2022 oneTargChoice and Conf!
                timeStampsShifted = thisParEvents.analogInfo.timeStampsShifted ./ double(thisParEvents.analogInfo.Fs);

                % do some concatenation in pldaps and events fields, in case the same par+block is split over multiple trellis files
                % NOTE: currently two (or more) runs of the same paradigm within a block will be pooled.
                % what about if we record tuning at the beg and end at the same location, but want to look at them separately? 
                % two options:
                % 1. mark as different sets from the beginning, spike sort independently
                % 2. split post-hoc based on times of spikes, for visualization/comparison

                nTr    = length(thisParEvents.Events.trStart);

                % add all the fields in events and pldaps to dataStruct
                % if event field is a time, shift it as needed
                fnames = fieldnames(thisParEvents.Events);         
                for f=1:length(fnames)
                    
                    if ismember(fnames{f},tEvs)
                        thisParEvents.(fnames{f}) = thisParEvents.Events.(fnames{f})  + timeStampsShifted(1);
                    end
                    dataStruct(sess).data.(paradigms{par}).events.(fnames{f})(1,currPos+1:currPos+nTr) = thisParEvents.Events.(fnames{f});
                end


                fnames = fieldnames(thisParEvents.pldaps);
                for f=1:length(fnames)
                    if strcmp(fnames{f},'unique_trial_number') && ~iscell(nsEvents.pldaps.unique_trial_number)
                        dataStruct(sess).data.(paradigms{par}).pldaps.unique_trial_number(currPos+1:currPos+nTr,:) = thisParEvents.pldaps.(fnames{f});
                    else
                        dataStruct(sess).data.(paradigms{par}).pldaps.(fnames{f})(1,currPos+1:currPos+nTr) = thisParEvents.pldaps.(fnames{f});
                    end
                end

                dataStruct(sess).data.(paradigms{par}).pldaps.blockNum2(1,currPos+1:currPos+nTr) = utf;
                currPos = currPos+nTr;

                % deal with single electrode recording neural data
                % 01-2023 this is no longer necessary - now I am sorting
                % single elec recordings with wave_clus via SpikeInterface,
                % and then exporting to Phy, so format is the same

                %{
                    %                 if contains(info.probe_type{1},'Single')
                    remoteDirSpikes = sprintf('/var/services/homes/fetschlab/data/%s/%s_neuro/%d/%s%ddots3DMP%04d/',subject,subject,info.date,subject,info.date,info.trellis_filenums(utf));
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
                        fprintf('Could not find or load waveforms file %d, id%04d\n',info.date,unique_trellis_files(utf))
                        continue
                    end

                    % rename/refactor some vars to match kilosort style
                    % mksort allows up to 4 units per ch - treat these
                    % clusters as the 'cluster ids', assuming that we've
                    % been consistent in labelling across separate files

                    % concatenate spike times across recordings of the same
                    % paradigm within a set (we will shift the times below)
                    sp.st   = [sp.st; (waveforms.spikeTimes')/1000]; % mksort waveform times are in ms it seems
                    sp.clu  = [sp.clu; waveforms.units'];

                    if utf==1
                        sp.cids = unique(sp.clu(sp.clu>0));

                        % manual ratings as SU/MU/unsorted (4 or 3 --> 2 for SU, 2 or 1 --> 1 for MU, 0 as noise clusters --> 3)
                        sp.cgs  = ceil(waveforms.ratings.ratings(sp.cids) / 2);
                        sp.cgs(sp.cgs==0) = 3;
                    end
                    thisParSpikes  = true(size(sp.st));
                    shiftSpikeTime = [shiftSpikeTime; timeStampsShifted(1)*ones(size(waveforms.spikeTimes'))];

                    %                 else % kilosort data
                end            
                %}

                % pick out spikes from the sp.st vector which can be linked to this paradigm's timeframe (with a reasonable buffer on either
                % side, e.g. 20secs), and what time to shift the spike times by (if any) so that they align with events again
                % Note. this shift is only necessary if multiple Trellis recordings were made for same location - these will
                % have been concatenated for kilosort sorting, but events will still be separate      

                timeLims = timeStampsShifted(1) + thisParEvents.Events.trStart([1 end]) + [-1 1]*20;

                thisFileSpikes = (sp.st >= timeLims(1) & sp.st < timeLims(2));
                thisParSpikes  = thisParSpikes | thisFileSpikes; % union

                % set the shifted time for each spike associated with a particular file, to reconcile with events.
                % NOTE: This is file- rather than paradigm-specific, as there may be more than one file for a given paradigm.
                % also note that this will miss spikes outside the timeLims above, but we don't really care about them
                % anyway because we select for thisParSpikes below
%                 shiftSpikeTime(thisFileSpikes) = timeStampsShifted(1);

                % SJ 01-2023 this was only necessary with mkSort, where I was sorting each trellis file separately, but
                % wanted to merge events and spikes for one paradigm

            end

            if isempty(sp.st), continue, end

            % shift the spike times now, because we have shifted the nsEvents too when storing them above.
            %sp.st = sp.st + shiftSpikeTime;

            if keepMU, inds = sp.cgs<=3;
            else,      inds = sp.cgs==2;
            end

            cids = sp.cids(inds);
            cgs  = sp.cgs(inds);

            dataStruct(sess).data.(paradigms{par}).units.cluster_id = cids;
            dataStruct(sess).data.(paradigms{par}).units.cluster_type = cgs;

            dataStruct(sess).data.(paradigms{par}).units.cluster_labels = {'MU','SU','UN'};

            if exist('unitInfo','var')
                try
                    keepUnits = ismember(unitInfo.cluster_id,sp.cids);
                    depth     = unitInfo.depth(keepUnits);
                    ch        = unitInfo.ch(keepUnits);
                    nspks     = unitInfo.n_spikes(keepUnits);

                    MDI_depth = info.depths{1}(theseFiles(1));
                    if contains(info.probe_type,'DBC')
                        probe = ['DBC' info.probe_ID{1}(1:5)];
                        ch_depth  = calcProbeChDepth(MDI_depth,depth,probe);
                    elseif contains(info.probe_type,'Single')
                        ch_depth = MDI_depth;
                    end


                    dataStruct(sess).data.(paradigms{par}).units.depth = ch_depth;
                    dataStruct(sess).data.(paradigms{par}).units.ch    = depth;

                catch
                    fprintf('Issue with unitInfo...\n')

                end

            end

            fprintf('Adding %d SU, %d MU, %d unsorted\n\n',sum(cgs==2),sum(cgs==1),sum(cgs==3))

            % add each unit's spikes to an entry in spiketimes cell
            for unit=1:sum(inds)
                theseSpikes = sp.clu==cids(unit) & thisParSpikes;
                %                 theseSpikes = sp.clu==cids(unit);
                dataStruct(sess).data.(paradigms{par}).units.spiketimes{unit} = sp.st(theseSpikes);
            end

        end


        filename = sprintf('%s%ddots3DMPevents_%d.mat',info.subject,info.date,unique_sets(u));
        folder = fullfile(localDir,'rec_events');

        % check if CSV already exists, overwrite set off, and all parsspecified
        if (overwriteEventSets || ~exist(fullfile(folder,filename),'file'))
            
            try
                S = createNeuralEvents_oneStruct(dataStruct(sess).data);
                fprintf('saving events file %s\n',filename)
                save(fullfile(folder,filename),'S');
            catch
                fprintf('could not save events file %s\n',filename);
            end
        end
    end        
end

file = [subject '_' num2str(dateRange(1)) '-' num2str(dateRange(end)) '_neuralData_' area '.mat'];

disp('saving...');
save([localDir(1:length(localDir)-length(subject)-7) file], 'dataStruct');