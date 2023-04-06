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

% added 04/2023
sess_info = readtable('/Users/stevenjerjian/Desktop/FetschLab/Analysis/RecSessionInfo.xlsx', sheet = subject);
sess_info.Properties.VariableNames = lower(sess_info.Properties.VariableNames);
sess_info.chs = table2cell(rowfun(@(x,y) x:y, sess_info(:,{'min_ch','max_ch'})));


% dataStruct = struct();
% sess = 0;
dataStruct = table2struct(sess_info);

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
        clear sp

%         sess = sess+1; % increment the row in dataStruct
        sess = find(sess_info.date == datetime(num2str(info.date),'InputFormat','yyyyMMdd') & sess_info.rec_set==unique_sets(u));

        remoteDirSpikes = sprintf('/var/services/homes/fetschlab/data/%s/%s_neuro/%d/%s%d_%d/',subject,subject,info.date,subject,info.date,unique_sets(u));
        mountDir = sprintf('/Volumes/homes/fetschlab/data/%s/%s_neuro/%d/%s%d_%d/',subject,subject,info.date,subject,info.date,unique_sets(u));

        if contains(info.probe_type{1},'Single')
            mountDir = [mountDir 'phy_WC/']; 
            continue
        end

        try
            disp(mountDir)
%             params = struct('excludeNoise',0);
            sp = loadKSdir(mountDir);
        catch
            sp.st = [];
            warning('dots3DMP:createSessionData:loadKSdir','Could not load kilosort sp struct for %d, set %d\n\n',info.date,unique_sets(u));
        end
        try
            unitInfo = getUnitInfo(mountDir, keepMU);
        catch
            warning('dots3DMP:createSessionData:getUnitInfo','Could not get cluster info for this ks file..file has probably not been manually curated\n')
        end

%         dataStruct(sess).date = info.date;
%         dataStruct(sess).info = info;
%         dataStruct(sess).set = unique_sets(u);
        
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
            
            % for each trellis file within a given set+paradigm, concatenate and store the events

            [unique_trellis_files,~,ii] = unique(info.trellis_filenums(theseFiles));
            
            thisParSpikes  = false(size(sp.st));

            % now loop over each trellis file within a particular paradigm
            currPos = 0;
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

                if ~isempty(sp.st)
                    % pick out spikes from the sp.st vector which can be linked to this paradigm's timeframe (with a reasonable buffer on either
                    % side, e.g. 20secs), and what time to shift the spike times by (if any) so that they align with events again
                    % Note. this shift is only necessary if multiple Trellis recordings were made for same location - these will
                    % have been concatenated for kilosort sorting, but events will still be separate

                    timeLims = timeStampsShifted(1) + thisParEvents.Events.trStart([1 end]) + [-1 1]*20;

                    thisFileSpikes = (sp.st >= timeLims(1) & sp.st < timeLims(2));
                    thisParSpikes  = thisParSpikes | thisFileSpikes; % union

                end

            end

            if isempty(sp.st), continue, end

            % shift the spike times now, because we have shifted the nsEvents too when storing them above.
            %sp.st = sp.st + shiftSpikeTime;

            if keepMU, inds = sp.cgs<=3;
            else,      inds = sp.cgs==2;
            end

            if exist('unitInfo','var')
                try
                    keepUnits = ismember(unitInfo.cluster_id,sp.cids)';
                    depth     = unitInfo.depth(keepUnits)';
                    ch        = unitInfo.ch(keepUnits)';
                    nspks     = unitInfo.n_spikes(keepUnits)';

                    % MDI_depth = info.depths{1}(theseFiles(1));
                    MDI_depth = dataStruct(sess).mdi_depth_um;
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

            inds = inds & ismember(depth, dataStruct(sess).chs);

            cids = sp.cids(inds);
            cgs  = sp.cgs(inds);

            dataStruct(sess).data.(paradigms{par}).units.cluster_id = cids;
            dataStruct(sess).data.(paradigms{par}).units.cluster_type = cgs;

            dataStruct(sess).data.(paradigms{par}).units.cluster_labels = {'MU','SU','UN'};

            fprintf('Adding %d SU, %d MU, %d unsorted\n\n',sum(cgs==2),sum(cgs==1),sum(cgs==3|cgs==0))

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

file = [subject '_' num2str(dateRange(1)) '-' num2str(dateRange(end)) '_neuralData.mat'];

disp('saving...');
save([localDir(1:length(localDir)-length(subject)-7) file], 'dataStruct');