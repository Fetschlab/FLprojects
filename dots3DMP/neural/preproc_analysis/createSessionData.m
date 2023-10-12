%% create SessionData

% create a struct with neural Data for each session
% see dots3DMP_NeuralPreProcessing for explanation of dataStruct structure

% RECORDING SET NOTES
%
% Successive recordings with the same putative units (i.e. at the same location) are concatenated where necessary, since Kilosort spike sorting on the combined data with kilosort 
% is much preferred, rather than trying to reconcile cluster ids post-hoc.
% The convention throughout the codebase is to refer to such 'combined' recordings as recording 'sets'.
% If recordings are concatenated, the de facto timestamps (0-len) are shifted according to the length of the concatenated set so that the range of events and spikes is
% matched for a given recording (nsEvents.analogData.timeStamps and .timeStampsShifted).
%
% e.g. if a set is comprised of two recordings, with timestamp ranges 0-->N and 0-->M, respectively the final timestamps will run from 0 --> N+M. 
% Spiketimes and events from recording 1 will be stored as is, spiketimes and events from recording 2 will be re-assigned to their existing value plus N.
% 
% If recording is left running across experiments, the set will be 1-to-1 matched to recording files, rendering the above shifting superfluous. 
% However, multiple PLDAPS files will be associated with one recording file, which requires other considerations.

% the loop structure is:
%   folder (date)
%       ---> set
%           ---> paradigm
%               ---> trellis files
%                   

%% initialize dataStruct

% added 04/2023 SJ
% externally defined excel file containing recording metadata
% some information is redundant with individual info mat files, but this
% data structure makes it easier to control which files become part of data
% struct (e.g. probe recordings only), and splits by area/probe

% this could also be achieved by a loop over info.probe, but might require
% some refactoring. Works fine for now.

sess_info = readtable(sess_info_file, sheet = subject);
sess_info.Properties.VariableNames = lower(sess_info.Properties.VariableNames);
sess_info.chs = table2cell(rowfun(@(x,y) x:y, sess_info(:,{'min_ch','max_ch'})));
sess_info = sess_info(logical(sess_info.is_good),:);

% dataStruct = struct(); sess = 0; % original
dataStruct = table2struct(sess_info);

%% main loop

% loop over each 'date' in folder list, then over unique sets
% each date/set is referenced against the sess_info sheet

% fields in 'events' which contain event times, this will be important later
tEvs   = {'trStart','fpOn','fixation','reward','stimOn','stimOff','saccOnset',...
    'targsOn','targHold','postTargHold','reward','breakfix','nexStart','nexEnd','return'};

clus_labels = {'MU','SU','UN'};


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

%         sess = sess+1; % increment the row in dataStruct % old
        sess = find(sess_info.date == datetime(num2str(info.date),'InputFormat','yyyyMMdd') & sess_info.rec_set==unique_sets(u));

        if sum(sess)==0
            fprintf('no date/set match in RecSessionInfo.xlsx\n')
            continue
        end
        
        remoteDirSpikes = sprintf('/var/services/homes/fetschlab/data/%s/%s_neuro/%d/%s%d_%d/',subject,subject,info.date,subject,info.date,unique_sets(u));
        mountDir = sprintf('/Volumes/homes/fetschlab/data/%s/%s_neuro/%d/%s%d_%d/',subject,subject,info.date,subject,info.date,unique_sets(u));

        % if recording was single electrode, sorting was done with SI, so sub-folder is phy_WC
        if contains(info.probe_type{1},'Single')
            mountDir = [mountDir 'phy_WC/']; 
            continue % skip these regardless for now
        end


        try
            disp(mountDir)
            sp = loadKSdir(mountDir);
        catch
            sp.st = [];
            error('dots3DMP:createSessionData:loadKSdir','Could not load kilosort sp struct for %d, set %d...Are you connected to the NAS?\n\n',info.date,unique_sets(u));
        end

        try
            unitInfo = readtable(fullfile(mountDir,'cluster_info.tsv'),'FileType','delimitedtext');
        catch
            error('dots3DMP:createSessionData:getUnitInfo','Could not get cluster info for this ks file..file has probably not been manually curated\n')
        end

        % old SJ 04/2023
%         dataStruct(sess).date = info.date;
%         dataStruct(sess).info = info;
%         dataStruct(sess).set = unique_sets(u);
        

        % loop over paradigms 
        % (NOTE that the logic here differs from how nsEvents is initially created on the experiment rig)
        for par=1:length(paradigms)
            
            theseFiles = find((ic'==unique_sets(u)) & ismember(lower(info.par),lower(paradigms{par})) & (~isnan(info.pldaps_filetimes)));
            if isempty(theseFiles), continue, end
                        
            % concatenate all the data and condition fields of PDS files for given paradigm which are marked as part of the same recording 'set'
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

            % now loop over each trellis file

            currPos = 0; % current position in each event field vector, for concat purposes
            for utf=1:length(unique_trellis_files)
                NSfilename  = sprintf('%s%ddots3DMP%04d_RippleEvents.mat',info.subject,info.date,unique_trellis_files(utf));
                
                % messed up with PDS files on this one, oops
                if strcmp(NSfilename, 'lucio20220719dots3DMP0008_RippleEvents.mat'), continue, end

                try
                    load(fullfile(localDir,NSfilename));
                    fprintf('adding data from %s (%s)\n\n', NSfilename, paradigms{par})
                catch
                    fprintf('Could not load %s, skipping...\n\n', NSfilename)
                    continue
                end
                

                % allPDS should now match nsEvents, so pull in relevant condition data from PLDAPS and sub-select trials from this paradigm
                [thisParEvents]   = nsEventConditions(nsEvents,allPDS,lower(paradigms{par})); % SJ added 08-22-2022 oneTargChoice and Conf!
                timeStampsShifted = thisParEvents.analogInfo.timeStampsShifted ./ double(thisParEvents.analogInfo.Fs);

                % add all the fields in events to dataStruct
                % if event field is a time, shift it as needed
                nTr    = length(thisParEvents.Events.trStart);
                fnames = fieldnames(thisParEvents.Events);         
                for f=1:length(fnames)
                    
                    if ismember(fnames{f},tEvs)
                        thisParEvents.(fnames{f}) = thisParEvents.Events.(fnames{f})  + timeStampsShifted(1); % SHIFT
                    end
                    for s = 1:length(sess)
                        dataStruct(sess(s)).data.(paradigms{par}).events.(fnames{f})(1,currPos+1:currPos+nTr) = thisParEvents.Events.(fnames{f});
                    end
                end

                % repeat for pldaps field
                fnames = fieldnames(thisParEvents.pldaps);
                for f=1:length(fnames)
                    for s = 1:length(sess)
                        if strcmp(fnames{f},'unique_trial_number') && ~iscell(nsEvents.pldaps.unique_trial_number)
                            dataStruct(sess(s)).data.(paradigms{par}).pldaps.unique_trial_number(currPos+1:currPos+nTr,:) = thisParEvents.pldaps.(fnames{f});
                        else
                            dataStruct(sess(s)).data.(paradigms{par}).pldaps.(fnames{f})(1,currPos+1:currPos+nTr) = thisParEvents.pldaps.(fnames{f});
                        end
                    end
                end
                currPos = currPos+nTr; % update


                if ~isempty(sp.st)
                    % mask sp.st for spikes which occured within this paradigm's timeframe (with a reasonable buffer on either
                    % side, e.g. 20secs). This is not fundamental, but saves storing the same set of spikes across all paradigms within each individual paradigm's field
                    
                    timeLims = timeStampsShifted(1) + thisParEvents.Events.trStart([1 end]) + [-1 1]*20;

                    thisFileSpikes = (sp.st >= timeLims(1) & sp.st < timeLims(2));
                    thisParSpikes  = thisParSpikes | thisFileSpikes;

                end

            end

            % save the units data, with relevant ch and depth information

            if isempty(sp.st), continue, end

            keepUnits = ismember(unitInfo.cluster_id,sp.cids)';
            depth     = unitInfo.depth(keepUnits)';
            ch        = unitInfo.ch(keepUnits)';
            nspks     = unitInfo.n_spikes(keepUnits)';

            % loop over sessions
            % a recording with two probes will have two separate rows 
            for s = 1:length(sess)

                % convert matlab datetime to string rep
                dataStruct(sess(s)).date = datestr(dataStruct(sess(s)).date);


                if contains(info.probe_type,'DBC')
                    ch_depth  = calcProbeChDepth(depth,dataStruct(sess(s)));
                elseif contains(info.probe_type,'Single')
                    ch_depth = MDI_depth;
                end

                if keepMU, inds = sp.cgs<=3;
                else,      inds = sp.cgs==2;
                end


                dataStruct(sess(s)).data.(paradigms{par}).units.depth = ch_depth(inds);
                dataStruct(sess(s)).data.(paradigms{par}).units.ch    = depth(inds);

                inds = inds & ismember(ch+1, dataStruct(sess(s)).chs);

                cids = sp.cids(inds);
                cgs  = sp.cgs(inds);

                dataStruct(sess(s)).data.(paradigms{par}).units.cluster_id = cids;
                dataStruct(sess(s)).data.(paradigms{par}).units.cluster_type = cgs;

                dataStruct(sess(s)).data.(paradigms{par}).units.cluster_labels = clus_labels(cgs);
                dataStruct(sess(s)).data.(paradigms{par}).units.npks = nspks(inds);

                if any(cgs==3)
                    fprintf('Adding %d SU, %d MU, %d unsorted\n\n',sum(cgs==2),sum(cgs==1),sum(cgs==3|cgs==0))
                else
                    fprintf('Adding %d SU, %d MU\n\n',sum(cgs==2),sum(cgs==1))
                end

                % finally, add each unit's spikes to an entry in spiketimes cell
                for unit=1:sum(inds)
                    theseSpikes = sp.clu==cids(unit) & thisParSpikes;
                    %                 theseSpikes = sp.clu==cids(unit);
                    dataStruct(sess(s)).data.(paradigms{par}).units.spiketimes{unit} = sp.st(theseSpikes);
                end
            end
        end

        %{ 
        % 04/06/2023 SJ defunct for now
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
        %}
    end        
end


% SAVE IT!!
file = [subject '_' num2str(dateRange(1)) '-' num2str(dateRange(end)) '_neuralData.mat'];
disp('saving...');
save([localDir(1:length(localDir)-length(subject)-7) file], 'dataStruct');