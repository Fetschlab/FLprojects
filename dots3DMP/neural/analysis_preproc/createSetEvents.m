%% create event files by set

% this is 

% fields in 'events' which contain event times, this will be important later
tEvs   = {'trStart','fpOn','fixation','reward','stimOn','stimOff','saccOnset',...
    'targsOn','targHold','postTargHold','reward','breakfix','nexStart','nexEnd','return'};

sflds = {'subject','date','pen','gridxy','probe_type','probe_ID'};
for n = 1:length(currentFolderList)
%     disp(currentFolderList{n})
    if isempty(strfind(currentFolderList{n},'20')) || contains(currentFolderList{n},'Impedance'); continue; end
    
    clear info
    load(fullfile(localDir,[subject currentFolderList{n} 'dots3DMP_info.mat']));
    fprintf('%s, %d of %d\n',currentFolderList{n},n,length(currentFolderList))
    
    % we want 1 row in dataStruct for each unique 'recording set'
    [unique_sets,~,ic] = unique(info.rec_group);
    
    for u=1:length(unique_sets)

        filename = sprintf('%s%ddots3DMPevents_%d.mat',info.subject,info.date,unique_sets(u));

        if (overwriteEventSets || ~exist(fullfile(folder,filename),'file'))
            mountDir = sprintf('/Volumes/homes/fetschlab/data/%s/%s_neuro/%d/%s%d_%d/',subject,subject,info.date,subject,info.date,unique_sets(u));

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
                    catch, fprintf('PDS file %s not found..are you connected to the NAS?\n',PDSfilenames{pf});
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

                    % pull in relevant condition data from PLDAPS and sub-select trials from this paradigm

                    [thisParEvents]   = nsEventConditions(nsEvents,allPDS,lower(paradigms{par})); % % SJ added 08-22-2022 oneTargChoice and Conf!
                    timeStampsShifted = thisParEvents.analogInfo.timeStampsShifted ./ double(thisParEvents.analogInfo.Fs);

                    % do some concatenation in pldaps and events fields, in case the same par+block is split over multiple trellis files

                    nTr    = length(thisParEvents.Events.trStart);

                    % add all the fields in events and pldaps to dataStruct
                    % if event field is a time, shift it as needed
                    fnames = fieldnames(thisParEvents.Events);
                    for f=1:length(fnames)

                        if ismember(fnames{f},tEvs)
                            thisParEvents.(fnames{f}) = thisParEvents.Events.(fnames{f})  + timeStampsShifted(1);
                        end
                        data.(paradigms{par}).events.(fnames{f})(1,currPos+1:currPos+nTr) = thisParEvents.Events.(fnames{f});
                    end

                    fnames = fieldnames(thisParEvents.pldaps);
                    for f=1:length(fnames)
                        if strcmp(fnames{f},'unique_trial_number') && ~iscell(nsEvents.pldaps.unique_trial_number)
                            data.(paradigms{par}).pldaps.unique_trial_number(currPos+1:currPos+nTr,:) = thisParEvents.pldaps.(fnames{f});
                        else
                            data.(paradigms{par}).pldaps.(fnames{f})(1,currPos+1:currPos+nTr) = thisParEvents.pldaps.(fnames{f});
                        end
                    end

                    data.(paradigms{par}).pldaps.blockNum2(1,currPos+1:currPos+nTr) = utf;
                    currPos = currPos+nTr;

                end

            end


            % check if CSV already exists, overwrite set off, and all parsspecified


            S = createNeuralEvents_oneStruct(data);
            folder = fullfile(localDir,'rec_events');
            fprintf('saving events file %s\n',filename)
            if ~exist(folder,'dir')
                mkdir(folder)
            end
            save(fullfile(folder,filename),'S');
        else
            fprintf('events file %s already exists...skipping\n',filename)

        end
    end        
end