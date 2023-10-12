% PLDAPS_preprocessing_nolocal


% ditch the overflow of rmfield calls in pdsCleanup 
% https://undocumentedmatlab.com/articles/rmfield-performance


today    = str2double(datestr(now,'yyyymmdd'));



subject = 'lucio';
paradigm = 'dots3DMP';

dateRange = 20220512:20230601;

% where dataStruct should be saved
localDir = '/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/';

opts = struct('eyemovement', 1);
data = pds_preprocessing_dots3DMP(subject, paradigm, dateRange, localDir, opts);

filename = [subject '_' num2str(dateRange(1)) '-' num2str(dateRange(end)) '_wEM.mat'];
save(fullfile(localDir, filename), 'data');


function data = pds_preprocessing_dots3DMP(subject, paradigm, dateRange, localDir, options)


    % TODO consider providing existing file to append to...

    if nargin < 5, options = struct(); end
    if nargin < 4, localDir = pwd; end

    opt_fields = {'nexonar', 'dotposition', 'eyemovement', 'reward', 'eventtiming', ...
        'useSCP', 'useVPN'};
    
    for f = 1:length(opt_fields)
        if ~isfield(options, opt_fields{f})
            options.(opt_fields{f}) = 0;
        end
    end

    if options.eyemovement || options.nexonar
        options.eventtiming = 1;
    end

    if options.eventtiming
        pt_names = {'timeTargFixEntered', 'timeConfTargEntered', 'timeToConfidence'};
    
        st_names = {'timeFpEntered', 'timeTargetOn', 'timeTargetOff', ...
            'timeMotionStateBegin', 'timeMotionDone', 'timeGoCue', ...
            'TimeTargEntered', 'motionStartTime', 'timeBreakFix', 'timeComplete', ...
            'delayToGoCue', 'delayToDots', 'holdAfterDotsOnset', 'timeLastFrame'};
    end
    bhv_names = {'choice', 'RT', 'PDW', 'TargMissed', 'oneTargChoice', 'oneTargConf'};


    % set directories
    if isunix
        remoteDir = ['/var/services/homes/fetschlab/data/' subject '/'];
        mountDir  = ['/Volumes/homes/fetschlab/data/' subject '/'];
    else
        % TODO
        % remoteDir = ...
    end

    localDir = fullfile(localDir, subject);
    if ~exist(localDir,'dir')
        mkdir(localDir);
    end


    % initialize data struct
    data = struct;
    data.filename = {};
    data.subj = {};
    data.date = [];
    data.choice = []; % initialize this one field, to keep counter of numtrials


    %% get data from NAS
    % requires password(s) to be entered if a auth key not available
    % https://www.howtogeek.com/66776/how-to-remotely-copy-files-over-ssh-without-entering-your-password/
    
    % VERSION 3.0: 10-06-2023 SJ
    % loading directly from mounted NAS each time might 
    
    % get file list from remote dir
    if ~options.useVPN
        ip_add = 'fetschlab@172.30.3.33'; % MBI 
    else
        ip_add = 'fetschlab@10.161.240.133'; % probably off campus, try proxy IP (requires VPN)
    end

    cmd = ['ssh ' ip_add ' ls ' remoteDir];

    [~,remoteFileList] = system(cmd);  % system(cmd, "-echo") % print output
    if any(strfind(remoteFileList,'timed out')); error(remoteFileList); end
    
    remoteFiles = splitlines(remoteFileList);
    remoteFiles = remoteFiles(contains(remoteFiles, subject) & ~contains(remoteFiles, '_'));
    [~, remoteFileNames, ~] = cellfun(@fileparts, remoteFiles, 'UniformOutput', false);

    dateStart = length(subject)+1;

    for f = 1:length(remoteFiles)
        thisDate = str2double(remoteFiles{f}(dateStart:dateStart+7)); % yyyymmdd date format
        thisPar  = remoteFiles{f}(dateStart+8:length(remoteFileNames{f})-4);
%         thisPar  = remoteFiles{f}(dateStart+8:length(remoteFiles{f})-8);

        if ~any(strcmp(data.filename, remoteFileNames{f})) && any(dateRange == thisDate) && strcmp(thisPar, paradigm) 
            % ok, we want this file

            T = length(data.choice);

            try

                tstart = tic;

                if options.useSCP
                    cmd = ['scp -r ' ip_add ':' remoteDir remoteFiles{f} ' ' localDir];
                    system(cmd, "-echo")
                    load(fullfile(localDir, remoteFiles{f}), '-mat', 'PDS')
                    % TODO, delete file afterwards, or save reduced
                    % version...
                else
                    load(fullfile(mountDir, remoteFiles{f}), '-mat', 'PDS');
                end


                % SJ added 02/2022, to generate 3DMP dots offline from trialSeeds, no
                % need to save online for storage space reasons
                if options.dotposition
                    try
                        if ~isfield(PDS.data{1}.stimulus,'dotX_3D')
                            [dotX_3D,dotY_3D,dotZ_3D,dotSize] = generateDots3D_offline(PDS);
                        end
                    catch
                        disp("offline dot generation did not work...")
                    end
                end

                % extract desired fields for each trial, append to struct
                for t = 1:length(PDS.data)


                    if isfield(PDS.data{t}.behavior,'choice') % and save out the data, excluding trials with
                        % missing data (choice is a good marker for this)
                        T = T+1; % increment trial counter

                        data.trialNum(T,1) = t;

                        data.filename{T,1} = remoteFileNames{f};
                        data.date(T,1) = thisDate;

                        if contains(subject,'human')
                            data.subj{T,1} = remoteFileNames{f}(dateStart(1)-3:dateStart(1)-1); % 3-letter code
                        else
                            data.subj{T,1} = subject;
                        end

                        % independent variables (conditions) are stored in PDS.conditions.stimulus
                        fnames = fieldnames(PDS.conditions{t}.stimulus);
                        for F = 1:length(fnames)
                            data.(fnames{F})(T,1) = PDS.conditions{t}.stimulus.(fnames{F});
                        end

                        % dependent variables (outcomes) stored in PDS.data.behavior
                        for F = 1:length(bhv_names)
                            data.(bhv_names{F})(T,1) = PDS.data{t}.behavior.(bhv_names{F});
                        end

                        % misc
                        try
                            data.iTrial(T,1) = PDS.conditions{t}.pldaps.iTrial;
                        catch
                            data.iTrial(T,1) = NaN;
                        end
                        
                        % options
                        if options.reward
                            fnames = fieldnames(PDS.data{t}.reward);
                            for F = 1:length(fnames)
                                data.(fnames{F})(T,1) = PDS.data{t}.reward.(fnames{F});
                            end
                        end

                        if options.eyemovement
                            try
                                data.ADCdata{T, 1} = PDS.data{t}.datapixx.adc.data;
                                dp_time = PDS.data{t}.datapixx.unique_trial_time(2);
                                data.ADCtime{T, 1} = PDS.data{t}.datapixx.adc.dataSampleTimes - dp_time;
                            catch
                                data.ADCdata{T, 1} = NaN;
                                data.ADCtime{T, 1} = NaN;
                            end
                        end

                        if options.eventtiming
                            for F = 1:length(pt_names)
                                if isfield(PDS.data{t}.postTarget, pt_names{F})
                                    data.(pt_names{F})(T,1) = PDS.data{t}.postTarget.(pt_names{F});
                                else
                                    data.(pt_names{F})(T,1) = NaN;
                                end
                            end

                            for F = 1:length(st_names)
                                if isfield(PDS.data{t}.stimulus, st_names{F})
                                    data.(st_names{F})(T,1) = PDS.data{t}.stimulus.(st_names{F});
                                else
                                    data.(st_names{F})(T,1) = NaN;
                                end
                            end
                        end

                        if options.dotposition
                            try
                                if ~isfield(PDS.data{t}.stimulus,'dotX_3D')
                                    PDS.data{t}.stimulus.dotX_3D = dotX_3D{t};
                                    PDS.data{t}.stimulus.dotY_3D = dotY_3D{t};
                                    PDS.data{t}.stimulus.dotZ_3D = dotZ_3D{t};
                                    PDS.data{t}.stimulus.dotSize = dotSize{t};
                                end
                                data.dotX_3D = PDS.data{t}.stimulus.dotX_3D;
                                data.dotY_3D = PDS.data{t}.stimulus.dotY_3D;
                                data.dotZ_3D = PDS.data{t}.stimulus.dotZ_3D;
                                data.dotSize = PDS.data{t}.stimulus.dotSize;
                            catch
                            end
                        end

                        if options.nexonar
                            disp('not yet implemented')
                        end
                

                    end
                end


                telapsed = toc(tstart);
                fprintf('%s, time taken: %.2fs', remoteFileNames{f}, telapsed)

            catch
                keyboard
                warning(['Processing issue, or could not load ' remoteFiles{f} '. File may be corrupt -- skipping']);
            end
         
            fprintf('\ncumulative trials = %d\n', T);
        end
    end
    disp('done')
end
        
    
    

