% PLDAPS_preprocessing_nolocal


% ditch the overflow of rmfield calls in pdsCleanup 
% https://undocumentedmatlab.com/articles/rmfield-performance

today    = str2double(datestr(now,'yyyymmdd'));

subject = 'lucio';
paradigm = 'dots3DMP';

dateRange = 20220512:20230601;

% where dataStruct should be saved
% localDir = '/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/';
localDir = uigetdir([], 'Choose directory to save data to');

opts = struct('eyemovement', 1);
data = pds_preprocessing_dots3DMP(subject, paradigm, dateRange, localDir, [], opts);

filename = [subject '_' num2str(dateRange(1)) '-' num2str(dateRange(end)) '_wEM.mat'];
save(fullfile(localDir, filename), 'data');


function [data] = pds_preprocessing_dots3DMP(subject, paradigm, dateRange, localDir, data, options)

    if nargin < 6, options = struct(); end
    if nargin < 4 || isempty(localDir), localDir = uigetdir([], 'Choose directory to save data to'); end

    use_existing_file = ~(nargin < 5 || isempty(data));

    if ~isfield(options, 'eventtiming')
        options.eventtiming = options.eyemovement || options.nexonar;
    end

    opt_fields = {'nexonar', 'dotposition', 'eyemovement', 'reward', 'eventtiming', 'useSCP', 'useVPN'};
    for f = 1:length(opt_fields)
        if ~isfield(options, opt_fields{f})
            options.(opt_fields{f}) = 0; % default = False
        end
    end
    

    if options.eyemovement
        % for aligning and cutting each datapixx trial to a 'time zero'
        align_event = 'timeGoCue';
        t_range = [-0.5, 1]; % seconds
        Fs = 1000;
        t_range = round(t_range*Fs);

        eyewin_samples = (t_range(1):t_range(2)-1);

    end
        
    if options.eventtiming
        pt_names = {'timeTargFixEntered', 'timeConfTargEntered', 'timeToConfidence'};
    
        st_names = {'timeFpEntered', 'timeTargetOn', 'timeTargetOff', ...
            'timeMotionStateBegin', 'timeMotionDone', 'timeGoCue', ...
            'TimeTargEntered', 'motionStartTime', 'timeBreakFix', 'timeComplete', ...
            'delayToGoCue', 'delayToDots', 'holdAfterDotsOnset', 'timeLastFrame'};
    else
        pt_names = {};
        st_names = {};
    end
    bhv_names = {'choice', 'RT', 'PDW', 'TargMissed', 'oneTargChoice', 'oneTargConf'};


    % set directories
    [~, hostname] = system('hostname');
    hostname = strip(hostname);

    if ismac % someone's local Mac
        remoteDir = fullfile('/var/services/homes/fetschlab/data/', subject);
        mountDir  = fullfile('/Volumes/homes/fetschlab/data/', subject);
    elseif strcmp(hostname, 'DESKTOP-JRJ0F9N')
        remoteDir = fullfile('Z:\fetschlab\data', subject);
        mountDir  = remoteDir;
    end

    localDir = fullfile(localDir, subject);
    if ~exist(localDir,'dir'), mkdir(localDir); end


    % initialize data struct
    if ~use_existing_file

        data.filename = {};
        data.subj = {};
        data.date = [];

%         ntr_est = 3e5;
%         data.filename = cell(ntr_est, 1);
%         data.subj = cell(ntr_est, 1);
%         data.date = nan(ntr_est, 1);
% 
%         all_flds = unique([bhv_names, pt_names, st_names, {'iTrial', 'trialNum'}]);
%        
%         for f = 1:length(all_flds)
%             data.(all_flds{f}) = nan(ntr_est, 1);
%             
%         end
%         if options.eyemovement
%             data.ADCdata = cell(ntr_est, 1);
%             data.ADCtime = cell(ntr_est, 1);
%         end

        T = 0;

    else
        % use existing fields..NEED TO ADD CHECKS BELOW TO ONLY ADD NEW
        % DATA IF FIELD EXISTS
        all_flds = fieldnames(data);
        T = find(isnan(data.date), 1, 'first');
    end


    %% get data from NAS
    % requires password(s) to be entered if a auth key not available
    % https://www.howtogeek.com/66776/how-to-remotely-copy-files-over-ssh-without-entering-your-password/
    
    % VERSION 3.0: 10-06-2023 SJ
    
    % get file list from remote dir
    if ~options.useVPN
        ip_add = 'fetschlab@172.30.3.33'; % MBI 
    else
        ip_add = 'fetschlab@10.161.240.133'; % probably off campus, try proxy IP (requires VPN)
    end

    
    if strcmp(hostname, 'DESKTOP-JRJ0F9N')
        remoteFiles = dir([mountDir, '/*.mat']); % ignore non .mat files (including old .PDS)
        remoteFiles = {remoteFiles.name}';
        options.useSCP = 0; % force to 0
        options.useVPN = 0; 
    else
        cmd = ['ssh ' ip_add ' ls ' remoteDir];
        [~,remoteFileList] = system(cmd, "-echo");  % system(cmd, "-echo") % print output
        if any(strfind(remoteFileList,'timed out')); error(remoteFileList); end
        remoteFiles = splitlines(remoteFileList);
        remoteFiles = remoteFiles(contains(remoteFiles, subject) & ~contains(remoteFiles, '_'));
    end
    [~, remoteFileNames, ~] = cellfun(@fileparts, remoteFiles, 'UniformOutput', false);

    dateStart = length(subject)+1;

    for f = 1:length(remoteFiles)
        thisDate = str2double(remoteFiles{f}(dateStart:dateStart+7)); % yyyymmdd date format
        thisPar  = remoteFiles{f}(dateStart+8:length(remoteFileNames{f})-4);

        if ~any(strcmp(data.filename, remoteFileNames{f})) && any(dateRange == thisDate) && strcmp(thisPar, paradigm) 
            % ok, we want this file

            try

                tstart = tic;
                 fprintf('Loading file: %s\n', remoteFileNames{f})


                if options.useSCP
                    % copy file to local machine
                    % to save reduced version(?)
                    cmd = ['scp -r ' ip_add ':' remoteDir remoteFiles{f} ' ' localDir];
                    system(cmd, "-echo")
                    load(fullfile(localDir, remoteFiles{f}), '-mat', 'PDS')
                    
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
                    
                    if isfield(PDS.data{t}.behavior,'choice') % excluding trials with missing data

                        T = T+1; % increment trial counter

                        data.trialNum(T) = t;

                        data.filename{T} = remoteFileNames{f};
                        data.date(T) = thisDate;

                        if contains(subject,'human')
                            data.subj{T} = remoteFileNames{f}(dateStart(1)-3:dateStart(1)-1); % 3-letter code
                        else
                            data.subj{T} = subject;
                        end

                        % independent variables (conditions) are stored in PDS.conditions.stimulus
                        fnames = fieldnames(PDS.conditions{t}.stimulus);
                        for F = 1:length(fnames)
                            data.(fnames{F})(T) = PDS.conditions{t}.stimulus.(fnames{F});
                        end

                        % dependent variables (outcomes) stored in PDS.data.behavior
                        for F = 1:length(bhv_names)
                            data.(bhv_names{F})(T) = PDS.data{t}.behavior.(bhv_names{F});
                        end

                        % misc
                        try
                            data.iTrial(T) = PDS.conditions{t}.pldaps.iTrial;
                        catch
                            data.iTrial(T) = NaN;
                        end
                        
                        % options
                        if options.reward
                            fnames = fieldnames(PDS.data{t}.reward);
                            for F = 1:length(fnames)
                                data.(fnames{F})(T) = PDS.data{t}.reward.(fnames{F});
                            end
                        end

                        if options.eventtiming
                            for F = 1:length(pt_names)
                                if isfield(PDS.data{t}.postTarget, pt_names{F})
                                    data.(pt_names{F})(T) = PDS.data{t}.postTarget.(pt_names{F});
                                else
                                    data.(pt_names{F})(T) = NaN;
                                end
                            end

                            for F = 1:length(st_names)
                                if isfield(PDS.data{t}.stimulus, st_names{F})
                                    data.(st_names{F})(T) = PDS.data{t}.stimulus.(st_names{F});
                                else
                                    data.(st_names{F})(T) = NaN;
                                end
                            end
                        end

                        if options.eyemovement
                            try
                                adc_data = PDS.data{t}.datapixx.adc.data;
                                dp_time = PDS.data{t}.datapixx.unique_trial_time(2);
                                adc_time = PDS.data{t}.datapixx.adc.dataSampleTimes - dp_time;

                                [~, zero_pos] = min(abs(adc_time - data.(align_event)(t)));

                                data.eyeX(T, 1:length(win_samples)) = adc_data(1, eyewin_samples+zero_pos);
                                data.eyeY(T, 1:length(win_samples)) = adc_data(2, eyewin_samples+zero_pos);
                                data.ADCtime(T, 1:length(win_samples)) = adc_time(eyewin_samples+zero_pos);

                            catch
                                data.eyeX(T, 1:length(eyewin_samples)) = NaN;
                                data.eyeX(T, 1:length(eyewin_samples)) = NaN;
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
                                data.dotX_3D{T} = PDS.data{t}.stimulus.dotX_3D;
                                data.dotY_3D{T} = PDS.data{t}.stimulus.dotY_3D;
                                data.dotZ_3D{T} = PDS.data{t}.stimulus.dotZ_3D;
                                data.dotSize{T} = PDS.data{t}.stimulus.dotSize;
                            catch
                            end
                        end

                        if options.nexonar
                            disp('not yet implemented')
                            % need to move over from original
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
        
    

