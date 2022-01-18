%% createDataStructure:
% load local files and create a single struct with all the trials and
% fields, excluding fieldExcludes (check struct data after using and add
% excludes as needed)

% VERSION 2.0: 08-30-19 CF
% must be run by PLDAPS_preprocessing

data = struct;
data.filename = {};
data.subj = {};
data.choice = []; % initialize this one field, you'll see why

fieldExcludes = {'leftEarly','tooSlow','fixFP','FPHeld','eyeXYs','corrLoopActive','goodtrial', ...
                 'timeTargDisappears','probOfMemorySaccade','leftTargR','leftTargTheta', ...
                 'rightTargR','rightTargTheta','audioFeedback','textFeedback','rewardDelay','reward'};

% now search localDir again for matching files and extract the desired variables from PDS
allFiles = dir(localDir);
if addNexonarDataToStruct, allNexFiles = dir([localDirNex '/*.mat']); end

for d = 1:length(dateRange)
    for f = 3:length(allFiles) % skip 1+2, they are are "." and ".."
        if contains(allFiles(f).name, subject) ... % check if the target file is in localDir
           && contains(allFiles(f).name, num2str(dateRange(d))) ...
           && contains(allFiles(f).name, paradigm) ...
           && ~contains(allFiles(f).name, 'icloud')
       
            disp(['loading ' allFiles(f).name]);
            
            try
                

                load([localDir allFiles(f).name],'-mat'); % load it. this is the time limiting step;
                                                          % will eventually change how data are saved to make this faster
                if exist('PDS','var')
                  
                    T = length(data.choice); % set trial counter to continue where left off (0, to start)
                    fprintf('\ncumulative trials processed = %d\n',T);
                    for t = 1:length(PDS.data) % loop over trials for this file,
                        if isfield(PDS.data{t}.behavior,'choice') % and save out the data, excluding trials with
                                                                  % missing data (choice is a good marker for this)
                            T = T+1; % increment trial counter

                            data.filename{T,1} = allFiles(f).name(1:end-4);
                            dateStart = strfind(allFiles(f).name,'20');
                            if contains(subject,'human')
                                data.subj{T,1} = allFiles(f).name(dateStart(1)-3:dateStart(1)-1); % 3-letter code
                            else
                                data.subj{T,1} = subject;
                            end
                            
                            % independent variables are stored in PDS.conditions.stimulus
                            fnames = fieldnames(PDS.conditions{t}.stimulus);
                            fnames(ismember(fnames,fieldExcludes)) = [];
                            for F = 1:length(fnames)
                                eval(['data.' fnames{F} '(T,1) = PDS.conditions{t}.stimulus.' fnames{F} ';']);
                            end
                                % EXCEPT! dot duration, and anything else in
                                % PDS.data.stimulus, i.e. things that are generated
                                % on each trial and not pre-configured. For now
                                % will need to hard-code them here.
                            if isfield(PDS.data{t}.stimulus,'dotDuration')
                                data.duration(T,1) = PDS.data{t}.stimulus.dotDuration;
                            end
                            if isfield(PDS.data{t}.stimulus,'dotPos')
                                data.dotPos{T,1} = PDS.data{t}.stimulus.dotPos;
                            end
                            
                            % to store all the reward amounts
                            % maybe don't need this, all reward analyses
                            % will have to be done within session first, so
                            % go to cleaned-up PDS files
%                             if isfield(PDS.data{t},'reward')
%                                 fnames = fieldnames(PDS.data{t}.reward);
%                                 fnames(ismember(fnames,fieldExcludes)) = [];
%                                 for F = 1:length(fnames)
%                                     eval(['data.reward.' fnames{F} '(T,1) = PDS.data{t}.reward.' fnames{F} ';']);
%                                 end
%                             end
                            
                            % SJ 10/2021 not true 'RT', need to replace
                            % with time that eyes leave choice target?
%                             if isfield(PDS.data{t},'postTarget')
%                                 try
%                                     data.confRT(T,1) = PDS.data{t}.postTarget.timeConfTargEntered - PDS.data{t}.postTarget.timeToConfidence;
%                                 catch
%                                     data.confRT(T,1) = NaN;
% %                                     data.confRT(T,1) = PDS.data{t}.postTarget.timeConfTargEntered - PDS.data{t}.postTarget.timeToConfidence;
%                                 end
%                             end
                            
                            % dependent variables are stored in PDS.data.behavior
                            fnames = fieldnames(PDS.data{t}.behavior);
                            fnames(ismember(fnames,fieldExcludes)) = [];
                            for F = 1:length(fnames)
                                % SJ 07-2020, correct defaults to logical but
                                % then gives error for NaN - use double instead
                                if strcmp(fnames{F},'correct'), data.correct(T,1) = 0; end
                                eval(['data.' fnames{F} '(T,1) = PDS.data{t}.behavior.' fnames{F} ';']);
                            end
                                                        
                            % noticed a couple extra things we need, not in either place -CF 02-2021
                            try
                                data.oneTargPDW(T,1) = PDS.data{t}.postTarget.markOneConf;
                            catch
                            end
                            try
                                data.delayToPDW(T,1) = PDS.data{t}.postTarget.delayToConfidence;
                            catch
                            end
                            
                        end
                    end
                    
                    % SJ 08-2021, adding Nexonar data in at this point
                    if addNexonarDataToStruct
                        matchingNexFile = cellfun(@(x) strcmp(x(1:25), allFiles(f).name(1:25)), {allNexFiles.name});
                        if sum(matchingNexFile)==1
                            load([localDirNex allNexFiles(matchingNexFile).name],'-mat'); % load the nexonar data
                            
                            nexPDS = dots3DMP_nexonarCleanUp(nex,PDS);
                            data.nexonar = nexPDS';

                        else
                            % no matching nexonar data (not recorded?)
                            % or more than one matching file (but this should be impossible)
                            disp(['could not find matching nexonar data for ' allFiles(f).name '...skipping'])
                        end
                    end
                    
                    clear PDS
                end
            
            catch me
                warning(['Processing issue, or could not load ' allFiles(f).name '. File may be corrupt -- skipping']);
            end

        end
    end
end


%% some bookkeeping: manually merge any vars whose names have changed, etc

if strcmp(subject(1:5),'human')
    data.conf = data.saccEndPoint;
else
    % do we still need this? SJ 07-2020
    if isfield(data,'postDecisionConfidence')
    if isfield(data,'PDW')
        if length(data.postDecisionConfidence)<length(data.PDW) && length(data.PDW)==length(data.choice)
            data.PDW(1:length(data.postDecisionConfidence)) = data.postDecisionConfidence;
            data = rmfield(data,'postDecisionConfidence');
        else
            error('unsure, diagnose by looking at data');
        end
    elseif length(data.postDecisionConfidence)==length(data.choice)
        data.PDW = data.postDecisionConfidence;
        data = rmfield(data,'postDecisionConfidence');
    else
        error('unsure, diagnose by looking at data');
    end
    end
end
if isfield(data,'saccEndPoint')
    data = rmfield(data,'saccEndPoint'); % either way, this gets removed
end

% SJ 07-2020
if isfield(data,'oneTargTrial')
    if isfield(data,'oneTargChoice')
        if length(data.oneTargTrial)<length(data.oneTargChoice) && length(data.oneTargChoice)==length(data.choice)
            data.oneTargChoice(1:length(data.oneTargTrial)) = data.oneTargTrial;
            data = rmfield(data,'oneTargTrial');
        else
            error('unsure, diagnose by looking at data');
        end
    elseif length(data.oneTargTrial)==length(data.choice)
        data.oneTargChoice = data.oneTargTrial;
        data = rmfield(data,'oneTargTrial');
    else
        error('unsure, diagnose by looking at data');
    end
end

% SJ 07-2020
if isfield(data,'oneConfTargTrial')
    if isfield(data,'oneTargConf')
        if length(data.oneConfTargTrial)<length(data.oneTargConf) && length(data.oneTargConf)==length(data.choice)
            data.oneTargConf(1:length(data.oneConfTargTrial)) = data.oneConfTargTrial;
            data = rmfield(data,'oneConfTargTrial');
        else
            error('unsure, diagnose by looking at data');
        end
    elseif length(data.oneConfTargTrial)==length(data.choice)
        data.oneTargConf = data.oneConfTargTrial;
        data = rmfield(data,'oneConfTargTrial');
    else
        error('unsure, diagnose by looking at data');
    end
end



disp('done.');

