% Dots_neuralAnalysis

% general script for extracting spike data and events, this can probably go
% back to preprocessing later
clear all
close all

subject = 'hanzo';
dateRange = 20201001:20201031;
% dateRange = 20201101:20201130;
% dateRange = 20201201:20201231;
% dateRange = 20210101:20210131;
% dateRange = 20210201:20210228;
% dateRange = 20210301:20210331;

% dateRange = 20201001:20210331;


folder = '/Users/chris/Documents/MATLAB/OpenEphys_data/';

% save currentFileList so you don't have to repeat the time consuming step
if sum(diff(dateRange)>1)==0
    file = [subject '_neuro_' num2str(dateRange(1)) '-' num2str(dateRange(end)) '.mat'];
elseif sum(diff(dateRange)>1)==1
    file = [subject '_neuro_' num2str(dateRange(1)) '-' num2str(dateRange(diff(dateRange)>1)) '+' num2str(dateRange(find(diff(dateRange)>1)+1)) '-' num2str(dateRange(end)) '.mat'];
else
    file = [subject '_neuro_' num2str(dateRange(1)) '---' num2str(dateRange(end)) '.mat'];
end


% temp
% file = 'fetsch_badProbe_continuous.dat';


load([folder file],'data');

data2 = struct;
data2.trialStart = [];
data2.motionStart = [];
data2.motionEnd = [];

%%
for m = 1:length(data.filename)
    
    slashes = strfind(data.filename{m},'/');
    currentFileName = data.filename{m}(slashes(end)+1:end);
    currentSubfolder = data.filename{m}(slashes(end-1)+1:slashes(end)-1);
    
    disp(['beginning folder /' currentSubfolder '...']);
    
    psortDataBase = Psort_read_psort(data.filename{m});               
    slash = strfind(psortDataBase.topLevel_data.file_fullPathOriginal,'/');
    if isempty(slash); slash = strfind(psortDataBase.topLevel_data.file_fullPathOriginal,'\'); end % Windows
    contFile = psortDataBase.topLevel_data.file_fullPathOriginal(slash(end)+1:end);
    contFullFile = fullfile(folder,subject,currentSubfolder,contFile);
    try
        disp('loading Continuous data...');
        [contData, contTS, contInfo] = load_open_ephys_data_faster(contFullFile);
        disp('done.');
    catch
        error(['file may be missing : ' contFullFile]);
    end
    
    data.spikeTimes{m} = contTS(logical(psortDataBase.topLevel_data.ss_index));

%     % sanity check: avg voltage at time of each spike should be
%     % a lot higher (or lower) than overall avg
%     figure;hist(contData,100);
%     title(num2str(median(contData)))
%     figure;hist(contData(logical(psortDataBase.topLevel_data.ss_index)),100);
%     title(num2str(median(contData(logical(psortDataBase.topLevel_data.ss_index)))));
    
    if data.run(m)==1
        evFile = 'messages.events';
    else
        evFile = ['messages_' num2str(data.run(m)) '.events'];
    end
    evFullFile = fullfile(folder,subject,currentSubfolder,evFile);
    
    try
        disp('loading Event data...');
        [evData, ~, ~] = load_open_ephys_data_faster(evFullFile);
%         [~, ~, ~, evData2] = load_open_ephys_data(evFullFile,1);
        disp('done.');
    catch
        error(['file may be missing : ' evFullFile]);
    end
            
    % events timestamps are, oddly enough, in units of 30,000ths of a sec;
    % div by 30k to convert to seconds, to align with continuous/sptimes
    data2.trialStart(end+1:end+length(evData.TrialStart)) = evData.TrialStart/30000;
    data2.motionStart(end+1:end+length(evData.MotionStart)) = evData.MotionStart/30000;
    data2.motionEnd(end+1:end+length(evData.MotionEnd)) = evData.MotionEnd/30000;
    
% %     if data.trialStart(1)>contTS(1) && data.trialStart(1)<contTS(end) && ...
% %        data.trialStart(end)>contTS(1) && data.trialStart(end)<contTS(end)
% %         data.valid(m,1) = 1;
% %     else
% %         data.valid(m,1) = 0;
% %         warning('Invalid timestamps (eg events data outside range of cont data) -- skipping this file');
% %         continue
% %     end
    
    
    % create spike rasters?
    

    keyboard
    
%     use filename to determine if mapping or expt, then
%     can also use this to figure out if it's MT or LIP mapping
    if ~all(isnan(evData2.Direction))
        keyboard
    end
    if ~all(isnan(evData2.TDirection))
        keyboard
    end


    % simple PSTH
    makePSTH_LIP_mapping
        
        
end
    







% TBD:
%
% createDataStructure
%
% 
% %% optional: save data struct to a mat file so you don't have to repeat the time consuming step
% file = [subject '_' num2str(dateRange(1)) '-' num2str(dateRange(end)) '.mat'];
% 
% data = rmfield(data,'dotPos'); % CAREFUL
% 
% save([localDir(1:length(localDir)-length(subject)-1) file], 'data');
% 
% % otherwise for larger files will need: 
% % save([localDir(1:length(localDir)-length(subject)-1) file], 'data','-v7.3');
% 
% 
% disp('done.');
% 
