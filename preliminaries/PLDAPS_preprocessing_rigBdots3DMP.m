% grabs PLDAPS data from NAS as specified, then loads desired variables
% from .PDS files into a simpler data structure

% SJ 02-2023 
% this is a spin-off of the original versoin for dots3DMP task, given
% sufficient differences and additional parameters
% e.g.flags to addNexonarData and addDotPosition

% requires password(s) to be entered if a auth key not available 
% see e.g.:
% https://www.howtogeek.com/66776/how-to-remotely-copy-files-over-ssh-without-entering-your-password/


%% decide which files to load

paradigm = 'dots3DMP';
today    = str2double(datestr(now,'yyyymmdd'));

subject = 'lucio';
% dateRange = 20220512:today; % RT
% dateRange = 20220512:20230601; % RT
dateRange = 20220512:20220513;

% subject = 'zarya';
% dateRange = 20221201:20230925;

% subject = 'human';
% dateRange = 20190625:20191231; % non-RT
% % dateRange = 20200213:today; % RT

dateStr = num2str(dateRange(1));
for d = 2:length(dateRange)
    dateStr = [dateStr newline num2str(dateRange(d))];
end

localDir = ['/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/' subject '/']; % this should be set outside of version-controlled script
% localDir = ['/Users/chris/Documents/MATLAB/PLDAPS_data/' subject '/'];
remoteDir = ['/var/services/homes/fetschlab/data/' subject '/'];

% to load files directly from mounted homes Volume (be careful not to save over to the original file!)
mountDir  = ['/Volumes/homes/fetschlab/data/' subject '/'];

addNexonarDataToStruct = 0; % SJ 08-2021
addDotPositionToStruct = 0; % SJ 01-2022
addEyeMovementToStruct = 1; % SJ 10-2023
saveRewardData = 1;         % SJ 06-2022

%% get PDS files from server -- DON'T FORGET VPN
% will skip files that already exist locally, unless overwrite set to 1

useSCP = 0; % 1 - secure copy of files to local folder, 0 - load files directly from mounted drive, save locally only after cleanup
useVPN = 0; % 1 - use proxy VPN (off campus), 0 - use 172 address
overwriteLocalFiles = 1; % set to 1 to always use the server copy
getDataFromServer % now also includes pdsCleanup to reduce file size and complexity

%% get Nexonar files from server

if addNexonarDataToStruct
    localDir = ['/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/' subject '/nexonar/'];
    remoteDir = ['/var/services/homes/fetschlab/data/' subject '/' subject '_nexonar/'];
    
    getDataFromServer % get nexonar data, this will skip pdsCleanup
    
    % re-assign localDirs for createDataStructure
    localDir = ['/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/' subject '/'];
    localDirNex = ['/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/' subject '/nexonar/'];
end

%% create data structure

createDataStructure

%% optional: save data struct to a mat file so you don't have to repeat the time consuming step

if sum(diff(dateRange)>1)==0
    file = [subject '_' num2str(dateRange(1)) '-' num2str(dateRange(end)) '.mat'];
elseif sum(diff(dateRange)>1)==1
    file = [subject '_' num2str(dateRange(1)) '-' num2str(dateRange(diff(dateRange)>1)) '+' num2str(dateRange(find(diff(dateRange)>1)+1)) '-' num2str(dateRange(end)) '.mat'];
else
    warning('Multiple breaks in date range, could overwrite existing file with different intervening dates!');
    disp('Press a key to continue');
    pause
    file = [subject '_' num2str(dateRange(1)) '---' num2str(dateRange(end)) '.mat'];
end
    
try
data = rmfield(data,'dotPos'); % CAREFUL
catch
end

disp('saving...');
save([localDir(1:length(localDir)-length(subject)-1) file], 'data');

% otherwise for larger files will need: 
% save([localDir(1:length(localDir)-length(subject)-1) file], 'data','-v7.3');

disp('done.');




%%

% create csv file from struct

