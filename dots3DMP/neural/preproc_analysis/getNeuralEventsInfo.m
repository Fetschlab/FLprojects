%% getNeuralEventsInfo;

% grabs the event timestamps and condition information, and info files for
% each recording session in the list

% requires password(s) to be entered if a auth key not available
% see e.g.:
% https://www.howtogeek.com/66776/how-to-remotely-copy-files-over-ssh-without-entering-your-password/

% 03-2022 SJ, separate to PLDAPS preprocessing, for spike sorted waveforms,
% timestamps, and RippleEvents


% TODO re-write this to match PLDAPS_preprocessing_nolocal - use striplines
% and no need to save files locally
 
if ~useVPN
    IPadd = 'fetschlab@172.30.3.33'; % MBI machine
else
    IPadd = 'fetschlab@10.161.240.133'; % probably off campus, try proxy IP (requires VPN)
end

% get file list from local dir
if ~exist(localDir,'dir')
    mkdir(localDir);
end
allFiles = dir(localDir);
if length(allFiles)>2
    localFileList = allFiles(3).name; % skip 1+2, they are are "." and ".."
    for n = 4:length(allFiles)
        localFileList = [localFileList newline allFiles(n).name];
    end
else
    localFileList = '';
end

% get folder list of recording sessions from remote dir
cmd = ['ssh ' IPadd ' ls ' remoteDir];
[~,remoteFolderList] = system(cmd, '-echo');
% [~,remoteFolderList] = system(cmd);
if any(strfind(remoteFolderList,'timed out')); error(remoteFolderList); end


% check each folder for a match to the desired date range, and download the RippleEvents and info files from the server
newlines = strfind(remoteFolderList,newline);
remoteFolders = cell(length(newlines),1);
m = 1; m2 = 1;
sess = 0;
for n = 1:length(newlines)
    if n==1
        remoteFolders{n} = remoteFolderList(1:newlines(1)-1);
    else
        remoteFolders{n} = remoteFolderList(newlines(n-1)+1:newlines(n)-1);
    end
    dateStart = strfind(remoteFolders{n},'20');
    if isempty(dateStart) || contains(remoteFolders{n},'Impedance'); continue; end
    
    thisDate = remoteFolders{n}(dateStart(1):dateStart(1)+7); % dateStart(1) just in case '20' appears in a timestamp
    
    if any(strfind(dateStr,thisDate))
        currentFolderList{m,1} = remoteFolders{n}; m=m+1;
        
        cmd = ['ssh ' IPadd ' ls ' remoteDir remoteFolders{n}];
        %         [~,thisFileList] = system(cmd, '-echo');
        [~,thisFileList] = system(cmd);
        
        newlines2 = strfind(thisFileList,newline);
        
        matFiles = cell(length(newlines2),1);
        for n2 = 1:length(newlines2)
            if n2==1
                matFiles{n2} = thisFileList(1:newlines2(1)-1);
            else
                matFiles{n2} = thisFileList(newlines2(n2-1)+1:newlines2(n2)-1);
            end
            [~,~,ext] = fileparts(matFiles{n2});
            if ~strcmp(ext,'.mat');  continue; end
            currentFileList{m2,1} = matFiles{n2}; m2=m2+1;
            
            if ~contains(localFileList,matFiles{n2}) || overwriteLocalFiles % always copy if overwrite option selected
                cmd = ['scp -r ' IPadd ':' remoteDir remoteFolders{n} '/' matFiles{n2} ' ' localDir];
                system(cmd,'-echo');
%                 system(cmd);
            else
                disp([remoteFolders{n} '/' matFiles{n2} ' exists locally, not copied']);
            end
            
        end
    end
end

