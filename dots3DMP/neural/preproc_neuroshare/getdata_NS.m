function [nsData] = getdata_NS(completeFilePath, dataType, dataChannel)

% dataType can be:
                    % Neural - 'Spike','LFP', 'Hi-Res', or 'Raw'; 
                    % Stim - 'Stim'
                    % Analog I/O - 'Analog 30k' or 'Analog 1k'
                    % Digital I/O - 'Digital'
                    
% dataChannel:
                   % A-1 starts at 1, A-2 starts at 33, 
                   % B-1 starts at 129
                   % Analog I/O starts at 10241
                   % Digital is 1 to 4 or 'parallel'
                   
                   
% Open the file and extract some basic information
[ns_status, hFile] = ns_OpenFileSJ(completeFilePath, 'single'); 

% Determine correct entityID for desired datastream
switch dataType    
    
    case {'Spike' 'LFP' 'Hi-Res' 'Raw' 'Analog 30k' 'Analog 1k'}
        EntityIndices = find([hFile.Entity(:).ElectrodeID] == dataChannel);         
        for i = 1:length(EntityIndices)       
            fileTypeNum = hFile.Entity(EntityIndices(i)).FileType;
            fileType = hFile.FileInfo(fileTypeNum).Type;            
            switch dataType
                case 'Spike' 
                    if strcmp('nev', fileType); entityID = EntityIndices(i); break; end 
                case {'LFP' 'Analog 1k'} 
                    if strcmp('ns2', fileType); entityID = EntityIndices(i); break; end 
                case 'Hi-Res'   
                    if strcmp('ns3', fileType); entityID = EntityIndices(i); break; end                
                case {'Raw' 'Analog 30k'} 
                    if strcmp('ns5', fileType); entityID = EntityIndices(i); break; end   
            end
        end

    case 'Stim'    
        entityID = find([hFile.Entity(:).ElectrodeID] == dataChannel + 5120);
        
    case 'Digital'
        if isnumeric(dataChannel)
            entityID = find(cellfun(@strcmpi, {hFile.Entity.Reason},...
                repmat({['Input Ch ' num2str(dataChannel)]},...
                size({hFile.Entity.Reason}))));
        else
            entityID = find(cellfun(@strcmpi, {hFile.Entity.Reason},...
                repmat({'Parallel Input'}, size({hFile.Entity.Reason}))));
        end        
end

% Extract channel info
[ns_RESULT, entityInfo] = ns_GetEntityInfo(hFile, entityID(end));

% Extract data
switch dataType            
    case 'Spike'
        for i = 1:entityInfo.ItemCount
            if i == 1; [ns_RESULT, nsSegmentSourceInfo] = ns_GetSegmentSourceInfo(hFile, entityID, i); end
            [ns_RESULT, spikeEventTime_s(i), spikeWindowData(:,i), sample_count(i), unit_id(i)] = ns_GetSegmentData(hFile, entityID, i);
        end
        % Sorted spikes
        uniqueUnits = unique(unit_id(unit_id~=0));
        SortedSpikeData = cell(1,1);
        for i = 1:length(uniqueUnits)
            SortedSpikeData{i} = spikeWindowData(:,unit_id==uniqueUnits(i));
        end
        nsData.spikeEventTime = spikeEventTime_s;
        nsData.spikeWindowData = spikeWindowData;
        nsData.sample_count = sample_count;
        nsData.unit_id = unit_id;
        
        nsData.SortedSpikeData = SortedSpikeData;

    case {'LFP' 'Hi-Res' 'Raw' 'Analog 30k' 'Analog 1k'}
        [ns_RESULT, analogInfo] = ns_GetAnalogInfo(hFile, entityID(end));     % analog info contains things like range and sampling rate
        
        TimeStamps = hFile.FileInfo(hFile.Entity(entityID).FileType).TimeStamps;
        numSamples = sum(TimeStamps(:,end));
        data = zeros(1,numSamples);
        startIndex = 1; 
        indexCount = TimeStamps(2,1);
        for i = 1:size(TimeStamps,2)                
            [~, ~, tempData] = ns_GetAnalogData(hFile, entityID, startIndex, indexCount);
            dataRange = TimeStamps(1,i) + (1:TimeStamps(2,i));
            analogInputData(dataRange) = tempData';
            clear tempData
            if i ~= size(TimeStamps,2) 
                startIndex = startIndex + TimeStamps(2,i);
                indexCount = TimeStamps(2,i+1);
            end
        end
        analogInputDataTime_s = (0:numSamples-1)' ./ analogInfo.SampleRate;
          
        nsData.numSamples = numSamples;
        nsData.analogData = analogInputData;
        nsData.analogTime = analogInputDataTime_s;
        nsData.analogInfo = analogInfo;
   
%     case 'Stim'   
%         % work through it backwards to setup variables first
%         for i = hFile.Entity(entityID).Count:-1:1      
%             [ns_RESULT, TimeStamp(i), Data(i,:), SampleCount] = ns_GetSegmentData(hFile, entityID, i);
%             if i == hFile.Entity(entityID).Count
%                 time = 0:1/3e4:TimeStamp(hFile.Entity(entityID).Count)+51/3e4;
%                 data = zeros(size(time));
%             end
%             Ts = find(time < TimeStamp(i),1,'last');
%             data(Ts+1:Ts+52) = Data(i,:);
%         end        
        
    case 'Digital'        
        % Get events and time stamps
        numCount = entityInfo.ItemCount;
        data = NaN(1, numCount); time = NaN(1, numCount); sz = NaN(1, numCount);
        for i = 1:numCount
            [~, time(i), data(i), sz(i)] = ns_GetEventData(hFile, entityID, i);
        end 
        
        nsData.data = data;
        nsData.time = time;
        nsData.sz = sz;
        
end

