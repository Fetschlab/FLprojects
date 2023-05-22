 function processTrellisData(info,createBinaryFiles,createEvents,checkOverwrite)
% processTrellis Data
%
% extract behavioural events, and online detected/sorted spikes
% extract, and save raw neural recordings in binary format for kilosort
% sorting

% SJ updated 01/2022
% SJ 04/2022 now called by createBinaryFiles script, which can be run from command line
% SJ 04/2022 - tweaked to handle single Trellis recording for multiple
% PLDAPS (will work either way - nsEvents still created for each Trellis recording)
% SJ 06/2022 using sequential fwrite and memory-mapping to handle large GB
% ns5 files, as rawData could not all be loaded into the workspace before
% being written to int16 file (RAM issues). file is closed when end of the
% 'set' is reached. fseek used to reset the file to the beginnig at the
% beginning of a set, so previous bin file will be overwritten! 
% SJ 08/2022 - 10/2022 various bug fixes. Addition of flags for
% creatingBinaryFiles, checking overwrite of existing binary files, and
% creating nsEvents

% This is the main script to run after a recording for extracting and
% formatting synced task and neural data.

%   a) create a RippleEvents structure containing timings of all
%       key task events, trial conditions, and behavioral outcomes
%   b) add a .spkData field to this events structure containing any
%   threshold crossings and spikes detected online (through Trellis hoops)
%   c) save raw continuous data (.ns5 format) as int16 binary for offline
%   sorting in Kilosort  (one binary file per recording set!)
%   d) analogData field for re-aligned spikes or nsEvents offline
% creation of nsEvents struct requires createEventStruct_dots3DMP_NS, and
% nested getData_NS

% REQUIRES
% - neuroshare toolbox
%       - getdata_NS to extract raw data
%       - createEventStruct_NS to re-organize raw data events into trials

% clear;clc

if nargin<4, checkOverwrite=0; end
if nargin<3, createEvents=1; end
if nargin<2, createBinaryFiles=1; end

ADCtoUV = 0.25; % from Grapevine

%%
% files should be in date directory, and on NAS
nevFiles = dir([info.filepath '*.nev']);
if isempty(nevFiles)
    disp('No files found...remember to move files into dataFile folder and then to NAS!')
end

[uFiles,ia,ic] = unique(info.trellis_filenums);

% still create one nsEvents for each trellis file, but concatenate trellis files from same recording 'set'/group to create one binary file for  Kilosort
% presumably if trellis is recording over multiple PLDAPs files, trellis ids will each corresponding to unique 'set', but this is agnostic to that
% - so rec_group is maintained as a separate variable even within the same session, this allows flexibility to record PLDAPS protocols in separate trellis files or the same one


% SJ 12-2022 this works fine if each Trellis file belongs to only one rec_set. Otherwise no!


% matFiles={};
for f = 1:length(uFiles)
    % select all the valid PLDAPS files coming from the same trellis recording
    theseFiles = find(ic==f & ~isnan(info.pldaps_filetimes)');
    if isempty(theseFiles), continue, end

    clear pldaps_filenames
    for pf=1:length(theseFiles)
        pldaps_filenames{pf} = sprintf('%s%d%s%04d',info.subject,info.date,info.par{theseFiles(pf)},info.pldaps_filetimes(theseFiles(pf)));
    end
    nev_filename = sprintf('%s%ddots3DMP%.04d.nev',info.subject,info.date,uFiles(f));
    [~,name,~] = fileparts(nev_filename);

    % create and save nsEvents, one per .ns5 file
    if createEvents
        fprintf('\n processing task events for file %s ...\n',nev_filename)
        nsEvents = createEventStruct_dots3DMP_NS_multiPDS(nev_filename,info.par(theseFiles),pldaps_filenames,info.filepath); % newer 04/2022
    end

    if ~isempty(info.chanlist)

        for ch=1:length(info.chanlist)
            thisChan = info.chanlist(ch);

            % online sorted spikes from chans of interest
            if isfield(info,'chanInterest')
                for chI = 1:length(info.chanInterest{f})
                    if info.chanInterest{f}(chI)==info.chanlist(ch)
                        fprintf('extracting online spikes from ch%d\n',thisChan)
                        nsEvents.spkData.data{chI} = getdata_NS(nsEvents.hdr.Ripple_filepath, 'Spike', thisChan);
                        nsEvents.spkData.chs(chI)  = thisChan;
                    end
                end
            end
        end


        % identify trellis files corresponding to one recording set (this allows for flexibility in mapping between trellis files and 'sets')
        set_ind = find(info.rec_group==info.rec_group(ia(f)));
        thisFile = find(set_ind==ia(f));

        % if at the first file in the rec set, set some vars, make the dir if needed, and set a bof start position in the binary file
        if thisFile==1
            startInd=1; endInd=0;

            if createBinaryFiles
                bin_folder   = sprintf('%s%d_%d',info.subject,info.date,info.rec_group(ia(f)));
%                 bin_folder   = sprintf('%s%d_%d_ks25',info.subject,info.date,info.rec_group(ia(f)));
                bin_filename = sprintf('%s%d_%d',info.subject,info.date,info.rec_group(ia(f)));
                if ~exist([info.filepath bin_folder],'dir'), mkdir([info.filepath bin_folder]); end

                if checkOverwrite && exist([info.filepath '\' bin_folder '\' bin_filename '.bin'],'file')
                    answer = questdlg('Binary file for this dataset already exists...would you like to proceed (and overwrite)?','Yes','No');
                    if strcmp(answer,'No')
                        return
                    end
                end 

                % rewind the binary file
                fid = fopen([info.filepath '\' bin_folder '\' bin_filename '.bin'],'w');
%                 fid = fopen([info.filepath '\' bin_folder '\' bin_filename '_test.bin'],'w'); % SJ 10/06/2022 testing

                fseek(fid,0,'bof');
            end
        end

        % 04/15/2022 trying to read entire ns5 file into memory is resulting in errors, they are too big
        % 06/2022 writing sequentially to file is the only way...

        fprintf('extracting raw data from file %s, %d of %d in set\n',[name '.ns5'],thisFile,numel(unique(info.trellis_filenums(set_ind))))
  
        % attempts to read in whole ns5 file in one go, doomed to failure
%         temp = read_nsx(fullfile(info.filepath,[name '.ns5']),'keepint',true);
%         data = openNSxHL(fullfile(info.filepath,[name '.ns5']));
%         NSxToHL(fullfile(info.filepath,[name '.ns5']))

        % so just read the header info, and then memory-map the raw data.
        % pos tracks the start of the actual data in .ns5 file (i.e. after header info), so that we can offset the memory-mapping to this point

        % actually modified this version to not read in the data at all,
        % just get the pos. keep data as int16
        [temp,pos] = read_nsx_SJ(fullfile(info.filepath,[name '.ns5']),'keepint',true); 

        nsEvents.analogInfo = temp.hdr;
        nsEvents.analogInfo.timeStampsShifted = nsEvents.analogInfo.timeStamps + endInd; % in samples, will divide by 30000 to get to seconds
        endInd = startInd-1 + nsEvents.analogInfo.nSamples;

        if createBinaryFiles
            fprintf('Saving to binary file %s in %s\n',bin_filename,bin_folder)
            % read data sequentially (in nMin-length at a time, arbitrary but reasonable number)
            nMin     = 2; %5
            samps    = double(nsEvents.analogInfo.Fs)*60*nMin;

            % option 1: sequential fread? 08-11-2022
%             fid_ns5 = fopen(fullfile(info.filepath,[name '.ns5']),'r');
%             fseek(fid_ns5,pos,'bof'); % go to the beginning of the data
%             while ~feof(fid_ns5)
%                 data = fread(fid_ns5,[nsEvents.analogInfo.nChans,samps],'*int16'); % this is the right orientation! SJ 08-11-2022
%                 fwrite(fid,data,'*int16')
%             end

%           % option 2: sequential using read_nsx
%             begsample = 1;
%             endsample = 0;
%             samps = double(nsEvents.analogInfo.Fs)*3;
%             while endsample < nsEvents.analogInfo.nSamples
%                 endsample = min(samps+begsample-1, nsEvents.analogInfo.nSamples);
%                 
%                 data = read_nsx(fullfile(info.filepath,[name '.ns5']),'begsample',begsample,'endsample',endsample,'keepint',true,'readdata',true);
%                 fwrite(fid,data,'*int16')
% 
%                 begsample = endsample+1;
%             end
            

            % option 3: memmapfile
            % only for memory-mapped file, because m.data is vector
            dataMins = ceil(nsEvents.analogInfo.nSamples/samps); % number of nMin segments in data
            dpts     = double(samps*nsEvents.analogInfo.nChans); % number of datapoints per segment            
            mmf        = memmapfile(fullfile(info.filepath,[name '.ns5']),'Offset',pos,'Format','int16');

            for mm=1:dataMins   
                if mm==dataMins
                    theseSamples = (dpts*(mm-1)+1):length(mmf.Data);
                else
                    theseSamples = (1:dpts)+dpts*(mm-1);
                end
                fprintf('writing segment %d of %d (time %d-%d)...\n',mm,dataMins,theseSamples(1),theseSamples(end))
%                 data = ADCtoUV*mmf.Data(theseSamples);
                data = mmf.Data(theseSamples);

                % data should be nChans x nSamples for kilosort
                % fwrite writes in column order, so this concatenates
                % samples for each channel correctly!
                data = reshape(data,nsEvents.analogInfo.nChans,[]);
                fwrite(fid,data,'int16'); 
            end

            clear mmf 

            % if we've reached the end of a set, close the file
            if thisFile==numel(set_ind)
                fprintf('closing binary file %s\n',bin_filename)
                fclose(fid);
            end

             % for debugging only, to view some data (need to load a subset
             % of the data only)
%              fid = fopen([info.filepath '\' bin_folder '\' bin_filename '.bin'],'r');
%              X = fread(fid,'int16');
%              fclose(fid);
        
        
        else % don't create binary files, create files for Wave_Clus instead

%             wcf = sprintf('%s%d_wc',info.subject,info.date);
%             if ~exist(wcf,'dir')
%                 fprintf('%s directory does not exist yet...creating it',wcf)
%                 mkdir([info.filepath wcf]);
%             end
% 
%             for ch=1:temp.hdr.nChans
%                 data = temp.data(ch,:);
%                 matfilename = [name '_' num2str(temp.hdr.label{ch}) '.mat'];
%                 sr = temp.hdr.Fs; % default sr is 30kHz anyway, but set it to data here to be foolproof
%                 save(fullfile(info.filepath,wcf,matfilename),'data','sr') % the variable needs to be called 'data'!
%                 matFiles{f} = matfilename;
%             end
            
           
        end
        startInd = endInd+1;

    else
        disp('channel list is empty...no neural data extracted')
    end
 
    if createEvents
        fprintf('saving nsEvents...\n')
        save([info.filepath name '_RippleEvents.mat'],'nsEvents');
    end
end

% run Wave_Clus
%{
if ~createBinaryFiles
    param.stdmin = 3.5;
    param.detection = 'neg';
    cd(fullfile(info.filepath,wcf));
    Get_spikes(matFiles,'par',param);

    spikes_files = dir([info.filepath wcf '/*_spikes.mat']);
    spikes_files = {spikes_files.name};
    Do_cluster0ing(spikes_files);
 end
 %}


%{ 
%%----------------
% OLD, deprecated
% assumes that trellis files and pldaps files are a one-to-one match!

for f=1:length(nevFiles)

    % ASSUME THAT EVENT FILES/PLDAPS FILES ARE MATCHED FOR SIMULTANEOUS RECORDINGS
    if contains(nevFiles(f).name,info.subject) && ~isnan(info.pldaps_filetimes(f))

        %         dateStart = strfind(nevFiles(f).name,'20');
        %         thisDate  = nevFiles(f).name(dateStart(1):dateStart(1)+7);

        dot = strfind(nevFiles(f).name,'.nev');
        %         filenum = str2double(nevFiles(f).name(dot-4:dot-1));

        pldaps_filename = sprintf('%s%d%s%04d',info.subject,info.date,info.par{f},info.pldaps_filetimes(f));
        savename = [nevFiles(f).name(1:dot-1) '_RippleEvents.mat'];
        %         savename = [pldaps_filename '_RippleEvents.mat'];

        fprintf('\n processing task events for file %s [%s] (%d of %d)...\n',nevFiles(f).name,sprintf('%04d',info.pldaps_filetimes(f)),f,length(nevFiles))
        nsEvents = createEventStruct_dots3DMP_NS(nevFiles(f).name,info.par{f},pldaps_filename,info.filepath); % newer 02/2022
        %         nsEvents = createEventStruct_NS(subject,thisDate,paradigm,filenum,pldaps_filetime,filepath);

        

        if ~isempty(info.chanlist)
            bin_filename = sprintf('%s%d_%d',info.subject,info.date,info.rec_group(f));
            if ~exist([info.filepath bin_filename],'dir')
                mkdir([info.filepath bin_filename])
            end

            set_ind = find(info.rec_group==info.rec_group(f));
            thisFile = find(set_ind==f);

            if thisFile==1, rawData=int16([]); startInd=1; endInd=0; end

            [~,name,~] = fileparts(nevFiles(f).name);
            clear temp
            %                 completeFilePath = fullfile(info.filepath,nevFiles(f).name);

            for ch=1:length(info.chanlist)
                thisChan = info.chanlist(ch);

                % worth getting online sorted data for chanInterest?
                if isfield(info,'chanInterest')
                    for chI = 1:length(info.chanInterest{f})
                        if info.chanInterest{f}(chI)==info.chanlist(ch)
                            fprintf('extracting online spikes from ch%d\n',thisChan)
                            nsEvents.spkData.data{chI} = getdata_NS(nsEvents.hdr.Ripple_filepath, 'Spike', thisChan);
                            nsEvents.spkData.chs(chI)  = thisChan;
                        end
                    end
                end

                %                 % get raw continuous data from all channels
                %                 % this is the slowest part, esp for long recordings
                %                 fprintf('loading raw data from chan %d\n',thisChan);
                %                 temp = getdata_NS(fullfile(info.filepath,[name '.ns5']), 'Raw', thisChan);
                %
                %                 if ch==1 % these should be the same for all chs in a simultaneous recording...
                %                     nsEvents.analogData.numSamples = temp.numSamples;
                %                     nsEvents.analogData.sampleRate = temp.analogInfo.SampleRate;
                % %                     nsEvents.analogData.sampleTimeRel = temp.analogTime;
                % %                     nsEvents.analogData.sampleTimeAbs = temp.analogTime + endInd / nsEvents.analogData.sampleRate;
                %                     nsEvents.analogData.sampleTimeRel = temp.analogTime([1 end]);
                %                     nsEvents.analogData.sampleTimeAbs = temp.analogTime([1 end]) + endInd / nsEvents.analogData.sampleRate;
                %
                %                     % re-assign the end and startInds, and add the data!
                %                     endInd = startInd-1 + nsEvents.analogData.numSamples;
                %                     fprintf('Start Ind %d, End Ind %d\n',startInd, endInd);
                %                 end
                %
                %                 try
                %                     rawData(ch,startInd:endInd) = int16(temp.analogData);
                %                 catch
                %                     keyboard
                %                 end
                %
                %                 if ch==length(info.chanlist)
                %                     startInd = endInd + 1;
                %                 end
            end

            % 04/15/2022 use read_nsx tool to read all channels at once!
            % usually run this in command line with no java to circumvent
            % RAM issues with loading and appending .ns5 data

            fprintf('extracting raw data from file %s, %d of %d in set\n',[name '.ns5'],thisFile,numel(set_ind))
            temp = read_nsx(fullfile(info.filepath,[name '.ns5']),'keepint',true);
            nsEvents.analogInfo = temp.hdr;
            nsEvents.analogInfo.timeStampsShifted = nsEvents.analogInfo.timeStamps + endInd; % in samples, will divide by 30000 to get to seconds

            endInd = startInd-1 + nsEvents.analogInfo.nSamples;
            rawData(1:nsEvents.analogInfo.nChans,startInd:endInd) = int16(temp.data);

            startInd = endInd+1;

            % if we've reached the end of a set, write rawData to an int16 binary file for Kilosort
            if thisFile==numel(set_ind)
                fprintf('saving binary file %s\n',bin_filename)
                fid = fopen([info.filepath '\' bin_filename '\' bin_filename '.bin'],'w');
                fwrite(fid,rawData,'int16');
                fclose(fid);
            end

            % for debugging
            %             fid = fopen([info.filepath '\' bin_filename '\' bin_filename '.bin'],'r');
            %             X = fread(fid,'int16');
            %             fclose(fid);


        else
            disp('channel list is empty...no neural data extracted')
        end
        % save Events
        fprintf('saving nsEvents...\n')
        save([info.filepath savename],'nsEvents');

    end
end



%}




