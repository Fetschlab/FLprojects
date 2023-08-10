 function processTrellisData(info,createBinaryFiles,createEvents,checkOverwrite)
% processTrellisData(info,createBinaryFiles,createEvents,checkOverwrite)
%
% WHAT THIS CODE DOES:
% This is the main workhorse after a recording for extracting and
% formatting synced task and neural data.
% extract behavioural events, and online detected/sorted spikes
% extract, and save raw neural recordings in binary format for kilosort
% sorting
%
% specifically,
%   a) create a RippleEvents structure containing timings of all
%       key task events, trial conditions, and behavioral outcomes
%   b) add a .spkData field to this events structure containing any
%   threshold crossings and spikes detected online (through Trellis hoops)
%   c) save raw continuous data (.ns5 format) as int16 binary for offline
%   sorting in Kilosort  (one binary file per recording set!)
%   d) add analogData field for re-aligned spikes or nsEvents offline
% creation of nsEvents struct requires createEventStruct_dots3DMP_NS, and
% nested getData_NS
%
% INPUTS:
% info file, generated from createTrellisInfo
% 
% Optional flags:
% createBinaryFiles - if false, skip binary file creation, default True
% createEvents - if false, skip nsEvents creation, default True
% checkOverwrite - if true, issue user warning if binary file already
% exists (since this function will overwrite), default False
%
% Having both createBinaryFiles and createEvents is probably redundant,
% because if both are set to false function does nothing...
%
% OUTPUTS:
% function does not return anything, but saves nsEvents and .bin files in
% relevant directories
%
%
% LOGS
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
% SJ 08/2023 - documentation updates
%
%
% DEPENDENCIES
% getdata_NS from Neuroshare toolbox to extract raw data
% read_nsx_SJ from nevutils-main for reading header information from binary
% files. Created a copy to not modify original.
% createEventStruct_NS to re-organize raw data events into trials

if nargin<4, checkOverwrite=0; end
if nargin<3, createEvents=1; end
if nargin<2, createBinaryFiles=1; end

% ADCtoUV = 0.25; % from Grapevine

%%
% files should be in date directory, and on NAS
nevFiles = dir([info.filepath '*.nev']);
if isempty(nevFiles)
    disp('No files found...remember to move files into dataFile folder, and then to NAS!')
end

[uFiles,ia,ic] = unique(info.trellis_filenums);

% still create one nsEvents for each trellis file, but concatenate trellis files from same recording 'set'/group to create one binary file for  Kilosort
% presumably if trellis is recording over multiple PLDAPs files, trellis ids will each corresponding to unique 'set', but this is agnostic to that
% - so rec_group is maintained as a separate variable even within the same session, this allows flexibility to record PLDAPS protocols in separate trellis files or the same one

% SJ 12-2022 this works fine if each Trellis file belongs to only one rec_set. Otherwise no!

% matFiles={};
for f = 1:length(uFiles)
    try

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

        % if at the first file in the rec set, set some vars, make the dir if needed
        % also set a bof start position in the binary file if 
        if thisFile==1
            startInd=1; endInd=0;

            if createBinaryFiles
                bin_folder   = sprintf('%s%d_%d',info.subject,info.date,info.rec_group(ia(f)));
                bin_filename = sprintf('%s%d_%d',info.subject,info.date,info.rec_group(ia(f)));
                if ~exist([info.filepath bin_folder],'dir'), mkdir([info.filepath bin_folder]); end

                if checkOverwrite && exist([info.filepath '\' bin_folder '\' bin_filename '.bin'],'file')
                    answer = questdlg('Binary file for this dataset already exists...would you like to proceed (and overwrite)?','Yes','No');
                    if strcmp(answer,'No')
                        warning("Aborting, no changes to binary file or nsEvents")
                        return
                    end
                end 

                % rewind the binary file
                fid = fopen([info.filepath '\' bin_folder '\' bin_filename '.bin'],'w');

                fseek(fid,0,'bof');
            end
        end

        fprintf('extracting raw data from file %s, %d of %d in set\n',[name '.ns5'],thisFile,numel(unique(info.trellis_filenums(set_ind))))
  
        % attempts to read in whole ns5 file in one go, various methods all
        % doomed to failure because files tend to be large
%         temp = read_nsx(fullfile(info.filepath,[name '.ns5']),'keepint',true);
%         data = openNSxHL(fullfile(info.filepath,[name '.ns5']));
%         NSxToHL(fullfile(info.filepath,[name '.ns5']))

        % Resolved here circa 06/2022 by just reading the header info, and then memory-mapping the raw data.

        % pos tracks the start of the actual data in .ns5 file (i.e. after
        % header info), so that we can offset the memory-mapping to this
        % point to create the binary file
        [temp,pos] = read_nsx_SJ(fullfile(info.filepath,[name '.ns5']),'keepint',true); 

        nsEvents.analogInfo = temp.hdr;
        nsEvents.analogInfo.timeStampsShifted = nsEvents.analogInfo.timeStamps + endInd; % in samples, will divide by 30000 to get to seconds
        endInd = startInd-1 + nsEvents.analogInfo.nSamples;

        if createBinaryFiles
            fprintf('Saving to binary file %s in %s\n',bin_filename,bin_folder)
            % read data sequentially (in nMin-length at a time, arbitrary
            % but reasonable number) Could eventually chunk based on
            % recording length, as kilosort does
            nMin     = 2; %5
            samps    = double(nsEvents.analogInfo.Fs)*60*nMin;

            % VARIOUS OPTIONS for reading data in chunks, all should do the
            % same thing although testing was not 100% foolproof

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
            % for memory-mapped file, because m.data is vector
            dataMins = ceil(nsEvents.analogInfo.nSamples/samps); % number of nMin segments in data
            dpts     = double(samps*nsEvents.analogInfo.nChans); % number of datapoints per segment            
            mmf      = memmapfile(fullfile(info.filepath,[name '.ns5']),'Offset',pos,'Format','int16');

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
                % fwrite writes in column order, so this concatenates samples for each channel correctly!
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
        
        
        else 
            
            % don't create binary files, create files for Wave_Clus instead
            % brief attempt to use wave_clus for single-channel sorting
            % eventually switched to using SpikeInterface, since we can use
            % a docker image of wave_clus without all this faff to create
            % local files, so this is deprecated. can be removed in future
            % commit.

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
%                 save(fullfile(info.filepath,wcf,matfilename),'data','sr') % the variable needs to be called 'data' for wave_clus to be happy!
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

    catch
        fprintf('something wrong with file %d',uFiles(f))
    end
end