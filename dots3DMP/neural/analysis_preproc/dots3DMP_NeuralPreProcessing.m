% Fetsch Lab
% S. J. Jerjian
% Created May 2022
% Last updated Feb 2023
%
% generate a struct containing spike sorted units and task data from rigB (dots3DMP recordings)
% user should specify the subject, dateRange, and paradigms
%
% organization of dataStruct is as follows:
%
% 1 - each row in dataStruct pertains to one recording 'session' (block/set) (in practice, this means each day)
% 2 - within a row, 
%       date and set
%       info field, with header information, 
%       data field
% data field is further split into fields for each paradigm recorded in that session. Each of these fields contains:
%
%   'events'    condition information for each trial, as well as the timing of key stimulus and behavioral events)
%   'pldaps'    useful pldaps information for each trial, largely bookkeeping
%   'units'     
%       'cluster_id'   :  1 x N ids, where N = number of clusters
%       'cluster_type' :  1 x N corresponding to cluster_ids above. By Kilosort convention,
%                               1 = MU, 2 = SU, (3 = noise, not included by default)
%                               if keepMU==1, MUs will be included, otherwise only SUs will be included in output struct
%       'spiketimes'   :  1 x N cell array containing spike times for each cluster, within the range of the given paradigm
%       'moreInfo'     :  1 x N cell array containing additional useful metadata      
% 
% **** Further notes for the user ****
% 
% 1. Each unique row of dataStruct pertains to an individual 'set' of
% recordings, i.e. data recorded at a particular location. 
% % A 'good' recording day with a full run of the task will likely only have one set, but 
% % in principle, there could be more than one set for a given day.
% 
% 2. Recordings are grouped by set, so multiple dots3DMP paradigms recorded in
% one set will be combined, whether the actual recording was done in
% multiple Trellis files or a single one - the "rec_group" group field in info struct is therefore critical!
% 
% 
% SJ 02-2023 select area for specific sessions
% SJ 01-2023 all recordings use phy format clusters
% SJ 08-2022 added in metadata (getUnitInfo.m)
% SJ 08-2022 cleanUp option - to sub-select desirable units
% SJ 06-2022 significant updates
%            Modified processing of mksort data, 
%            + shifting of timestamps of multiple recordings. 
%            Switched dataCell from {} to dataStruct () struct format.
%
%
% Feature list to implement: 
% - include metadata about recorded units (WORK IN PROGRESS)
% - add option to load and append to existing dataStruct instead of generating from scratch at each run
% - consider whether any other info merits its own field in top level of
% dataStruct
% - MDI and stopper depth? distance to first contact, 'actual'
% distances between contacts
%
% I should add to nsEventConditions to fix up the cohInd so that >0.5 coh
% is always '2' in the index? and ves is always low coh, or lowest within
% the set

%% decide which files to load

% which paradigms do we care about?
paradigms = {'dots3DMPtuning','dots3DMP','RFMapping','VesMapping'};

% at the end, will discard any units whose gross spike rate is <minRate spikes
% per second, or the number of trials of any of the unique stimulus

subject = 'lucio';
area    = 'PIVC';

dateRange = [20220727:20220804 20220913:20220923 20221003 20230220 20230221 20230223]; % PIVC

% this is not sustatinable in the case of dual recordings...i just need a
% way to mark units as coming from an area
% dateRange = dots3DMP_getSessionDates(subject, area);

dateStr = num2str(dateRange(1));
for d = 2:length(dateRange)
    dateStr = [dateStr newline num2str(dateRange(d))];
end

%%

keepMU = 1;           % include all SU and MU, by default, do it, can always remove them later
useSCP = 1;
useVPN = 0;
overwriteLocalFiles = 0; % set to 1 to always use the server copy
overwriteEventSets = 0;

% SJ 04-2022
% download associated PDS data files (we'll need this for some cleanup)
% actually, just specify the mountDir, don't bother downloading the files

% localDir = ['/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/' subject '/']; 
% remoteDir = ['/var/services/homes/fetschlab/data/' subject '/'];
% getDataFromServer;

mountDir  = ['/Volumes/homes/fetschlab/data/' subject '/'];    
PDSdir = mountDir; % reassign for later

localDir = ['/Users/stevenjerjian/Desktop/FetschLab/Analysis/data/' subject '_neuro/'];
remoteDir = ['/var/services/homes/fetschlab/data/' subject '/' subject '_neuro/'];
%     mountDir  = ['/Volumes/homes/fetschlab/data/' subject '/' subject '_neuro/'];

%% grab the task events and info files
% this time we will locally download these (_RippleEvents and info files)

getNeuralEventsInfo; 

%% create the dataStruct

createSessionData;   

%% create one .mat file for events from all paradigms for a given set
% this will be useful for pure behavior analyses (if desired)
% and importing into python

% paradigms = {'dots3DMPtuning','dots3DMP','RFMapping','VesMapping'};
% createSetEvents;


%% exclude cells which were not adequately recorded in ALL fundamental experiments

runCleanUp = 0; % don't do this here, do it later before analysis


if runCleanUp
    % inputs to dots3DMP_NeuralStruct_runCleanUp
    parSelect  = {'dots3DMPtuning','dots3DMP','RFmapping','VesMapping'};
    minRate    = 5;
    minTrs     = 5;

    dots3DMP_NeuralStruct_runCleanUp(dataStruct,parSelect,minRate,minTrs);
end