% Fetsch Lab
% S. J. Jerjian
% Created May 2022
% Last updated August 2022
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
% **** A few notes for the user ****
% 
% Each unique row of dataStruct pertains to an individual 'set' of
% recordings, i.e. data recorded at a particular location. 
% % A 'good' recording day will likely only have one set, but 
% % in principle, there could be more than one set for a given day.
% Recordings are grouped by set, so multiple dots3DMP paradigms recorded in
% one set will be combined, whether the actual recording was done in
% multiple Trellis files or a single one - the "rec_group" group field in info struct is therefore critical!
% The timestamps for concatenated recordings are shifted according to the
% length of the overall data so that the range of events and spikes is
% matched for a given recording (nsEvents.analogData.timeStamps).
% i.e. if a recording is the first in the set, it's timestamps should start from around 0, whereas 
% % if it is later in the set they will start from some other time >0. The
% spiketimes will be matched accordingly, so the range of spiketimes should
% approximately match the range of e.g. nsEvents.Events.trStart.
%
% SJ 08-2022 added in metadata (getUnitInfo.m)
% SJ 08-2022 cleanUp option - to sub-select desirable units
% SJ 06-2022 significant updates
%            Fixed issues with processing of mksort data, 
%            + shifting of timestamps of multiple recordings. 
%            Switched dataCell from {} to dataStruct () struct format.
%
% Future improvements: 
% - include metadata about recorded units - depth, ch etc. (WORK IN
% PROGRESS 1) need to standardize for non-kilosort recordings too 2) need to store MDI and stopper depth)
% 
% - add option to load and append to existing dataStruct instead of
% generating from scratch at each run
% - consider whether any other info merits its own field in top level of
% dataStruct - MDI and stopper depth? distance to first contact, 'actual'
% distances between contacts
% could put this into the kilosort chanmap directly...instead of the
% unitless 1:32 currently used

% clear all
close all
addpath(genpath('/Users/stevenjerjian/Desktop/FetschLab/Analysis/codes/'))

%% decide which files to load

% which paradigms do we care about?
paradigms = {'dots3DMPtuning','dots3DMP','RFMapping','VesMapping'};

% at the end, will discard any units whose gross spike rate is <minRate spikes
% per second, or the number of trials of any of the unique stimulus

runCleanUp = 0;

% inputs to dots3DMP_NeuralStruct_runCleanUp
parSelect  = {'dots3DMPtuning','dots3DMP'}; 
minRate    = 5;
minTrs     = 5;

subject = 'lucio';

% dateRange = [20220223:20220331 20220512:20220531];
% dateRange = 20220805:20220902;
% dateRange = 20220901:20220902;
dateRange = 20220913;

dateStr = num2str(dateRange(1));
for d = 2:length(dateRange)
    dateStr = [dateStr newline num2str(dateRange(d))];
end

%%
keepMU = 1;           % include all SU and MU
useSCP = 1;
useVPN = 1;
overwriteLocalFiles = 0; % set to 1 to always use the server copy

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

%% exclude cells which were not adequately recorded in ALL fundamental experiments

if runCleanUp
    dots3DMP_NeuralStruct_runCleanUp(dataStruct,parSelect,minRate,minTrs);
end