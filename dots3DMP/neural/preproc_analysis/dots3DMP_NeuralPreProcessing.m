% Fetsch Lab
% S. J. Jerjian
% Created May 2022
% Last updated May 2023
%
% This script generates a large dataStruct containing spike sorted units from rigB (dots3DMP recordings)
% 
% Included sessions (and the row structure of dataStruct) is specified in RecSessionInfo.xlsx
%
% organization of dataStruct:
%
% 1 - each row in dataStruct pertains to one recording 'session'
% (block/set) (in practice, this means each day) from a given probe
% each row contains relevant metadata about recording, and a 'data' field
% 2 - each data field contains subfields for each recorded paradigm
% paradigm can be one of:
% 'VesMapping','RFMapping','dots3DMPtuning','dots3DMP'
%
% each paradigm contains:
%
%   'events'    condition information for each trial, as well as the timing of key stimulus and behavioral events)
%   'pldaps'    useful pldaps information for each trial, largely bookkeeping
%   'units'     
%       'cluster_id'   :  1 x N ids, where N = number of clusters
%       'cluster_type' :  1 x N corresponding to cluster_ids above. By Kilosort convention,
%                               1 = MU, 2 = SU, (3 = noise, not included by default)
%                               if keepMU==1, MUs will be included, otherwise only SUs will be included in output struct
%       'cluster_label':  1 x N cluster labels (MU, SU, UN = 'unsorted')
%       'spiketimes'   :  1 x N cell array containing spike times for each cluster, within the range of the given paradigm
% 
% See createSessionData for understanding how the dataStruct is created 
% 
%
% Updates
%
% SJ 04-2023 extensive improvements to documentation
%            fixed issues with mismatched unitInfo - using readtable instead of textscan to read cluster_info files
%            restructured to use recording metadata spreadsheet for defining dataStruct 
% SJ 01-2023 all recordings use kilosort/phy (no more mksort - spikeInterface for single electrode recordings)
% SJ 08-2022 added in metadata (getUnitInfo.m)                          
%            cleanUp option - to sub-select desirable units
% SJ 06-2022 shifting of timestamps of multiple recordings. 
%            switched dataCell from {} to dataStruct () struct format.

clear;clc;close all

%% USER INPUTS

clear;clc

paradigms = {'dots3DMPtuning','dots3DMP','RFMapping','VesMapping'};

subject = 'lucio';
dateRange = 20220512:20230602;

keepMU = 1;           % include all SU and MU, by default, do it, can always remove them later
useSCP = 1;             % use scp to copy PDS files, or just load directly from mount
useVPN = 0;             % use off-campus proxy (vpn), or MBI address
overwriteLocalFiles = 0; % set to 1 to always use the server copy


%%

dateStr = num2str(dateRange(1));
for d = 2:length(dateRange)
    dateStr = [dateStr newline num2str(dateRange(d))];
end

% SJ 04-2022
% these will be user/OS specific...
localDir = ['/Users/stevenjerjian/Desktop/FetschLab/Analysis/data/' subject '_neuro/'];

PDSdir  = ['/Volumes/homes/fetschlab/data/' subject '/'];    
remoteDir = ['/var/services/homes/fetschlab/data/' subject '/' subject '_neuro/'];
mountDir  = ['/Volumes/homes/fetschlab/data/' subject '/' subject '_neuro/'];

%% grab the task events and info files
% this time we will locally download these (_RippleEvents and info files)
% TODO probably unnecessary, could just read them directly from mount
getNeuralEventsInfo; 

%% create the dataStruct

sess_info_file = '/Users/stevenjerjian/Desktop/FetschLab/Analysis/info/RecSessionInfo.xlsx';
createSessionData;   


%% exclude cells which were not adequately recorded in ALL fundamental experiments

% SJ don't do this in pre-processing, do within analysis script

% parSelect  = {'dots3DMPtuning','dots3DMP','RFmapping','VesMapping'};
% minRate    = 5; % min average firing rate
% minTrs     = 5; % unit must have spikes in at least 5 trials per condition
% 
% dots3DMP_NeuralStruct_runCleanUp(dataStruct,parSelect,minRate,minTrs);
