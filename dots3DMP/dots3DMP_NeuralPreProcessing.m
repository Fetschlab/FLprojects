% last updated June 2022
% SJ
%
% this creates a big struct containing neural activity and information from
% rigB (dots3DMP recordings)
% user should specify the subject, dateRange, and paradigms to draw from
%
% organization of dataCell is as follows:
%
% 1 - each row in dataCell pertains to one recording 'session' (block/set) (in practice, this means each day)
% 2 - within a row, 
%       info field, with header information, 
%       data field
% data field is further split into fields for each paradigm, with each
% paradigm field containing:
%   'events' condition information for each trial, as well as the timing of key stimulus and behavioral events)
%   'pldaps' useful pldaps information for each trial, largely bookkeeping
%   'cluster_id' 1 x N ids from Kilosort, where N = number of clusters
%   'cluster_type' 1 x N, 1 - MU, 2 - SU, 3 - noise (3 is skipped out, MU
%   will also be skipped if keepMU below == 0)
%   'spiketimes' 1 x N cell array containing spike times in seconds for
%   each cluster, within the range of the given paradigm
%
% Note that each unique row of dataCell pertains to an individual 'set' of
% recordings - there could be more than one for a given day.
% Recordings are grouped by set, so multiple dots3DMP paradigms recorded in
% one set will be combined, whether the actual recording was done in
% multiple Trellis files or a single one - rec_group in info is therefore
% critical!
% The timestamps for concatenated recordings are shifted according to the
% length of the overall data so that the range of events and spikes is
% matched for a given recording (nsEvents.analogData.timeStamps).
% i.e. if a recording is the first in the set, it's timestamps should start
% from near 0, if it is later in the set they will start from some other
% time.
%
%
% TODO 
% - add option to append to existing dataCell instead of creating new one
% from scratch? but need to be careful about pars list too
% - add more useful metadata to info?
% - key fields from info should get their own field in dataCell rows?
% - exclusion criteria for cells
%   - cluster_type, but also numSpikes/spikerate, presence in all
%   paradigms if selected for?
%
%
%
% SJ 06-2022 significant updates
%           Fixed issues with processing of mksort data, timestamps of
%           multiple recordings. 
%           Switched dataCell from {} to () struct format.

clear all
close all
addpath(genpath('/Users/stevenjerjian/Desktop/FetschLab/Analysis/codes/'))

%% decide which files to load

% which paradigms do we care about
paradigms = {'dots3DMPtuning','dots3DMP','RFMapping','VesMapping'};
% paradigms = {'RFMapping'};
% paradigms = {'dots3DMP'};

% currently, units will be included if they are recorded in *any* of
% paradigms, future version will provide the option to only include units
% that are in *all* paradigms selected for

subject = 'lucio';

% dateRange = [20220223:20220331 20220512:20220531];
dateRange = 20220801:20220802;

dateStr = num2str(dateRange(1));
for d = 2:length(dateRange)
    dateStr = [dateStr newline num2str(dateRange(d))];
end

%%
% set to 0 for testing, smaller size data, SU only
% set to 1 to keep all MUs
keepMU = 1; 
useSCP = 1;
useVPN = 0;
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

getNeuralEventsInfo; % grab the task events and info files
createSessionData;   % create the dataCell
