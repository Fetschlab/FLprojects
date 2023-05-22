% generate .mat file with recording info

clear info
createBinFile = 1; % set to 0 to skip (just create nsEvents)

addpath('C:\Program Files (x86)\Ripple\Trellis\Tools\neuroshare');
addpath(genpath('C:\Users\fetschlab\Documents\MATLAB\dots3DMP'));

%% run all
subject = 'lucio';
folder = ['\\172.30.3.33\homes\fetschlab\data\' subject '\' subject '_neuro\'];
files  = dir(folder);

for f = 13:length(files)
    clear info
    if strcmp(files(f).name(1:2),'20') && files(f).isdir
        info_filename = [subject files(f).name 'dots3DMP_info.mat'];

        try
            load(fullfile(folder,files(f).name,info_filename));
        catch
            continue
        end
        if length(info.chanlist)==32
            processTrellisData(info,1,0,0);
        end
    end
end

%% START ACTUAL OFFLINE PROCESSING HERE
% one cell per recording date (penetration)

%% 2023-01-31, pen 58

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20230131;
info.pen         = 58;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'30227_32ch_65um_80mm'};
info.gridxy              = {[2, 3]};
info.chanlist            = [1:32]; 

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1218 1231 1305 1542];
info.trellis_filenums    = [1 1 1 1];
info.par                 = {'RFmapping','dots3DMPtuning','dots3DMPtuning','dots3DMP'};
info.rec_group           = [1 1 1 1];
info.comments            = {'','.2 & .6 coh','nice behavior. high bets on high vis coh 0','6 reps per cond'};
info.chanInterest        = {[],[],[],[],[],[],[],[]};        % chans of interest for online thresholded spikes, not really using this anymore since can't reliably keep track of discrim quality with many chs

% >=1 value, different across probes (recording info)
info.depths              = {[10200 10200 10200 10200]}; % MDI depth
info.cellcomments        = {{'10,31,32 bad/broken','','6 to 7, then maybe gone...5 huge','5 got smaller'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,1);

%% 2023-01-27, pen 57

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20230127;
info.pen         = 57;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'10034_32ch_65um_80mm'};
info.gridxy              = {[3, 2]};
info.chanlist            = [1:32]; 

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1233 1245 1319 1353 1435 1530 1545 1555];
info.trellis_filenums    = [1 1 1 1 1 1 1 1];
info.par                 = {'RFmapping','dots3DMPtuning','dots3DMP','dots3DMP','dots3DMP','dots3DMP','dots3DMP','dots3DMP'};
info.rec_group           = [1 1 1 1 1 1 1 1];
info.comments            = {'','','','','','','',''};
info.chanInterest        = {[],[],[],[],[],[],[],[]};        % chans of interest for online thresholded spikes, not really using this anymore since can't reliably keep track of discrim quality with many chs

% >=1 value, different across probes (recording info)
info.depths              = {[10600 10600 10600 10600 10600 10600 10600 10600]}; % MDI depth
info.cellcomments        = {{'','','','','16 still v good, stable session','','',''}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,1);

%% 2023-01-20, pen 56

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20230120;
info.pen         = 56;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'10034_32ch_65um_80mm'};
info.gridxy              = {[3, 3]};
info.chanlist            = [1:32]; 

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1200 1213 NaN 1253 1323 1357];
                                       %1248 - ignore
info.trellis_filenums    = [1 1 1 1 2 2];
info.par                 = {'RFmapping','dots3DMPtuning','dots3DMP','dots3DMP','dots3DMPtuning','dots3DMP'};
info.rec_group           = [1 1 1 1 2 2];
info.comments            = {'','8 reps','ignore - eye tracking for PDW was off','started well, but cells changed, stopped early','8 reps','all mods, decent behavior'};
info.chanInterest        = {[],[14],[14],[],[11],[11]};        % chans of interest for online thresholded spikes, not really using this anymore since can't reliably keep track of discrim quality with many chs

% >=1 value, different across probes (recording info)
info.depths              = {[12000 12000 12000 12000 11800 11800]}; % MDI depth
info.cellcomments        = {{'','','ignore','lost cells, stopped early','good','some changes about halfway, e.g. 11 moved to 10, and then disappeared'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,1);

%% 2023-01-17, pen 55

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20230117;
info.pen         = 55;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'10034_32ch_65um_80mm'};
info.gridxy              = {[3, 1]};
info.chanlist            = [1:32]; 

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1457 1515 1540 1730];
info.trellis_filenums    = [1 1 1 1];
info.par                 = {'RFmapping','dots3DMPtuning','dots3DMP','RFmapping'};
info.rec_group           = [1 1 1 1];
info.comments            = {'a lot of brfix during stim, esp early','','R bias on vis, got a bit better. Low bet freq a bit low?',''};
info.chanInterest        = {[],[],[],[]};        % chans of interest for online thresholded spikes, not really using this anymore since can't reliably keep track of discrim quality with many chs

% >=1 value, different across probes (recording info)
info.depths              = {[11500 11500 11500 11500]}; % MDI depth
info.cellcomments        = {{'','13 getting bigger','13 continues getting bigger, then stable. two nice cells on 14. lots of activity 13-28',''}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,1);

%% 2022-12-15, pen 54

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20221215;
info.pen         = 54;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'10034_32ch_65um_80mm'};
info.gridxy              = {[3, 2]};
info.chanlist            = [1:32]; 

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1443 1456 1518 1530 1610 1700];
info.trellis_filenums    = [1 1 1 1 1 1];
info.par                 = {'RFmapping','dots3DMPtuning','dots3DMP','dots3DMP','dots3DMP','dots3DMP'};
info.rec_group           = [1 1 1 1 1 1];
info.comments            = {'6 thetas','','vis only','vis','comb, v good','ves. slight L bias'};
info.chanInterest        = {[],[],[],[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recor ding info)
info.depths              = {[10500 10500 10500 10500 10500 10500]}; % MDI depth
info.cellcomments        = {{'11 very nice','','','11 got a bit smaller, 2 moved to 1','','10/11 now'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,createBinFile);

%% 2022-12-09, pen 53

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20221209;
info.pen         = 53;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'10034_32ch_65um_80mm'};
info.gridxy              = {[3, 2]};
info.chanlist            = [1:32]; 

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1218 1224 1248 1344 1411 1433 1451];
info.trellis_filenums    = [1 1 1 1 1 1 2];
info.par                 = {'dots3DMPtuning','dots3DMPtuning','dots3DMP','dots3DMP','dots3DMP','dots3DMPtuning','RFmapping'};
info.rec_group           = [1 1 1 1 1 1 2];
info.comments            = {'forgot zero hdg','with zero hdg','comb only. good','forgot to turn off skewed reward for Z training...this could be why he was just betting high all the time','pretty much all high bets','ves+vis only',''};
info.chanInterest        = {[],[],[],[],[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[10800 10800 10800 10800 10800 10800 10800]}; % MDI depth
info.cellcomments        = {{'','weird crackling when platform moves','shift up e.g. 18 now on 17','','','','crackling again'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,createBinFile,0);

%% 2022-12-06, pen 52

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20221206;
info.pen         = 52;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'10034_32ch_65um_80mm'};
info.gridxy              = {[3, 3]};
info.chanlist            = [1:32]; 

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1431 1445 1541];
info.trellis_filenums    = [1 1 1];
info.par                 = {'RFmapping','dots3DMPtuning','dots3DMP'};
info.rec_group           = [1 1 1];
info.comments            = {'stopped early, lots of brfix today. Eye tracker issue?','','all 3 mods, poor on ves'};
info.chanInterest        = {[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[11150 11150 11150]}; % MDI depth
info.cellcomments        = {{'','',''}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,createBinFile);


%% 2022-10-03, pen 51

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20221003;
info.pen         = 51;
info.gridtype    = 'AP15_angled';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'10034_32ch_65um_80mm'};
info.gridxy              = {[3, 3]};
info.chanlist            = [1:32]; 

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1248 1302 1322 1333 1338 1357 1501];
info.trellis_filenums    = [1 1 1 1 1 1 1];
info.par                 = {'VesMapping','dots3DMPtuning','dots3DMP','dots3DMP','dots3DMP','dots3DMP','VesMapping'};
info.rec_group           = [1 1 1 1 1 1 1];
info.comments            = {'no FP','salvaged','ves only','short','','2+3, 0.3/0.7 coh','repeat'};
info.chanInterest        = {[],[],[],[],[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[12400 12400 12400 12400 12400 12400 12400]}; % MDI depth
info.cellcomments        = {{'','','','','','some drift 21-20, 18-17 etc',''}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,createBinFile);

%% 2022-09-23, pen 50

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220923;
info.pen         = 50;
info.gridtype    = 'AP15_angled';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'10034_32ch_65um_80mm'};
info.gridxy              = {[4, 1]};
info.chanlist            = [1:32]; 

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1509 1516 1537 1544 1623 1629 1654 1715 1738];
info.trellis_filenums    = [1 1 1 1 2 2 2 2 2];
info.par                 = {'VesMapping','dots3DMPtuning','dots3DMP','dots3DMP','VesMapping','dots3DMPtuning','dots3DMP','dots3DMP','VesMapping'};
info.rec_group           = [1 1 1 1 2 2 2 2 2];
info.comments            = {'','','','','','','','',''};
info.chanInterest        = {[],[],[],[],[],[],[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[12400 12400 12400 12400 10000 10000 10000 10000 10000]}; % MDI depth
info.cellcomments        = {{'','','','','','','','',''}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,createBinFile);

%% 2022-09-16, pen 48

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220916;
info.pen         = 48;
info.gridtype    = 'AP15_angled';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS

% >=1 value, fixed across session
info.probe_type          = {'Single electrode'};
info.probe_ID            = {'FHC'};
info.gridxy              = {[4, 0]};
info.chanlist            = [3]; 

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1644 1653 1714 1723 1741 1758];
info.trellis_filenums    = [1 2 2 3 4 5];
info.par                 = {'dots3DMPtuning','dots3DMP','dots3DMP','RFmapping','dots3DMPtuning','VesMapping'};
info.rec_group           = [1 1 1 1 1 1];
info.comments            = {'1+2 only, 0.9 coh','Ves only, started bad','short block','light was on in rig...','all 3 mods, 0.7coh only',''};
info.chanInterest        = {[3],[3],[3],[3],[3],[3]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[11050 11050 11050 11050 11050 11050]}; % MDI depth
info.cellcomments        = {{'huge v clean cell, decent 2nd one','','RF near center?','','cell started to get a bit smaller?','responding to up/down?'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,0);

%% 2022-09-13, pen 47

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220913;
info.pen         = 47;
info.gridtype    = 'AP15_angled';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'Single electrode'};
info.probe_ID            = {'FHC'};
info.gridxy              = {[4, 1]};
info.chanlist            = [3]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1551 1611];
info.trellis_filenums    = [1 2];
info.par                 = {'dots3DMPtuning','dots3DMP'};
info.rec_group           = [1 1];
info.comments            = {'all 3 mods, 0.9 coh','ves only, good'};
info.chanInterest        = {[3],[3]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[11100 11100]}; % MDI depth
info.cellcomments        = {{'2 cells? posisbly similar','more separate, ves responsive. lost cell ~150trials in'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,0);

%% 2022-09-02, pen 45

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220902;
info.pen         = 45;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'10034_32ch_65um_80mm'};
info.gridxy              = {[3, 1]};
info.chanlist            = [1:32]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1351 1407 1427 1507 1608];
info.trellis_filenums    = [1 1 1 1 1];
info.par                 = {'RFmapping','dots3DMPtuning','dots3DMP','dots3DMP','dots3DMPtuning'};
info.rec_group           = [1 1 1 1 1];
info.comments            = {'','all 3 mods','vis only, good','ves+comb. some R bias','repeat tuning'};
info.chanInterest        = {[],[],[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[8900 8900 8900 8900 8900]}; % MDI depth
info.cellcomments        = {{'13+28 nice','now 12/27','','2 units on ch3? now 11/26','11 smaller again'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
% processTrellisData(info,createBinFile);

%% 2022-09-01, pen 44

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220901;
info.pen         = 44;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'10034_32ch_65um_80mm'};
info.gridxy              = {[3, 3]};
info.chanlist            = [1:32]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1401 1423 1504 1547];
info.trellis_filenums    = [1 1 1 1];
info.par                 = {'dots3DMPtuning','dots3DMP','dots3DMP','dots3DMPTUNING'};
info.rec_group           = [1 1 1 1];
info.comments            = {'all 3 mods','vis only, good','ves+comb','repeat tuning'};
info.chanInterest        = {[],[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[10500 10500 10500 10500]}; % MDI depth
info.cellcomments        = {{'','','',''}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,createBinFile);
%% 2022-08-19, pen 43

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220819;
info.pen         = 43;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'10034_32ch_65um_80mm'};
info.gridxy              = {[3, 2]};
info.chanlist            = [1:32]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1330 1345 1407 1422 NaN 1432 NaN 1450 1558 1602];
info.trellis_filenums    = [1 1 1 1 1 1 1 1 1 1];
info.par                 = {'RFmapping','dots3DMPtuning','dots3DMP','dots3DMP','dots3DMP','dots3DMP','dots3DMP','dots3DMP','VesMapping','dots3DMPtuning'};
info.rec_group           = [1 1 1 1 1 1 1 1 1 1];
info.comments            = {'','brfix more than usual','ves only','hardly any low bets','ignore','more low bets','ignore','better. slightly low PDW on R','no fix req','good!'};
info.chanInterest        = {[],[],[],[],[],[],[],[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[10000 10000 10000 10000 10000 10000 10000 10000 10000 10000]}; % MDI depth
info.cellcomments        = {{'','','','','','','','','',''}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,createBinFile);

%% 2022-08-16, pen 42

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220816;
info.pen         = 42;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'10034_32ch_65um_80mm'};
info.gridxy              = {[2, 2]};
info.chanlist            = [1:32]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1239 1252 1305 1321 1413 NaN 1432 1437 1526 1541];
info.trellis_filenums    = [1 1 1 1 2 2 2 2 2 2];
info.par                 = {'RFmapping','dots3DMPtuning','dots3DMP','dots3DMP','dots3DMPtuning','dots3DMP','dots3DMP','dots3DMP','dots3DMP','RFmapping'};
info.rec_group           = [1 1 1 1 2 2 2 2 2 2];
info.comments            = {'','','ves only','vis+comb. Vis RTs faster than usual? - minRewardTime change?','','ignore - tried all 3 mods','also short','','good',''};
info.chanInterest        = {[],[],[],[],[],[],[],[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[10000 10000 10000 10000 11000 11000 11000 11000 11000 11000]}; % MDI depth
info.cellcomments        = {{'','','very clean cell on 4','stopped early, lots changed','','','','','10+32 seem to be responding','32 very big now, multiple responses to RFmapping'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,createBinFile);

%% 2022-08-15, pen 41

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220815;
info.pen         = 41;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'10034_32ch_65um_80mm'};
info.gridxy              = {[4, 1]};
info.chanlist            = [1:32]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1308 1331 1407 1422 1526 1543];
info.trellis_filenums    = [1 2 3 3 3 3];
info.par                 = {'dots3DMPtuning','dots3DMPtuning','dots3DMPtuning','dots3DMP','dots3DMP','RFmapping'};
info.rec_group           = [1 2 3 3 3 3];
info.comments            = {'','','','2+3, low coh vis PDW a bit flat','ves only, good',''};
info.chanInterest        = {[],[],[],[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[8850 9200 9000 9000 9000 9000]}; % MDI depth
info.cellcomments        = {{'','16+26 good at start, both got smaller','','','15 got nicer?',''}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,createBinFile);

%% 2022-08-12, pen 40

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220812;
info.pen         = 40;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'10034_32ch_65um_80mm'};
info.gridxy              = {[3, 2]};
info.chanlist            = [1:32]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1313 1335 1355 1521 1615];
info.trellis_filenums    = [1 1 1 1 1];
info.par                 = {'dots3DMPtuning','dots3DMP','dots3DMP','dots3DMP','RFmapping'};
info.rec_group           = [1 1 1 1 1];
info.comments            = {'all 3 mods','ves only','2+3, comb PDW very nice','delta! comb only. some L bias','lots of brfix at start'};
info.chanInterest        = {[],[],[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[10500 10500 10500 10500 10500]}; % MDI depth
info.cellcomments        = {{'14+18 nice','clearly ves responsive','','18 got a bit smaller','14+18 both clearly modulated!'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,createBinFile);
%% 2022-08-08, pen 39

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220808;
info.pen         = 39;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'10034_32ch_65um_80mm'};
info.gridxy              = {[3, 3]};
info.chanlist            = [1:32]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1332 1349 1411 1429];
info.trellis_filenums    = [1 1 1 1];
info.par                 = {'RFmapping','dots3DMPtuning','dots3DMP','dots3DMP'};
info.rec_group           = [1 1 1 1];
info.comments            = {'','all 3 mods','ves only','2+3, 2 cohs. Changing reward structure slightly from ves block to get low bets'};
info.chanInterest        = {[],[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[11250 11250 11250 11250]}; % MDI depth
info.cellcomments        = {{'14,28,11,29 ok. bad artefacts','','cells still stable. possibly ves?','noise artefacts several times during recording - 1:08, 1:20, 1:48, 1:58'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,createBinFile);

%% 2022-08-05, pen 38

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220805;
info.pen         = 38;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'10034_32ch_65um_80mm'};
info.gridxy              = {[3, 3]};
info.chanlist            = [1:32]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1408 1421 1442 1527 1541 1555 1603 1706];
info.trellis_filenums    = [1 1 1 2 2 2 2 2];
info.par                 = {'RFmapping','dots3DMPtuning','dots3DMP','dots3DMPtuning','dots3DMP','dots3DMP','dots3DMP','dots3DMPtuning'};
info.rec_group           = [1 1 1 2 2 2 2 2];
info.comments            = {'','0.7 coh only, all 3 mods','2+3 only, good PDW','1+2 only, 0.2+0.6 coh','ves only','2+3, more high bets','good PDW','more brfixes, but good'};
info.chanInterest        = {[],[],[],[],[],[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[10550 10550 10550 10740 10740 10740 10740 10740]}; % MDI depth
info.cellcomments        = {{'','','best cells got much smaller, stopped early. Some other good ones though','','lots of activity','','','several nice'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,createBinFile);

%% 2022-08-04, pen 37

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220804;
info.pen         = 37;
info.gridtype    = 'AP15_angled';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'Single electrode'};
info.probe_ID            = {'FHC'};
info.gridxy              = {[2, 2]};
info.chanlist            = [3]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1633];
info.trellis_filenums    = [1];
info.par                 = {'dots3DMPtuning'};
info.rec_group           = [1];
info.comments            = {'tried 8 reps per cond, but stopped early. all 3 mods'};
info.chanInterest        = {[3]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[12600]}; % MDI depth
info.cellcomments        = {{'small MU, got even smaller/less clear, could have been ves. SNR was poor today / unstable signal'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,0);
%% 2022-08-02, pen 36

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220802;
info.pen         = 36;
info.gridtype    = 'AP15_angled';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'Single electrode'};
info.probe_ID            = {'FHC'};
info.gridxy              = {[4,3]};
info.chanlist            = [3]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1300 1312];
info.trellis_filenums    = [1 2];
info.par                 = {'dots3DMPtuning','VesMapping'};
info.rec_group           = [1 1];
info.comments            = {'1+2 only','no fixation req. Ripple Events were not set to record!'};
info.chanInterest        = {[3],[3]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[7900 7900]}; % MDI depth
info.cellcomments        = {{'stable big cell'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,0);

%% 2022-08-01, pen 35

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220801;
info.pen         = 35;
info.gridtype    = 'AP15_angled';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'Single electrode'};
info.probe_ID            = {'FHC'};
info.gridxy              = {[3,1]};
info.chanlist            = [3]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1317];
info.trellis_filenums    = [1];
info.par                 = {'dots3DMPtuning'};
info.rec_group           = [1];
info.comments            = {'all 3 mods, 0.7 coh'};
info.chanInterest        = {[3],[3]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[11275]}; % MDI depth
info.cellcomments        = {{'at least 2 cells'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,0);

%% 2022-07-29, pen 34 - didn't record

%% 2022-07-27, pen 33

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220727;
info.pen         = 33;
info.gridtype    = 'AP15_angled';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'Single electrode'};
info.probe_ID            = {'FHC'};
info.gridxy              = {[2,3]};
info.chanlist            = [3]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1406 1530];
info.trellis_filenums    = [1 2];
info.par                 = {'dots3DMPtuning','dots3DMPtuning'};
info.rec_group           = [1 2];
info.comments            = {'1+2 only, 0.7 coh',''};
info.chanInterest        = {[3],[3]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[5450 13970]}; % MDI depth
info.cellcomments        = {{'good clean cell','sharp cell,got smaller'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,0);

%% 2022-07-20, pen 32

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220720;
info.pen         = 32;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'Single electrode'};
info.probe_ID            = {'FHC'};
info.gridxy              = {[5, 1]};
info.chanlist            = [3]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1048 1115 NaN 1125 1222 1250 1305];
info.trellis_filenums    = [1:7];
info.par                 = {'RFmapping','dots3DMPtuning','dots3DMPtuning','dots3DMPtuning','RFmapping','RFmapping','dots3DMPtuning'};
info.rec_group           = [1 1 1 1 2 3 3];
info.comments            = {'3 R, 3 ap','1+2, coh=1','ignore','coh=0.6','','2ap 3R',''};
info.chanInterest        = {[3],[3],[3],[3],[3],[3],[3]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[7600 7600 7600 7600 12550 12340 12340]}; % MDI depth
info.cellcomments        = {{'neautiful cell, very active','cell still good, wider MU also','','2 MU, small','cells better on the way up, probably the same cells as 05','very similar now, hopefully can separate!'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,0);

%% 2022-07-19, pen 30

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220719;
info.pen         = 31;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'Single electrode'};
info.probe_ID            = {'FHC'};
info.gridxy              = {[4, 3]};
info.chanlist            = [3]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1157 1201 1231 NaN 1300 1311 1338];
info.trellis_filenums    = [2:8];
info.par                 = {'dots3DMPtuning','dots3DMPtuning','RFmapping','RFmapping','RFmapping','dots3DMPtuning','RFmapping'};
info.rec_group           = [1 1 2 2 2 2 3];
info.comments            = {'stopped early','','finally realized problem with brfix - protVis was not being updated properly, this might mean its association with condition information is messed up from whenever the first brfix is inserted',...
    'room light was on for the first few trials','mouse stim pos was on, ignore','brfix a lot towards the end, finally killed the cell','2 apertures, 3R. finally fixed brfix issue after this session'};
info.chanInterest        = {[3],[3],[3],[3],[3],[3],[3]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[5300 5300 9535 9535 9535 9535 12260]}; % MDI depth
info.cellcomments        = {{'1 cell died, shifted','one cell still decent','good cell, smaller one underneath','same cell','cell still v good','finally killed the cell, near the end','large sharp cell + smaller MU'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,0);

%% 2022-07-15, pen 30

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220715;
info.pen         = 30;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'Single electrode'};
info.probe_ID            = {'FHC'};
info.gridxy              = {[3, 0]};
info.chanlist            = [3]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1252 NaN 1321 1339 1351 1417 1432 1505];
info.trellis_filenums    = [1 2 3 4 5 6 6 7];
info.par                 = {'RFmapping','RFmapping','RFmapping','dots3DMPtuning','dots3DMP','dots3DMPtuning','dots3DMPtuning','RFmapping'};
info.rec_group           = [1 2 2 2 2 2 2 3];
info.comments            = {'3 reps, salvaged again','ignore','3 apertures, 4 reps','coh=1','bad ves, short','1+2','3, same trellis file','salvaged'};
info.chanInterest        = {[3],[3],[3],[3],[3],[3],[3],[3]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[5400 6810 6810 6810 6810 6810 6810 11920]}; % MDI depth
info.cellcomments        = {{'wide unit, poss 2','','v clean cell','still good','','','cell finally got a bit smaller','cell became clean'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,0);


%% 2022-07-14, pen 29

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220714;
info.pen         = 29;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'Single electrode'};
info.probe_ID            = {'FHC'};
info.gridxy              = {[3, 2]};
info.chanlist            = [3]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1350];
info.trellis_filenums    = [1];
info.par                 = {'RFmapping'};
info.rec_group           = [1];
info.comments            = {'2 reps per cond only, salvaged PDS'};
info.chanInterest        = {[3]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[14200]}; % MDI depth
info.cellcomments        = {{'wide cell, got smaller towards the end'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,0);

%% 2022-06-15, pen 27

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220615;
info.pen         = 27;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'Single electrode'};
info.probe_ID            = {'FHC'};
info.gridxy              = {[3, 5]};
info.chanlist            = [3]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1502 1512 1528 1606 1626];
info.trellis_filenums    = [1 2 3 4 5];
info.par                 = {'dots3DMPtuning','dots3DMP','dots3DMP','dots3DMP','dots3DMPtuning'};
info.rec_group           = [1 1 1 2 2];
info.comments            = {'lots of brfix initially, coh 1 only','1 only, 20 per cond. initially lots of low bets','1+2+3, 0.7 coh, stopped early','L bias (ves)',''};
info.chanInterest        = {[3],[3],[3],[3],[3]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[9100 9100 9100 9120 9120]}; % MDI depth
info.cellcomments        = {{'very clean','consistent, probably a small MU too','','cell got a bit smaller but still decent',''}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,0);

%% 2022-05-26, pen 26

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220526;
info.pen         = 26;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'10034_32ch_65um_80mm'};
info.gridxy              = {[1, 0]};
info.chanlist            = [1:32]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1255 1310 1444 1452 1528 1541];
info.trellis_filenums    = [1 1 1 1 1 1];
info.par                 = {'dots3DMPtuning','dots3DMP','dots3DMP','dots3DMP','dots3DMPtuning','RFMapping'};
info.rec_group           = [1 1 1 1 1 1];
info.comments            = {'1+2 only','all 3 mods. some R bias on high coh. Forgot to set numReps so stopped early','stuck in zero coh corr loop, need to fix this','all high bet on vis...','',''};
info.chanInterest        = {[],[],[],[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[9400 9400 9400 9400 9400 9400]}; % MDI depth
info.cellcomments        = {{'','','','','',''}};
% this is effectively becoming chanInterest at this stage for probe recordings, can probably
% skip it unless there is something relevant about signal in general

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,1);

%% 2022-05-23, pen 25

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220523;
info.pen         = 25;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'10034_32ch_65um_80mm'};
info.gridxy              = {[1, 1]};
info.chanlist            = [1:32]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1214 1228 1249 1303 1325 1352 1502];
info.trellis_filenums    = [1 1 1 1 1 1 1];
info.par                 = {'RFMapping','dots3DMPtuning','dots3DMP','dots3DMP','dots3DMP','dots3DMP','dots3DMPtuning'};
info.rec_group           = [1 1 1 1 1 1 1];
info.comments            = {'just 4 thetas','all mods','all mods, ves PDW drifting low','ves only','vis + comb, more high bets','tried to encourage more low bets, but overall similar','all 3 mods'};
info.chanInterest        = {[],[],[],[],[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[10200,10200,10200,10200,10200,10200,10200]}; % MDI depth
info.cellcomments        = {{'','','','','','',''}};
% this is effectively becoming chanInterest at this stage, can probably
% skip it unless there is something relevant about signal in general

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,createBinFile);

%% 2022-05-20, pen 24

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220520;
info.pen         = 24;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'10034_32ch_65um_80mm'};
info.gridxy              = {[1, 1]};
info.chanlist            = [1:32]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1305 1316 1348 1411 1605];
info.trellis_filenums    = [1 2 2 2 2];
info.par                 = {'RFMapping','RFMapping','dots3DMPtuning','dots3DMP','dots3DMPtuning'};
info.rec_group           = [1 2 2 2 2];
info.comments            = {'ignore, testing mode','2 apertures, 2R, 8 theta','all 3 mods','all 3. more low bets on ves, especially at the start',''};
info.chanInterest        = {[],[],[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[9900 9600 9600 9600 9600]}; % MDI depth
info.cellcomments        = {{'','','','',''}};
% this is effectively becoming chanInterest at this stage, can probably
% skip it unless there is something relevant about signal in general

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,createBinFile);

%% 2022-05-17, pen 23

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220517;
info.pen         = 23;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'10034_32ch_65um_80mm'};
info.gridxy              = {[2, 5]};
info.chanlist            = [1:32]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1151 1158 1206 1237 1249 1335 1353 1402];
info.trellis_filenums    = [1 1 1 1 1 1 1 1];
info.par                 = {'VesMapping','dots3DMPtuning','dots3DMP','dots3DMP','dots3DMP','dots3DMP','dots3DMPtuning','VesMapping'};
info.rec_group           = [1 1 1 1 1 1 1 1];
info.comments            = {'no FP, horiz plane only','1+2,0.7 coh only','1 only','2+3, mostly betting high','changed rew/pens','lots more low bets, esp on vis','just 0.7 coh again','no FP'};
info.chanInterest        = {[],[],[],[],[],[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[9100 9100 9100 9100 9100 9100 9100 9100]}; % MDI depth
info.cellcomments        = {{'ch1,ch30 clear units, + others','ch31, weirdly regular firing','huge cell on ch10','ch7 got smaller, large sharp unit on ch21 + 32','2 cells on 5, huge on 16+24, 7 got smaller','','',''}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,createBinFile);

%% 2022-05-13, pen 22

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220513;
info.pen         = 22;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'10034_32ch_65um_80mm'};
info.gridxy              = {[3, 3]};
info.chanlist            = [1:32]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1152 1207 1220 1251 1350 1358 1406];
info.trellis_filenums    = [1 1 1 1 1 1 1];
info.par                 = {'VesMapping','dots3DMPtuning','dots3DMP','dots3DMP','dots3DMP','dots3DMP','dots3DMP'};
info.rec_group           = [1 1 1 1 1 1 1];
info.comments            = {'no fix req, but FP shown.','1+2','0.05 1-targ, Ves only. Good','Vis+Comb, stopped early as increasingly betting high','Still almost all PDW high','Mostly high bet','More low bets, changed incentives'};
info.chanInterest        = {[],[],[],[],[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[8900 8900 8900 8900 8900 8900 8900]}; % MDI depth
info.cellcomments        = {{'','','','','','',''}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,createBinFile);

%% 2022-05-12, pen 21

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220512;
info.pen         = 21;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'10034_32ch_65um_80mm'};
info.gridxy              = {[3, 3]};
info.chanlist            = [1:32]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1238 1252 1525];
info.trellis_filenums    = [1 1 1];
info.par                 = {'dots3DMPtuning','dots3DMP','RFMapping'};
info.rec_group           = [1 1 1];
info.comments            = {'','All 3 mods. brfix a lot, R biased on choices, PDW good','15 trials, brfix a lot, stopped early'};
info.chanInterest        = {[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[8100 8100 8100]}; % MDI depth
info.cellcomments        = {{'','lots of activity! Ref switch was the issue...',''}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,createBinFile);

% 2022-05-05, pen 20

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220505;
info.pen         = 20;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'10034_32ch_65um_80mm'};
info.gridxy              = {[2, 3]};
info.chanlist            = [1:32]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1441 1454 1517];
info.trellis_filenums    = [1 1 1];
info.par                 = {'dots3DMPtuning','dots3DMP','dots3DMP'};
info.rec_group           = [1 1 1];
info.comments            = {'no IMU rec. 0.2+0.7coh','2 only','1+3, PDW steeper',};
info.chanInterest        = {[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[10100 10100 10100]}; % MDI depth
info.cellcomments        = {{'lots of active chs but not clear SU it seems','activity decreased a bit? ch11 seems to be task-responsive',''}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,createBinFile);

%% 2022-04-08, pen 19

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220408;
info.pen         = 19;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'8714_32ch_65um_80mm'};
info.gridxy              = {[6, 3]};
info.chanlist            = [1:32]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1345 1420 1435 1506 1645];
info.trellis_filenums    = [1 2 3 4 5];
info.par                 = {'dots3DMPtuning','dots3DMPtuning','dots3DMP','dots3DMP','dots3DMPtuning'};
info.rec_group           = [1 2 2 2 2];
info.comments            = {'','','','','',''};
info.chanInterest        = {[],[],[],[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[6300 6500 6500 6500 6500]}; % MDI depth
info.cellcomments        = {{'','','','','',''}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,createBinFile);

%% 2022-04-07, pen 18

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220407;
info.pen         = 18;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'8714_32ch_65um_80mm'};
info.gridxy              = {[6, 3]};
info.chanlist            = [1:32]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1305 1325 1338 1402 1414 1445];
info.trellis_filenums    = [1 2 3 4 5 6];
info.par                 = {'dots3DMPtuning','dots3DMPtuning','dots3DMP','dots3DMPtuning','dots3DMP','dots3DMP'};
info.rec_group           = [1 1 1 2 2 2];
info.comments            = {'1+2 only','1+2 only','Vis only. Betting high, aborted','1+2 only','1 only, good!','2+3, good! 0 hdgs always rewarded right, BUG!!, 0.1 high target hidden'};
info.chanInterest        = {[],[],[],[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[7200 7200 7200 7700 7700 7700]}; % MDI depth
info.cellcomments        = {{'few units on mapped chs 12,21,27. Probably going for reward','','12 got smaller','12-31 activity','v clear on 12 and 21, 27 gets smaller','12,21,22,31 all good'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,createBinFile);

%% 2022-04-05, pen 17

clear info
% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220405;
info.pen         = 17;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'8714_32ch_65um_80mm'};
info.gridxy              = {[2, 3]};
info.chanlist            = [1:32]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1311 1325 1421];
info.trellis_filenums    = [1 2 3];
info.par                 = {'dots3DMPtuning','dots3DMP','dots3DMP'};
info.rec_group           = [1 1 1];
info.comments            = {{'1+2 only, 0.2+0.7 coh','Vis only,0.2+0.7,5% 1-target PDW, betting low ~10%. Slight R bias','Ves + Comb, 0.2+0.7. ~10% low bet again, Low and High Coh Comb look similar'}};
info.chanInterest        = {[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[9400 9400 9400]}; % MDI depth
info.cellcomments        = {{'ch30 possibly more than one unit','ch5 as well, ch30 maybe bigger until the end','ch30 smaller again, definitely task responsive'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,createBinFile);

%% 2022-04-01, pen 16
% 1 value, fixed across session
clear info

info.subject     = 'lucio';
info.date        = 20220401;
info.pen         = 16;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'8714_32ch_65um_80mm'};
info.gridxy              = {[4, 2]};
info.chanlist             = [1:32];

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1406];
info.trellis_filenums    = [1];
info.par                 = {'dots3DMPtuning'};
info.rec_group           = [1];
info.comments            = {{'1+2 only, 0.7 coh only'}};
info.chanInterest        = {[]};        % chans of interest for online thresholded spikes, can be different lengths

% >=1 value, different across probes (recording info)
info.depths              = {[7005]} ; % MDI depth
info.cellcomments        = {{'nothing really there'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,0);

%% 2022-03-14, pen 14
% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220324;
info.pen         = 14;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'Single electrode'};
info.probe_ID            = {'FHC'};
info.gridxy              = {[3, 0]};
info.chanlist             = 3;

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1432,1456];
info.trellis_filenums    = [1 2];
info.par                 = {'dots3DMPtuning','dots3DMPtuning'};
info.rec_group           = [1 2];
info.comments            = {{'all mods, 0.7 coh only','all mods, 0.7 coh only'}};
info.chanInterest        = {[3],[3]};        % chans of interest for online thresholded spikes, can be different lengths

% >=1 value, different across probes (recording info)
info.depths              = {[9929,10632]}; % MDI depth
info.cellcomments        = {{'very active, multiple cells. probably MU','also very active, probably MU'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
% save([info.filepath savefilename],'info','-mat')
processTrellisData(info,0); % no bin file for single ch recording

%% 2022-03-14, pen 11
% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220314;
info.pen         = 11;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'Single electrode'};
info.probe_ID            = {'AlphaOmega'};
info.gridxy              = {[2,-2]};
info.chanlist             = 3;

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1400,1528];
info.trellis_filenums    = [1 2];
info.par                 = {'RFMapping','dots3DMPtuning'};
info.rec_group           = [1 2];
info.comments            = {{'tried to sample space, 30+20apertureD','all mods, ~80 trials then forgot to stop Trellis recording!'}};
info.chanInterest        = {[3],[3]};        % chans of interest for online thresholded spikes, can be different lengths

% >=1 value, different across probes (recording info)
info.depths              = {[4966,7924]}; % MDI depth
info.cellcomments        = {{'very sharp, clear cell eventually, clear response to stim','probably MU,wider'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,0);

%% 2022-03-11, pen 10 
% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220311;
info.pen         = 10;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'Single electrode'};
info.probe_ID            = {'FHC'};
info.gridxy              = {[6,3]};
info.chanlist             = 3;

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1315,1522,1540,1551];
info.trellis_filenums    = [1 2 3 4];
info.par                 = {'dots3DMPtuning','dots3DMPtuning','dots3DMPtuning','dots3DMPtuning'};
info.rec_group           = [1 2 2 2];
info.comments            = {{'4trs per condition, 0.2 low coh','full set, started second PLDAPs at end, ignore last few trials (combined)','maybe ignore','started recording after PLDAPS'}};
info.chanInterest        = {[3],[3],[3],[3]};        % chans of interest for online thresholded spikes, can be different lengths for different recordings

% >=1 value, different across probes (recording info)
info.depths              = {[6076,11733,11733,11870]}; % MDI depth
info.cellcomments        = {{'good cell','','','cell better isolated, same as before?'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,0);

%% 2022-03-08, pen 9
% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220308;
info.pen         = 9;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'Single electrode'};
info.probe_ID            = {'FHC'};
info.gridxy              = {[4,4]};
info.chanlist             = 3;

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1441,1502];
info.trellis_filenums    = [1,2];
info.par                 = {'dots3DMPtuning','dots3DMP'};
info.rec_group           = [1 1];
info.comments            = {{'full protocol','short block'}};
info.chanInterest        = {[3],[3]};        % chans of interest for online thresholded spikes, can be different lengths

% >=1 value, different across probes (recording info)
info.depths              = {[12454,12454]}; % MDI depth
info.cellcomments        = {{'very clean cell','initially good, cell response changed a bit later'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,0);

%% 2022-03-07, pen 8
% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220307;
info.pen         = 8;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'Single electrode'};
info.probe_ID            = {'FHC'};
info.gridxy              = {[2,6]};
info.chanlist             = 3;

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1424];
info.trellis_filenums    = [1];
info.par                 = {'dots3DMPtuning'};
info.rec_group           = [1];
info.comments            = {{'90 trials'}};
info.chanInterest        = {[3]};        % chans of interest for online thresholded spikes, can be different lengths

% >=1 value, different across probes (recording info)
info.depths              = {[8591]}; % MDI depth
info.cellcomments        = {{'spike sorted online, ok. one mid-sized unit, and a larger one, possibly responding to reward'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,0);

%% 2022-03-01, pen 5
% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220301;
info.pen         = 5;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'Single electrode'};
info.probe_ID            = {'AlphaOmega'};
info.gridxy              = {[2,3]};
info.chanlist             = 1;

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1425, 1445, 1505, 1519, 1540, 1614];
info.trellis_filenums    = [1 2 3 4 5 7];
info.par                 = {'dots3DMPtuning','dots3DMPtuning','dots3DMPtuning','dots3DMP','dots3DMP','dots3DMPtuning'};
info.rec_group           = [1 2 3 3 3 4];
info.comments            = {{'','','','extra +1.5 hdg,vis+comb, mostly betting high','ves only, good',''}};  
info.chanInterest        = {1,1,1,1,1,1};        % chans of interest for online thresholded spikes, can be different lengths

% >=1 value, different across probes (recording info)  
info.depths              = {[5099, 5496, 6079, 6079, 6079, 7506]}; % MDI depth
info.cellcomments        = {{'not much, some cell bursts','one big cell low FR, got smaller, another smaller one','','cell very stable','got a bit smaller','huge cell'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
% save([info.filepath savefilename],'info','-mat')
processTrellisData(info,0);

%% 2022-02-23, pen 4

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220223;
info.pen         = 4;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type           = {'Single electrode'};
info.probe_ID            = {'AlphaOmega'};
info.gridtype             = 'odd_sym';
info.gridxy               = {[1,1]};
info.chanlist             = 1;

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1438, 1451, NaN, 1515, NaN, 1600];
info.trellis_filenums    = 1:6;
info.par                 = {'dots3DMPtuning','RFMapping','RFMapping','RFMapping','dots3DMPtuning','dots3DMPtuning'};
info.rec_group           = [1 1 2 3 4 5];
info.comments            = {{'','','ignore','forgot to stop the Trellis recording for a while','ignore',''}}; 
info.chanInterest        = {1,1,1,1,1,1};         % chans of interest for online thresholded spikes, can be different lengths

% >=1 value, different across probes (recording info)
info.depths              = {[7583, 7583, NaN, 7910, NaN, 8719]}; 
info.cellcomments        = {{'probably lost SU, then got another one?','','lost the cell','','',''}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,1);

%% 2022-02-18, pen 3

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220218;
info.pen         = 3;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type           = {'Single electrode'};
info.probe_ID             = {'FHC'};
info.gridxy               = {[0,0]};
info.chanlist             = 1;

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1622,NaN];
info.trellis_filenums    = [1 2];
info.par                 = {'dots3DMP','RFMapping'};
info.rec_group           = [1 2];
info.comments            = {{'R biased choices','pldaps recording wasn''t set up'}};
info.chanInterest        = {1,1};         % chans of interest for online thresholded spikes, can be different lengths

% >=1 value, different across probes (recording info)
info.depths              = {[9772,9521]}; % MDI depth
info.cellcomments        = {{'maybe an SU, bad online sorting, not super responsive','clearer unit'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
% save([info.filepath savefilename],'info','-mat')
processTrellisData(info,0);


%% TESTING , FAKE RECORDINGS
if 0

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220615;
info.pen         = NaN;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS

% >=1 value, fixed across session
info.probe_type          = {'NaN'};
info.probe_ID            = {'NaN'};
info.gridxy              = {[NaN,NaN]};
info.chanlist             = 1:32;

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [NaN];
info.trellis_filenums    = [1 1 1 2];
info.par                 = {'RFMapping','dots3DMPtuning','dots3DMP','dots3DMP'};
info.rec_group           = [1 1 1 2];
info.comments            = {{''}};
info.chanInterest        = {[],[]};        % chans of interest for online thresholded spikes, can be different lengths for different recordings

% >=1 value, different across probes (recording info)
info.depths              = {[NaN]}; % MDI depth
info.cellcomments        = {{'',''}};

% savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
% save([info.filepath savefilename],'info','-mat')
processTrellisData(info,0);
end