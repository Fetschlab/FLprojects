% plot raster and psth plot for individual units

close all; clc

%% 1. LOAD THE DATASTRUCT

subject = 'lucio';
date_range = 20231003:20231017;

% datapath = 'C:\Users\yhaile2\Documents\CODE_Projects\3DMP_Data\Lucio\Neural/';
datapath = '/Users/stevenjerjian/Desktop/FetschLab/Analysis/data/';

data_filename = sprintf('%s_%d-%d_neuralData.mat',subject,date_range(1),date_range(end));
load(fullfile(datapath,data_filename));

%% 2. Define the conditions, and how to plot them

%{
% MY ORIGINAL VERSION
hdgs   = [-12 -6 -3 -1.5 0 1.5 3 6 12];
mods   = [1 2 3];
cohs   = 1; % coherence 'index' rather than true value
deltas = [0];

% generate the full list of conditions
[hdg,modality,coh,delta,~] = dots3DMP_create_trial_list(hdgs,mods,cohs,deltas,1,0); 
conds = [modality, coh, hdg, delta];
condlabels = {'modality','coherenceInd','heading', 'delta'};

% row, col (not yet implemented), would allow to specify row,col arrangement of subplots
% row = 'modality';
% col = 'coherenceInd';

% hue aka color - what condition variable do we want to plot across 
hue = 'heading';
hs_col = contains(condlabels, hue);

% and define a custom colormap for the hue
hue_colors = cbrewer('div','RdBu',length(hdgs));
hue_colors(hdgs==0,:) = [0 0 0]; % force zero hdg to be black

% set the text for plot titles
cond_titles = {'Ves', 'Vis', 'Comb'};
%}

%{
% this cell block is what needs work to plot different deltas
% you could pick just one heading e.g. 0, set deltas = [-3, 0, 3],
% and then select deltas as the 'hue'. then run all the subsequent cells
% the same.

% entires in condlabels (including hue) have to correspond to a column in conds
% and a field in events. so you would add an extra column/entry specifying your
preceding trial(s) condition in order to make separate subplots or hues.

% e.g.

hdgs = 0;
deltas = [-3 0 3];

% generate the full list of conditions
[hdg,modality,coh,delta,~] = dots3DMP_create_trial_list(hdgs,mods,cohs,deltas,1,0); 
conds = [modality, coh, hdg, delta];
condlabels = {'modality','coherenceInd','heading', 'delta'};

hue = 'delta';
hs_col = contains(condlabels, hue);
hue_colors = 'bcg';

%}

% we need to specify modality 3, otherwise ves and vis trials are included
% when vesox2==0 and delta==0
mods = 3; 
deltas = [-3 0 3];
vesox2 = [0 1];
hdgs = [-1 1];

% Use ndgrid to create a grid of all combinations
[deltaGrid, vesox2Grid, modGrid, hdgsGrid] = ndgrid(deltas, vesox2, mods, hdgs);
% Reshape the grids to get a list of unique combinations
conds = [deltaGrid(:), vesox2Grid(:), modGrid(:), hdgsGrid(:)]; 
condlabels = {'delta', 'VesOx2', 'modality', 'heading'};

hue = 'VesOx2';
style = 'heading';

hue_colors = 'kr'; % black for vesox2==0, red for vesox2==1
cond_titles = {'-delta', '0delta', '+delta'};


%% 3. Define the alignment and range of PSTHs

expt = 'dots3DMP';

align_event = {'stimOn', 'saccOnset'};
other_events = {{'fpOn','fixation','stimOff','postTargHold'}, {'stimOn','postTargHold'}};
tmin = [-0.5, -0.5]; % seconds
tmax = [1.5, 1.5];

binsize = 0.05; % 50ms bins

% apply 250ms window smoothing kernel (boxcar)
sm_method = 'boxcar';
sm_width = 0.25; % smooth


%% 4. pick a unit

u = 4;
sess = 1;

%% define some extra things

% assign for easier reference later
spiketimes = dataStruct(sess).data.(expt).units.spiketimes{u};
events     = dataStruct(sess).data.(expt).events;

% assign column for VesOx2
events.VesOx2 = PTE_genLI(events, {'modality&modality&modality'},{[0,1,2]},{[3,1,1]});
events.heading = sign(events.heading);

% 'unique' conditions are what we want on the subplots,
% so use condition columns that aren't hue or style
hs_cols = ismember(condlabels, {hue, style});
[uconds,~,~] = unique(conds(:, ~hs_cols), 'rows'); 

% extract the trial conditions from this session (Ntrs x length(condlabels))
condsUnit = nan(length(events.trStart),length(condlabels));
for cond=1:length(condlabels)
    condsUnit(:,cond) = events.(condlabels{cond});
end

% and figure out how wide each sub-subplot should be
% just a hack to make the sub-subplot widths prop. to the time duration shown
Nalign = length(align_event);
if Nalign > 1
    align_lens = tmax - tmin;
    sp_prc = align_lens / sum(align_lens);
end

%% let's plot

[p,~] = numSubplots(size(uconds,1));
f=figure('position',[100 100 1000 350*p(1)],'color','w');

for uc = 1:size(uconds, 1)

    % get all the trials for this condition
    I = all(condsUnit(:,~hs_cols)==uconds(uc,:),2) & events.goodtrial';

    % now we select only the trials for this condition that also have a
    % hue/style value matching the unique hues/styles in our grid
    % (this allows you to exclude certain hue/style conditions that are in the data,
    % if you want. e.g. if you set the hdgs sign to only be [-1, 1], 0 will get
    % dropped here)
    I = I & ismember(condsUnit(:, hs_cols), conds(:, hs_cols), 'rows');

    if sum(I)==0, continue, end

    cond_axis = subplot(p(1), p(2), uc);
    ca_pos = get(cond_axis, 'Position'); % for later

    % now loop over the alignment events
    for iae = 1:Nalign

        align       = events.(align_event{iae});
        eventsnames = [align_event{iae} other_events{iae}];

        % store the times of the 'other' events on the trials of interest
        otherevs = cell(1,length(other_events{iae}));
        for ioe=1:length(other_events{iae})
            otherevs{ioe} = events.(other_events{iae}{ioe});
        end

        % and make sure we also select the relevant trials for these
        clear oe
        for ioe=1:length(other_events{iae})
            oe{ioe} = otherevs{ioe}(I);
        end

        % determine the axis widths for each alignment by interval length
        hh = cond_axis;
        if Nalign > 1
            if iae==1
                p1 = ca_pos(1);
            else
                p1 = p1 + ca_pos(3)*1.3*sp_prc(1:iae-1);
            end
            hh=subplot('position',[p1 ca_pos(2)+0.05 ca_pos(3)*sp_prc(iae)*0.92 ca_pos(4)-0.05]); %axis off;
        end
        hold on;

        % rh_plot now uses current subplot hh, splitting it in two for raster and PSTH
        [~,~,hhist(iae, uc), hrast] = rh_plot(spiketimes, align(I), 'tmin',tmin(iae),'tmax',tmax(iae),...
            'bin', binsize , 'sm_method', sm_method, 'sm_width', sm_width, 'plotsingletrialevents',1,...
            'hue_style', condsUnit(I, hs_cols), 'hue_colors', hue_colors, 'otherevents',oe,'event_names',eventsnames);

        if iae > 1
            hhist(iae, uc).YAxis.Visible = 'off';

        else
            hrast.YLabel.String = 'trial #';
            hhist(iae, uc).YLabel.String = 'spikes s^{-1}';
            hrast.Title.String = cond_titles{uc};
        end


    end

end


