% plot raster and psth plot for individual units


% 
%
% NOTE: you will need rh_plot function and cbrewer toolbox from FLutils
% repository to run this as is.
%
%
% Plot cosmetics need some work, happy to help with that (e.g. axes,
% titles, legend/keys). (and I will update the rh_plot on github)
%
% rh_plot is borrowed from my PhD lab and has some very unintelligible (academic) variable
% names, i've kept it around just for plotting purposes.
% FLutils also has a function called trial_psth which similarly calculates
% a PSTH for a single units spikes and has (much better) documentation on how I do this, 
% It's more designed for building population firing rates matrix, but it doesn't return rasters 
% (or any plots for that matter). Also probably needs testing again, and it
% doesn't have any in-built smoothing options.
%
%
% 
%

% EDIT variables in cells 1-4 to pass in data/change what is plotted.

%% 1. LOAD THE DATASTRUCT

subject = 'lucio';
date_range = 20231003:20231017;

datapath = '/Users/stevenjerjian/Desktop/FetschLab/Analysis/data/';
data_filename = sprintf('%s_%d-%d_neuralData.mat',subject,date_range(1),date_range(end));
load(fullfile(datapath,data_filename));

%% 2. Define the conditions, and how to plot them

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
hue_col = contains(condlabels, hue);

% and define a custom colormap for the hue
hue_colors = cbrewer('div','RdBu',length(hdgs));
hue_colors(hdgs==0,:) = [0 0 0]; % force zero hdg to be black

% set the text for plot titles
cond_titles = {'Ves', 'Vis', 'Comb'};


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
hue_col = contains(condlabels, hue);
hue_colors = 'bcg';


%}


%% 3. Define the alignment and range of PSTHs

expt = 'dots3DMP';

align_event = {'stimOn', 'saccOnset'};
other_events = {{'fpOn','fixation','stimOff','postTargHold'}, {'stimOn','postTargHold'}};
tmin = [-0.5, -0.5]; % seconds
tmax = [1.5, 1.5];


%% 4. pick a unit

u = 3;
sess = 2;

%% define some extra things

% assign for easier reference later
spiketimes = dataStruct(sess).data.(expt).units.spiketimes{u};
events     = dataStruct(sess).data.(expt).events;

% 'unique' conditions are what we want on the subplots, so use all the columns barring hue
[uconds,ia,ic1] = unique(conds(:,~hue_col),'rows'); 

% extract the trial conditions from this session
condsUnit = nan(length(events.trStart),length(condlabels));
for cond=1:length(condlabels)
    condsUnit(:,cond) = events.(condlabels{cond});
end

% and figure out how wide each sub-subplot should be
if length(align_event) > 1
    align_lens = tmax - tmin;
    sp_prc = align_lens / sum(align_lens);
end

%% let's plot



[p,~] = numSubplots(size(uconds,1));
f=figure('position',[100 100 1000 350*p(1)],'color','w');

for uc = 1:size(uconds, 1)

   ct = cond_titles{uc};
%     ct = sprintf('%s=%s', condlabels{1}, num2str(uconds(uc, 1)));
%     for c = 2:size(condsUnit(:,~hue_col), 2)
%         ct = [ct ', ' condlabels{c}, uconds(uc, c)];
%     end

    % get all the trials for this condition, that have a 'hue' within the
    % hue condition column 
    I = all(condsUnit(:,~hue_col)==uconds(uc,:),2) & events.goodtrial'; 
    I = I & ismember(condsUnit(:, hue_col), conds(:, hue_col));
    if sum(I)==0, continue, end

    cond_axis = subplot(p(1), p(2), uc);
    ca_pos = get(cond_axis, 'Position'); % for later

    % now loop over the alignment events
    for iae = 1:length(align_event)

        align       = events.(align_event{iae});
        eventsnames = [align_event{iae} other_events{iae}];

        % store the times of the 'other' events on the trials of interest
        otherevs = cell(1,length(other_events{iae}));
        for ioe=1:length(other_events{iae})
            otherevs{ioe} = events.(other_events{iae}{ioe});
        end

        clear oe
        for ioe=1:length(other_events{iae})
            oe{ioe} = otherevs{ioe}(I);
        end

        % determine the axis widths for each alignment by interval length
        hh = gca;
        if length(align_event) > 1
            if iae==1
                p1 = ca_pos(1);
            else
                p1 = p1 + ca_pos(3)*1.2*sp_prc(1:iae-1);
            end
            hh=subplot('position',[p1 ca_pos(2)+0.05 ca_pos(3)*sp_prc(iae)*0.95 ca_pos(4)-0.05]); %axis off;
        end
        hold on;

        [~,~,hhist(iae, uc),hrast] = rh_plot(spiketimes,align(I),'tmin',tmin(iae),'tmax',tmax(iae),...
            'histmethod','boxcar','bin',0.05,'hist_smbins',5,'plotsingletrialevents',1,...
            'split',condsUnit(I,hue_col),'splitcolors',hue_colors,'otherevents',oe,'eventsnames',eventsnames,...
            'title',ct);

        if iae>1
            hhist(iae, uc).YAxis.Visible = 'off';

        end


    end

end


