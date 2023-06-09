

% rh plots from unitStruct

clear;clc

subject   = 'lucio';
dateRange = 20220512:20230602;

dataPath = '/Users/stevenjerjian/Desktop/FetschLab/Analysis/data/lucio_neuro_datasets';
dataFileName = sprintf('%s_%d-%d_neuralData.mat',subject,dateRange(1),dateRange(end));
load(fullfile(dataPath,dataFileName));


%%

par = 'dots3DMP';
unitStruct = unitStruct_from_dataStruct(dataStruct, par);

mods = [1 2 3];
% cohs = [0.2 0.6];
cohs = [1 2]; % use inds instead, real-val cohs were sometimes different
deltas = 0;

% use actual headings vector, and then remove/ignore NaNs
% should get this vector from matching behavioral dataset eventually
% or generate behavioral dataset from neural dataStruct

switch par
    case 'dots3DMPtuning'
        hdgs = [-90 -60 -45 -25 -22.5 -12 0 12 22.5 25 45 60 90];
        hdgs([1 end]) = []; % drop +/- 90, weird platform motion
    case 'dots3DMP'
    hdgs = [-12 -6 -3 -1.5 0 1.5 3 6 12];
end

N    = length(hdgs);
splitcols = cbrewer('div','RdBu',N*2);
splitcols = splitcols([1:floor(N/2) end-floor(N/2):end],:);
splitcols(hdgs==0,:) = [0 0 0];

% generate the whole list of conditions we want to include
[hdg,modality,coh,~,~] = dots3DMP_create_trial_list(hdgs,mods,cohs,deltas,1,0); 
condsTuning = [modality,coh,hdg];
condlabels  = {'modality','coherenceInd','heading'}; % must be consistent order with conds lists below


alignEvent = {'stimOn'};
% otherEvents = {{'fpOn','fixation','stimOff'}};
otherEvents = {{'fpOn','fixation','stimOff','postTargHold'}};
tmin = [-1.5];
tmax = [3];


cohCol = contains(condlabels,'coherence');
hdgCol = contains(condlabels,'heading');

[uconds,ia,ic1] = unique(condsTuning(:,~hdgCol),'rows'); 
[p,~] = numSubplots(size(uconds,1));

condtitles = {'Ves','Vis-L','Vis-H','Comb-L','Comb-H'};

save_plots = 1;

%%

f=figure('position',[100 100 1000 700],'color','w');
set(f,'PaperUnits','centimeters');
set(f,'PaperPositionMode', 'manual');
x=f.Position(3); y=f.Position(4);
set(f,'PaperSize',[x,y],'PaperPosition',[0,0,x,y])


for u = 1:length(unitStruct)
    
    clf
    clear hhist ym

    events   = unitStruct(u).events;
    spktimes = unitStruct(u).spiketimes{:};

    condsUnit = nan(length(events.trStart),length(condlabels));
    for cond=1:length(condlabels)
        condsUnit(:,cond) = events.(condlabels{cond});
    end

    % force cohInd to match coh % do this in nsEvents?
    if length(unique(condsUnit(:,cohCol)))==1 && strcmp(condlabels{cohCol},'coherenceInd')
        condsUnit(:,cohCol) = double(events.coherence>=0.5)+1;
%         fprintf('Session %d: only 1 coh above 0.5, changing "coherence index" from 1 to 2 for consistency\n',ses)
    end

    
    unitcols = splitcols(ismember(hdgs,unique(condsUnit(:,hdgCol))),:);

    % this won't work for multiple alignEvents...

    for iae=1:length(alignEvent)
        otherevs = cell(1,length(otherEvents{iae}));

        align       = events.(alignEvent{iae});
        eventsnames = [alignEvent{iae} otherEvents{iae}];

        % modify here to break up e.g. stimOn_all
        for ioe=1:length(otherEvents{iae})
            otherevs{ioe} = events.(otherEvents{iae}{ioe});
        end

        nspk     = length(spktimes);

        if nspk==0, continue, end

        tt = sprintf('%s, Unit %d, ClusGrp=%d, ClusID=%d, nspk=%d', unitStruct(u).date, u, ...
            unitStruct(u).cluster_type, unitStruct(u).cluster_id, nspk);


        ym = nan(size(uconds,1),1);
        for uc = 1:size(uconds,1)

            subplot(p(1),p(2),uc);

            I = all(condsUnit(:,~hdgCol)==uconds(uc,:),2) & events.goodtrial'; % all trials from this condition vector
            I = I & ismember(condsUnit(:,hdgCol),hdgs);


            if sum(I)==0, continue, end

            clear oe
            for ioe=1:length(otherEvents{iae})
                oe{ioe} = otherevs{ioe}(I);
            end

            [~,~,hhist(uc),hrast] = rh_plot(spktimes,align(I),'tmin',tmin(iae),'tmax',tmax(iae),...
                'ratprc', 0.6, 'histmethod','boxcar','bin',0.05,'hist_smbins',5,'plotsingletrialevents',0,...
                'split',condsUnit(I,hdgCol),'splitcolors',unitcols,'otherevents',oe,'eventsnames',eventsnames,...
                'title',condtitles{uc});
                            
            ym(uc) = hhist(uc).YLim(2);

        end

        ht=suptitle(tt);
        ht.FontSize = 16;

        for s = 1:length(hhist)
            try
                hhist(s).YLim(2) = max(nanmax(ym), 10);
            catch
                disp('couldnt adjust y-axis height')
            end
        end

    end

    if save_plots

        subplot(p(1),p(2),p(1)*p(2))
        sc = 3; len = 3;

        hdgVec = hdgs;
        hdgXY = len .* [sind(hdgVec.*sc); cosd(hdgVec.*sc)];
        textVec = hdgXY .* [1.05; 1.1];
        startPoint = [0 0];
        % axis([-6 6 -6 6]);
        axis([-1 1 -1 1]*4); axis square
        hold on
        for h=1:length(hdgVec)
            plot(startPoint(1)+[0; hdgXY(1,h)],startPoint(2)+[0; hdgXY(2,h)],'linew',2,'color',splitcols(h,:));

            [th,r] = cart2pol(startPoint(1)+textVec(1,:),startPoint(2)+textVec(2,:));
            th = rad2deg(th);

            if hdgVec(h)>0, ha = 'left'; ra = th(h); va = 'middle';
            elseif hdgVec(h)<0, ha = 'right'; ra = th(h)-180; va = 'middle';
            else , ha = 'center'; ra = 0; va = 'bottom';
            end
            if hdgVec(h)==0 || abs(hdgVec(h))>2
                text(startPoint(1)+textVec(1,h),startPoint(2)+textVec(2,h),num2str(hdgVec(h)),'fontsize',14,'horizo',ha,'verti',va,'rotation',ra,'color',splitcols(h,:),'fontweight','bold')
            end
        end
        plot(startPoint(1)+[0; 0],startPoint(2)+[0; len],'k--')
        set(gca,'Visible','Off');

        printfig(f,sprintf('au_rh_%s_%s',par,datestr(now,'mm-dd-yyyy')),'ps',[],1);

    end

end