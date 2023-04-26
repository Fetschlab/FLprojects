function [au, conds, nUnits] = dots3DMP_FRmatrix_fromDataStruct(dataStruct,par,timingInfo,aconds,condlabels,opts)
% au = DOTS3MP_FRMATRIX_FROMDATASTRUCT(dataStruct,par,timingInfo,conds,condlabels,opts)
%
% 
% conds and condlabels together allow extraction of trial-averaged firing rates (within an interval and/or binned around a certain event), for any
% condition grouping
%
% par       : 
% timingInfo: struct containing necessary fields for alignment and binning
%           - alignEvent    :
%           - otherEvents   :
%           - tStart        :
%           - tEnd          :
%           - binSize       : e.g. 0.05 (50ms)
%           - overlap       : e.g. 0.5
% conds      : c * n matrix, with one row per unique condition, and n columns to cover all the settings for one condition e.g. mod, coh, hdg
% condlabels : c x 1 cell array of strings, should match the order of conds, and refer to fields in dataStruct.events e.g. {'modality','coherence','heading'}
% opts       :


if ~isfield(opts,'keepMU'), opts.keepMU = 1; end
if ~isfield(opts,'smoothFR'), opts.smoothFR = 0; end
                
if opts.smoothFR==1 && ~isfield(opts,'convKernel')
    opts.convKernel = fspecial('average', [1 20]); % N bins wide boxcar (acausal, centered)
end

if ~isfield(opts,'keepSingleTrials'), opts.keepSingleTrials = 0; end
if ~isfield(timingInfo,'overlap'), timingInfo.overlap = 0; end
if ~isfield(opts, 'collapse_conds'), opts.collapse_conds = false(1, size(aconds,2)); end

% set defaults for timingInfo

alignEvent  = timingInfo.alignEvent;
otherEvents = timingInfo.otherEvents;
tStart      = timingInfo.tStart;
tEnd        = timingInfo.tEnd;
binSize     = timingInfo.binSize;
overlap     = timingInfo.overlap;

%% now loop over the dataStruct



clear au
clear condlist
unitInd = 0;
for ses = 1:length(dataStruct) 
    if ~isfield(dataStruct(ses).data,par) || ~isfield(dataStruct(ses).data.(par),'units')
        continue
    end
    temp = dataStruct(ses).data.(par);

    % kluge for now 11/06/22
    % stimOn sent from PLDAPS to Ripple is the 'motion' state begin, not
    % the motionOnsetLatency when actual motion begins...
    if ~isfield(temp.events, 'motionOn')
        temp.events.motionOn = temp.events.stimOn + 0.2994;
    end
    
    if opts.keepMU, unit_inds = find(temp.units.cluster_type<=3);
    else,           unit_inds = find(temp.units.cluster_type==2); % SU only
    end

    Ntr = length(temp.events.trStart);
    condlist = nan(Ntr,length(condlabels));

    % get the list of conditions on each trial for *this* session
    for c=1:length(condlabels)
        if strcmp(par,'VesMapping') && contains(condlabels{c}(1:7),'heading')
            condlist(:,c) = cell2mat(temp.events.(condlabels{c}));
        else
            condlist(:,c) = temp.events.(condlabels{c});
        end
    end


    % === this entire thing should be done in nsEventConditions really ===
    % (if we have heading as a condition), force very small ~eps headings to zero (some discrepancies)
    hdgCol = strcmp(condlabels,'heading');
    if any(hdgCol)
        condlist(abs(condlist(:,hdgCol))<0.01,hdgCol) = 0;
    end

    % force cohInd to match coh
    modCol = strcmp(condlabels,'modality');
    cohCol = contains(condlabels,'coherence');
    if length(unique(condlist(:,cohCol)))==1 && strcmp(condlabels{cohCol},'coherenceInd')
        condlist(condlist(:,modCol)>1,cohCol) = double(temp.events.coherence(condlist(:,modCol)>1)>=0.5)+1;
        %fprintf('Session %d: only 1 coh above 0.5, changing "coherence index" from 1 to 2 for consistency\n',ses)
    end

    % THEN force ves to be lowest coh in the set, kluge
    if any(cohCol) && any(modCol)
        condlist(condlist(:,modCol)==1,cohCol) = min(condlist(:,cohCol));
    end
    % =======
    
    % remove any unwanted trials - brfix, or condition not in the conds list
    isInCondList = false(1,Ntr);
    condI        = nan(1,Ntr); 
    for itr=1:Ntr
        [isInCondList(itr),condI(itr)] = ismember(condlist(itr,:),aconds,'rows');
    end
    goodtrs = temp.events.goodtrial & isInCondList;
    
    % condI stores the condition of each trial by its index in conds
    condI    = condI(goodtrs);
    condlist = condlist(goodtrs,:);

    % re-assign condI after collapsing across collapse_conds columns
    % we've already stored goodtrs here, so we should be ok
    [conds, ~, ic] = unique(aconds(:, ~opts.collapse_conds), 'rows', 'stable');
    condI = ic(condI);
%     conds = aconds;

    % startInd (+1) is the starting ind for units in the current session, unitInd will go from 1:length(unit_inds) for each alignment event
    startInd = unitInd;

    for iae=1:length(alignEvent)

        % for all of this, remember that alignment for time-res psth (and relative time of otherEvents) will
        % ALWAYS ALWAYS ALWAYS be to the first column in alignEvent
        % so, depending on what is passed in, set ae variable appropriately

        % if just one entry is given, align to this event in events
        if length(alignEvent{iae})==1
            ae = temp.events.(alignEvent{iae}{1})';
        else % if 2 events are given, then align to first one (tStart), and go up to tEnd from second one
             ae = [temp.events.(alignEvent{iae}{1}); temp.events.(alignEvent{iae}{2})]';
        end

        % if a third argument is given, re-assign ae as this proportion
        % away from the first event, towards the second event
        if length(alignEvent{iae})==3
            ae = ae(:,1) + (ae(:,2)-ae(:,1))*alignEvent{iae}{3};
        end
        ae = ae(goodtrs,:);

        % calculate times of 'otherEvents' relative to main alignment event
        oe_times = nan(size(ae,1),length(otherEvents{iae}));
        for ioe=1:length(otherEvents{iae})
            try
                oe_times(:,ioe) = temp.events.(otherEvents{iae}{ioe})(goodtrs)' - ae(:,1);
            catch
%                 fprintf('Could not calculate relative time of %s event\n',otherEvents{iae}{ioe});
            end
        end

        for ic = 1:size(conds,1)
            au.times.evTimes_bySession{iae}(ic,:,ses) = nanmedian(oe_times(condI==ic,:),1);
        end

        unitInd = startInd; % back to the 'top' (i.e. startInd), so that we are adding to the same unit index for each alignment event
        
        for u = 1:length(unit_inds)
            unitInd = unitInd+1; % ADD 1 to the INDEX HERE

%             fprintf('unit %d, ae %d\n',unitInd,iae)

            au.hdr.unitDate{unitInd} = dataStruct(ses).date;
            au.hdr.unitSet(unitInd)  = dataStruct(ses).rec_set;
            au.hdr.area{unitInd}     = dataStruct(ses).brain_area;
            au.hdr.unitID(unitInd)   = temp.units.cluster_id(unit_inds(u));
            au.hdr.unitType(unitInd) = temp.units.cluster_type(unit_inds(u));

            
            % calculate time-resolved and average firing rate within desired intervals
            [fr, x, fr_mean, durs, aeI]  = trial_psth(temp.units.spiketimes{unit_inds(u)},ae,'tStart',tStart(iae),'tEnd',tEnd(iae),'binSize',binSize,'overlap',overlap);

            if opts.smoothFR
                % extend time intervals a bit to avoid edge effects of conv
                t_ext = length(opts.convKernel)/2 * binSize;
                [fr, x1, ~, ~, ~]  = trial_psth(temp.units.spiketimes{unit_inds(u)},ae,'tStart',tStart(iae)-t_ext,'tEnd',tEnd(iae)+t_ext,'binSize',binSize,'overlap',overlap);

                [~,pos(1)] = min(abs(x1-x(1)));
                [~,pos(2)] = min(abs(x1-x(end)));

                fr = smoothRaster ( fr , opts.convKernel );
                fr = fr(:,pos(1):pos(2)); % cut back to size

            end

            % now get activity across trials of each condition in conds

            cond_ntrs   = zeros(size(conds,1),1);               % number of trials per cond       nconds x 1
            condFR      = nan(size(conds,1),size(fr,2));        % time-res rate during interval   nconds x nbins
            condFR_mean = nan(size(conds,1),size(fr_mean,2));   % mean rate during interval       nconds x 1
            condFR_sem  = condFR_mean;                          % sem rate                        nconds x 1
            durs_median = condFR_mean;                          % median

            trialFR      = cell(size(conds,1),1);  
            trialFR_mean = cell(size(conds,1),1); 

            au.times.evTimes_byUnit{iae}(:,:,unitInd) = nan(size(conds,1),length(otherEvents{iae}));

            for ic = 1:size(conds,1)
                
                if sum(condI==ic)

                    trialFR{ic}      = fr(condI==ic,:);
                    trialFR_mean{ic} = fr_mean(condI==ic);

                    % average for this condition
                    cond_ntrs(ic)    = sum(condI==ic);
                    
                    condFR(ic,:)     = nanmean(trialFR{ic},1);
                    %                 condFR_sem(ic,:) = nanstd(trialFR{ic},[],1)/sqrt(cond_ntrs(ic));

                    condFR_mean(ic)  = nanmean(trialFR_mean{ic},1);
                    condFR_sem(ic)   = nanstd(trialFR_mean{ic},[],1)/sqrt(cond_ntrs(ic));

                    % store other event times for each unit (for PSTH plot marking)
                    au.times.evTimes_byUnit{iae}(ic,:,unitInd) = nanmedian(oe_times(condI==ic,:),1);

                    durs_median(ic)  = nanmedian(durs(condI==c)); % and median duration
                end
            end


            

            % these could be slightly different lengths across units
            % so keep as cells and then we will 'matricise' them outside
            % the loop. maybe there is a more efficient way to do this
            temp_xvec{iae}{unitInd} = x;
            temp_psth{iae}{unitInd} = condFR;
            

            au.data.FRmean{iae}(:,unitInd) = condFR_mean;
            au.data.FRsem{iae}(:,unitInd) = condFR_sem;
            au.data.durs{iae}(:,unitInd)  = durs_median; 
            au.data.condntrs{iae}(:,unitInd) = cond_ntrs;

            if opts.keepSingleTrials
                au.data.FRtrial{iae}(:,unitInd) = trialFR;
                au.data.FRtrialmean{iae}(:,unitInd) = trialFR_mean;
            end
        end
    end
end

% finally, transfer the firing rate of each unit for each condition into a
% big matrix (conds x timepoints x units)
% set the start and end points using whichever xvec is the longest
% ... in practice, there might just be one or two bins difference on either
% side across sessions, and this is only really necessary if alignEvent
% involves two events with variable interval between them
% e.g. motionOn to RT+fixed_time

for iae=1:length(alignEvent)
    [maxlenx,pos] = max(cellfun(@length, temp_xvec{iae}));

    au.data.PSTHs{iae} = nan(size(temp_psth{iae}{1},1),maxlenx,unitInd);
    au.times.xvec{iae} = temp_xvec{iae}{pos};

    for u=1:unitInd
        this_psth = temp_psth{iae}{u};
        if aeI == 1
            st = 1;
            en = size(this_psth,2);
        elseif aeI == 2
            st = maxlenx - size(this_psth,2) + 1;
            en = maxlenx;
        end
        
        au.data.PSTHs{iae}(:,st:en,u) = this_psth;

    end
end

nUnits = unitInd;
    

