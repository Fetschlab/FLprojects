% load/plot data from Yared struct

clear
close all
load('/Users/chris/Downloads/lucio_20230331-20230701_neuralData_clean.mat');
dbstop if error


%%
plotTuningFits = 0;
plotTuningPSTH = 0;
plotTaskPSTH = 1;

excludes = []; % to exclude non-selective cells, etc.

extRaster = [50 100]; % extend raster this much before stim onset (1) and after stim offset (2)
latency = [100 -100]; % boundaries of spike rate window relative to stim onset (1) and offset (2); latency(2) must be < extRaster(2)!
convKernel = fspecial('average', [1 40]); % N ms wide boxcar (acausal, centered)

% some manual plotting lims
psth_xlim1 = 1500; % endpoint for align stim on
psth_xlim2 = [-500 100]; % start and endpoint for align RT

maxHdg = 1.5; % maximum heading (absolute value) for inclusion in CP+confP and choice/conf-cond PSTH
% should be 0 or close to 0! unless testing code / sanity checks, or calcuating residuals

minRT = 400; % short RT can muck things up in several ways

dots3DMP_neuralAnalysis_prelims



%%
%%%%%%%%%%%%%%%%
% BEGIN MAIN ANALYSIS LOOP: calc rasters, rates/counts, and psths for each unit
%%%%%%%%%%%%%%%%
C=1; % counter for valid units, across all sessions

for n = 1 %:length(newDataStruct)
    disp(['session ' num2str(n)]);
    
    for c = 1:nU(n)
        
        disp(['cell ' num2str(c)]);
        disp(['total cell count ' num2str(C)]);

        if ~isempty(excludes)
            if ismember([n c],excludes,'rows')
                continue
            end
        end
        
        disp(num2str([n c]));
        
    %*****************************
    % TUNING
    %*****************************
        ev = newDataStruct(n).data.dots3DMPtuning.events;
        un = newDataStruct(n).data.dots3DMPtuning.units;
        good = logical(ev.goodtrial);

        mStart = ev.motionOn(good)';
        mEnd = ev.stimOff(good)';
        dur = round((mEnd-mStart)*1000);
        mod = ev.modality(good)';
        mods = unique(mod);
        coh = ev.coherence(good)';
        cohs = unique(coh);
            if length(cohs)>3; keyboard; end
            if length(cohs)==3
%                 warning('three cohs?');
                coh(coh==cohs(2)) = cohs(3); % pool largest two cohs, at least for tuning
            end
        cohs = unique(coh);
        hdg = ev.heading(good)';
        hdgs = unique(hdg);
        
        spikeTimes = un.spiketimes{c};
        startTrial = mStart-extRaster(1)/1000;
        endTrial = mEnd+extRaster(2)/1000;
        trialLen = endTrial - startTrial; % in seconds
        
        nTr = length(mStart);
        spRate_tun = nan(nTr,1);
        
        for t = 1:nTr % now includes only goodtrials                
            if round(trialLen(t)*1000)<size(raster_tun,3) 
                raster_tun(C,t,1:round(trialLen(t)*1000)) = 0; %#ok<SAGROW> % init bins to zero but keep nans after extRaster(2) cutoff
                spkInd = spikeTimes >= startTrial(t) & spikeTimes <= endTrial(t);
                spkTimes_thisTr = spikeTimes(spkInd) - startTrial(t);                
                if isempty(spkTimes_thisTr)
                    spRate_tun(t) = 0;
                else
                    x = 0:0.001:ceil(spkTimes_thisTr(end)*1000)/1000; % bin spikes only up to the last one; after that will be zeros, followed by NaNs, with a cutoff determined above (see comment)
                    raster_tun(C,t,1:length(x)-1) = histcounts(spkTimes_thisTr,x); %#ok<SAGROW>
                    a = latency(1)+extRaster(1);
                    b = a + dur(t) - latency(1) + latency(2);
                    spRate_tun(t) = sum(raster_tun(C,t,a:b)) / ((b-a+1)/1000);
                end
            end % else, trial exceeds the max length specified in prealloc, likely invalid, loss of eye signal etc (will stay nans)
            
%             a = latency(1)+extRaster(1);
%             b = a + dur(t) - latency(1) + latency(2);
        end

        if sum(isnan(spRate_tun(length(dir)-5:length(dir)))) > 3; keyboard; end % possible mismatch in ntrials given by dir, vs. spRate

        % plot tuning curves & PSTHs from tuning block

        clr{1,1} = 'k';
        clr{2,1} = 'm'; clr{2,2} = 'r';
        clr{3,1} = 'c'; clr{3,2} = 'b';
        clr2 = {'g','r'};
        mlabel{1} = 'ves'; mlabel{2} = 'vis'; mlabel{3} = 'comb'; 
        for m = 1:length(mods)
%             for o = 1:length(cohs) % keep it simple, show only highest coh
                o = 1;
                if m==1
                    I = mod==1;
                    if o==2; continue; end
                else
%                     I = mod==mods(m) & coh==cohs(o);
                    I = mod==mods(m) & coh==max(cohs); % keep it simple, show only highest coh
                end
                [FRmean, FRse, ~] = calc_mean(spRate_tun(I),hdg(I),hdgs);
                peak = hdgs(FRmean==max(FRmean)); peak = peak(1);
                    % vonMises params: b = [ampl, kappa, theta, offset]
                guess = [max(FRmean)-min(FRmean), 1.5, peak, min(FRmean)];
                [beta, ~] = fminsearch(@(x) tuning_vonMises_err(x,hdg(I),spRate_tun(I)), guess);
                prefhdg_fit(n,c,m,o) = beta(3);
                    % significance of tuning by anova
                pvalTuning(n,c,m,o) = anovan(spRate_tun(I),hdg(I),'display','off');

                if plotTuningFits
                    figure(n*10+c); set(gcf, 'Color', 'w', 'Position', [1000 500 410 860], 'PaperPositionMode', 'auto');
                    suptitle(['session ' num2str(n) ', cell ' num2str(c)]);
                    hdgs_fine = hdgs(1):0.1:hdgs(end);
                    subplot(3,1,m);
                    errorbar(hdgs,FRmean,FRse,[clr{m,o} 'o-']); hold on;
                    plot(hdgs_fine,tuning_vonMises(beta,hdgs_fine),[clr{m,o} '--']);
                    if m==2; ylabel('firing rate (sp/s)'); end
                    if m==3; xlabel('heading (deg)'); end
                    title(mlabel{m});
                    changeAxesFontSize(gca, 16, 16);
                    yl = ylim; xl = xlim;
                    text(xl(1)+0.02*(xl(2)-xl(1)),yl(1)+0.08*(yl(2)-yl(1)),sprintf('p = %1.2f',pvalTuning(n,c,m,o)));
                    keyboard
                end
                
                % fit a line to middle few headings, to establish a slope-around-zero criterion
                J = I & abs(hdg)<20;
                [P,S] = polyfit(hdg(J),spRate_tun(J),1);
                se = sqrt(diag(inv(S.R)*inv(S.R')).*S.normr.^2./S.df);
                t = abs(P(1)/se(1));
                slope(C,m) = P(1);
                pval_ttest_slope(C,m) = 2*(1-tcdf(t,S.df)); % two-tailed
                
                % identify pref direction across all headings [for referencing CPs below]
                % here we want to avoid using task data, because Zaidel.
                % so just compare spike rates on R vs L trials during tuning block
% %                 pref(m) = double(nanmean(spRate_tun(hdg>0)) > nanmean(spRate_tun(hdg<0)));
                % OR use slope of line fit to the middle ~3 hdgs
                pref(m) = (sign(P(1))+1)/2;
                
                
                if plotTuningPSTH
                    figure(n*100+c); set(gcf, 'Color', 'w', 'Position', [400 500 800 800], 'PaperPositionMode', 'auto');                                
                    suptitle(['session ' num2str(n) ', cell ' num2str(c)]);
                    subplot(3,1,m);
                    if sum(FRmean)>5
                        prefHdg = hdgs(FRmean==max(FRmean)); prefHdg = prefHdg(1);
                        nullHdg = hdgs(FRmean==min(FRmean)); nullHdg = nullHdg(1);
                        psth = calc_mean(squeeze(raster_tun(C,:,:)), hdg, [prefHdg;nullHdg])*1e3;
                        psth = smoothRaster(psth, convKernel);
                        tAxis = -extRaster(1) : size(squeeze(raster_tun(C,:,:)),2) - extRaster(1) - 1;
                        tAxis = tAxis + length(convKernel)/2+1; % offset so that filter is 'causal' (spikes do not affect spike density backward in time)
                        YLim = [min(min(psth))*0.9 max(max(psth))*1.1];
                        XLim = [tAxis(1) 2000];
                        hold on;
                        for p = 1:size(psth,1)
                            plot(tAxis, psth(p,1:length(tAxis)),'Color', clr2{p},'LineWidth',2);
                        end
                        set(gca, 'XLim', XLim, 'YLim', YLim, 'TickDir','out'); box off; % 'xtick', 0:100:400, 'xticklabel', {'0','','200','','400'},
                        if m==1; ylabel('ves'); end
                        if m==2; ylabel(sprintf('Firing rate (spikes/s)\nvis')); end
                        if m==3; ylabel('comb'); xlabel('Time from stim on (ms)'); end
                        changeAxesFontSize(gca, 16, 16);
                    end
                end               
%             end
        end
        


        
    %*****************************
    % TASK (exp)
    %*****************************
    
        ev = newDataStruct(n).data.dots3DMP.events;
        good = logical(ev.goodtrial');
        un = newDataStruct(n).data.dots3DMP.units;
        
%         mStart = ev.stimOn(good)';
%         mEnd = ev.stimOff(good)';
        % OR
        mStart = ev.motionOn(good)';
        mEnd = ev.postTargHold(good)';
        
        dur = round((mEnd-mStart)*1000);
        mod = ev.modality(good)';
        mods = unique(mod);
        coh = ev.coherence(good)';
        cohs = unique(coh);
            if length(cohs)>3; keyboard; end
            if length(cohs)==3
                warning('three cohs?');
                coh(coh==cohs(2)) = cohs(3); % pool largest two cohs, at least for tuning; this shouldn't happen often
            end
        cohs = unique(coh);                
        hdg = ev.heading(good)';
        hdgs = unique(hdg);
        
        spikeTimes = un.spiketimes{c};
%         invalidTr = dur > 2000; % exclude super long RTs which must have been invalid 
        startTrial = mStart-extRaster(1)/1000;
        endTrial = mEnd+extRaster(2)/1000;
        trialLen = endTrial - startTrial; % in seconds
        
        nTr = length(mStart);
        spRate_exp = nan(nTr,1);
        spCountFixedWin = nan(nTr,1);
        
        for t = 1:nTr
%             if dur>1400; continue; end % exclude super long RTs which must have been invalid 
            % raster, align stimOn (& spike counts/rates)
            if round(trialLen(t)*1000)<size(raster_stimOn,3)
                raster_stimOn(C,t,1:round(trialLen(t)*1000)) = 0; % init bins to zero but keep nans after RT+extRaster(2)
                spkInd = spikeTimes >= startTrial(t) & spikeTimes <= endTrial(t);
                spkTimes_thisTr = spikeTimes(spkInd) - startTrial(t);                
                if isempty(spkTimes_thisTr) % is there a possibility that no spikes means invalid trial that should be NaN? need to consider this
                    spRate_exp(t) = 0;
                    spCount_exp(t) = 0;
                else
                    x = 0:0.001:ceil(spkTimes_thisTr(end)*1000)/1000;
                    raster_stimOn(C,t,1:length(x)-1) = histcounts(spkTimes_thisTr,x);
                    a = latency(1)+extRaster(1); % this will be just 'latency' because raster has been offset by extRaster
                    b = a + dur(t) - latency(1) + latency(2);
                    spRate_exp(t) = sum(raster_stimOn(C,t,a:b)) / ((b-a+1)/1000);
%                     spCountFixedWin(t) = sum(raster_stimOn(C,t,a:a+minRT-latency(1)-20)); % try to avoid any interference with end of trial
                end
            end % else, trial exceeds the max length specified in prealloc (will stay nans)
                        
            % raster, align RT (stimOff)
            if round(trialLen(t)*1000)<size(raster_RT,3)
                raster_RT(C,t,end-round(trialLen(t)*1000):end) = 0; % init bins to zero but keep nans before stimOn-extRaster(1)
                spkTimes_thisTr = flipud(endTrial(t)-spikeTimes(spkInd)); % flipud instead, makes hist easier
                if ~isempty(spkTimes_thisTr) && spkTimes_thisTr(end)>0.001
                    x = 0:0.001:ceil(spkTimes_thisTr(end)*1000)/1000;
                    raster_RT(C,t,end-length(x)+2:end) = fliplr(histcounts(spkTimes_thisTr,x));
                end            
            end            
        end
        
        
% %     % identify pref choice targ [for referencing CPs below]
% %         % two options here: fit of mapping data, or just compare spike
% %         % rates on R and L trials during task (try both, but I think
% %         % mapping data makes the most sense in terms of what the readout
% %         % 'knows'; task data can be biased by, well, biases in RT or conf)        
% % %         pref = double(prefDir_fit(C)<90 || prefDir_fit(C)>270);
% %         % OR
% % %         I = dataCell{n}.Exp.openEvents.direction==0 & dataCell{n}.Exp.openEvents.coherence > 0.25;   % high coh but some sess might exclude 0.512, 
% % %         J = dataCell{n}.Exp.openEvents.direction==180 & dataCell{n}.Exp.openEvents.coherence > 0.25; % plus we get more s/n this way
% % %         pref = double(nanmean(spRate_exp(I)) > nanmean(spRate_exp(J)));

    % see above; now defined using tuning data, not task data!


%************************************************************************
% CP, confP, and choice/conf-conditioned PSTHs
%************************************************************************
        
        % grab behavioral data
        choice = ev.choice(good)';
        pdw = ev.PDW(good)';
        RT = round(ev.RT(good)*1000)';
        if max(choice)==2 % convert choices to 0:1
            choice = choice-1;
        end
        if any(isnan(choice))
            keyboard
        end
        if any(isnan(pdw))
            keyboard
        end
        
        % PSTHs now go first, so we can adaptively select the window for
        % spike rate for CP based on the peak time (NOT DONE)
% % %         % *****ON HOLD*****
% % %         % start with all headings, to find peak time (for spRate)
% % %         % but still enforce a min RT
% % %         J = RT >= minRT;
% % %         % aligned stim-on
% % %         clear psth
% % %         % first, irrespective of choice/PDW, to subtract for residuals
% % %         wideKernel = fspecial('average', [1 100]);
% % %         psthAll = smoothRaster(nanmean(squeeze(raster_stimOn(C,J,:)))*1e3, wideKernel);
% % %         % *****ON HOLD*****

        for m = 1:length(mods)
        
        % for conditioned PSTH, use only small/zero hdgs
        % or else there's a confound (later we'll pool across hdg)
        J = mod==mods(m) & abs(hdg)<=maxHdg & RT >= minRT;

        % aligned stim-on
        clear psth
        % first, irrespective of choice/PDW, to subtract for residuals
        psthAll = smoothRaster(nanmean(squeeze(raster_stimOn(C,J,:)))*1e3, convKernel);
        % now condition on each outcome
        I = J & choice==pref(m) & pdw==1; % pref-high
        psth(1,:) = smoothRaster(nanmean(squeeze(raster_stimOn(C,I,:)))*1e3, convKernel);        
        I = J & choice==pref(m) & pdw==0; % pref-low
        psth(2,:) = smoothRaster(nanmean(squeeze(raster_stimOn(C,I,:)))*1e3, convKernel);
        I = J & choice~=pref(m) & pdw==1; % null-high
        psth(3,:) = smoothRaster(nanmean(squeeze(raster_stimOn(C,I,:)))*1e3, convKernel);        
        I = J & choice~=pref(m) & pdw==0; % null-low
        psth(4,:) = smoothRaster(nanmean(squeeze(raster_stimOn(C,I,:)))*1e3, convKernel);
        
        rMax = max(max(psth(:,1:psth_xlim1)));
        psthNorm_stimOn_one{m}(C,:) = psth(1,:)/rMax;
        psthNorm_stimOn_two{m}(C,:) = psth(2,:)/rMax;
        psthNorm_stimOn_three{m}(C,:) = psth(3,:)/rMax;
        psthNorm_stimOn_four{m}(C,:) = psth(4,:)/rMax;
        % and the residuals
        psthResid_stimOn_one{m}(C,:) = (psth(1,:)-psthAll)/rMax;
        psthResid_stimOn_two{m}(C,:) = (psth(2,:)-psthAll)/rMax;
        psthResid_stimOn_three{m}(C,:) = (psth(3,:)-psthAll)/rMax;
        psthResid_stimOn_four{m}(C,:) = (psth(4,:)-psthAll)/rMax;
        
        if plotTaskPSTH && c<5
            figure(n*1000); set(gcf, 'Color', 'w', 'Position', [900 425 600 1000], 'PaperPositionMode', 'auto');
            subplot(4,1,c);
            tAxis = -extRaster(1) : psth_xlim1;
            tAxis = tAxis + length(convKernel)/2+1; % offset for causal
            clr = {'g-', 'g--', 'r-', 'r--'};            
            YLim = [min(min(psth(:,1:length(tAxis))))*0.9 max(max(psth(:,1:length(tAxis))))*1.1];
            XLim = [tAxis(1) psth_xlim1];
            hold on;
            for p = 1:size(psth,1)
                plot(tAxis, psth(p,1:length(tAxis)), clr{p},'LineWidth',2);
            end
            set(gca, 'XLim', XLim, 'YLim', YLim, 'TickDir','out'); box off;
            xlabel('Time from motion on (ms)');
            ylabel('Firing rate (spikes/s)');
            changeAxesFontSize(gca, 20, 20);
            if c==nU(n)
                legend('pref high','pref low','null high','null low','Location','Northeast'); legend('boxoff');
            end
        end

        
        % aligned RT
        clear psth
        % first, irrespective of choice/PDW, to subtract for residuals
        psthAll = smoothRaster(nanmean(squeeze(raster_RT(C,J,:)))*1e3, convKernel);
        % now condition on each outcome
        I = J & choice==pref(m) & pdw==1; % pref-high
        psth(1,:) = smoothRaster(nanmean(squeeze(raster_RT(C,I,:)))*1e3, convKernel);
        I = J & choice==pref(m) & pdw==0; % pref-low
        psth(2,:) = smoothRaster(nanmean(squeeze(raster_RT(C,I,:)))*1e3, convKernel);
        I = J & choice~=pref(m) & pdw==1; % null-high
        psth(3,:) = smoothRaster(nanmean(squeeze(raster_RT(C,I,:)))*1e3, convKernel);
        I = J & choice~=pref(m) & pdw==0; % null-low
        psth(4,:) = smoothRaster(nanmean(squeeze(raster_RT(C,I,:)))*1e3, convKernel);
        psthNorm_RT_one{m}(C,:) = psth(1,:)/rMax;
        psthNorm_RT_two{m}(C,:) = psth(2,:)/rMax;
        psthNorm_RT_three{m}(C,:) = psth(3,:)/rMax;
        psthNorm_RT_four{m}(C,:) = psth(4,:)/rMax;
        % and the residuals
        psthResid_RT_one{m}(C,:) = (psth(1,:)-psthAll)/rMax;
        psthResid_RT_two{m}(C,:) = (psth(2,:)-psthAll)/rMax;
        psthResid_RT_three{m}(C,:) = (psth(3,:)-psthAll)/rMax;
        psthResid_RT_four{m}(C,:) = (psth(4,:)-psthAll)/rMax;

        if plotTaskPSTH && c<5
            figure(n*10000); set(gcf, 'Color', 'w', 'Position', [1000 125 600 1000], 'PaperPositionMode', 'auto');
            subplot(4,1,c);
            tAxis = -(size(raster_RT,3)-extRaster(2)-1) : extRaster(2);
            tAxis = tAxis + length(convKernel)/2+1; % offset for causal
            XLim = psth_xlim2;
            YLim = [min(min(psth(:,end-(psth_xlim2(2)-psth_xlim2(1)):end)))*0.9 max(max(psth(:,end-(psth_xlim2(2)-psth_xlim2(1)):end)))*1.1];
            hold on;
            for p = 1:size(psth,1)
                plot(tAxis, psth(p,:), clr{p},'LineWidth',2);
            end
            set(gca, 'XLim', XLim, 'YLim', YLim, 'TickDir','out', 'YAxisLocation', 'Right'); box off;
            xlabel('Time from saccade (ms)');
            changeAxesFontSize(gca, 20, 20);
        end
        
        
        % Next, conventional CP from full-trial spike rates, referenced to
        % preferred direction (as defined by the task trials themselves),
        % rather than Haefner-style which is to reference to a fixed choice.
        % We use rates instead of counts because in RT task, using counts
        % would overweight long-RT trials. Also, the choice-conditioned
        % count distributions may not be equal under the null hypothesis,
        % if RT isn't perfectly symmetric for the two choices.
        
        % which data?
%         respForCP = spCountFixedWin;
        respForCP = spRate_exp;
%         respForCP = spCountPeak; % NOT READY YET
        
        % there are 4 choices available, so CP and ConfP each come in 2 flavors:
        I = J & choice==pref(m) & pdw==1;
            X = respForCP(I); % pref-high        
        I = J & choice~=pref(m) & pdw==1;
            Y = respForCP(I); % null-high
        try
            CP_high(C,m) = rocN(X(~isnan(X)),Y(~isnan(Y)),100);
        catch
            CP_high(C,m) = NaN;
        end
                
        I = J & choice==pref(m) & pdw==0;
            X = respForCP(I); % pref-low        
        I = J & choice~=pref(m) & pdw==0;
            Y = respForCP(I); % null-low
        try
            CP_low(C,m) = rocN(X(~isnan(X)),Y(~isnan(Y)),100);
        catch
            CP_low(C,m) = NaN;
        end
        
        I = J & choice==pref(m);
            X = respForCP(I); % pref-all        
        I = J & choice~=pref(m);
            Y = respForCP(I); % null-all
        try
            CP_all(C,m) = rocN(X(~isnan(X)),Y(~isnan(Y)),100);
        catch
            CP_all(C,m) = NaN;
        end
        
        I = J & choice==pref(m) & pdw==1;
            X = respForCP(I); % pref-high        
        I = J & choice==pref(m) & pdw==0;
            Y = respForCP(I); % pref-low
        try
            ConfP_pref(C,m) = rocN(X(~isnan(X)),Y(~isnan(Y)),100);
        catch
            ConfP_pref(C,m) = NaN;
        end
        
        I = J & choice~=pref(m) & pdw==1;
            X = respForCP(I); % null-high        
        I = J & choice~=pref(m) & pdw==0;
            Y = respForCP(I); % null-low
        try
            ConfP_null(C,m) = rocN(X(~isnan(X)),Y(~isnan(Y)),100);
        catch
            ConfP_null(C,m) = NaN;
        end
        
        I = J & pdw==1;
            X = respForCP(I); % all-high        
        I = J & pdw==0;
            Y = respForCP(I); % all-low
        try
            ConfP_all(C,m) = rocN(X(~isnan(X)),Y(~isnan(Y)),100);
        catch
            ConfP_all(C,m) = NaN;
        end        
        
% % grab an example pair of histograms, eg for talk figure
% if n==17 && c==3
%     keyboard
%     figure; set(gcf, 'Color', 'w', 'Position', [800 100 450 300], 'PaperPositionMode', 'auto');
%     [Nx,~] = histcounts(X,0:5:70);
%     [Ny,E] = histcounts(Y,0:5:70);
%     bar(E(2:end)-2.5,[Nx;Ny]',.9); box('off');
%     set(gca,'YLim',[0 12],'ytick',0:2:12,'XLim',[0 70],'xtick',0:10:70,'TickDir','out');
%     legend('high bet','low bet'); legend('boxoff')
%     xlabel('Firing rate (sp/s)'); ylabel('Count');
%     changeAxesFontSize(gca, 16, 16);
% end
        


%         % Next, repeat the above but separately for each coh (within each
%         % cell), then take a weighted average across most cohs. This was an
%         % original strategy of Britten et al: even with a strong sensory
%         % signal, the choice-conditioned distributions (within a coh) are
%         % overlapping under the null hypothesis. Of course, usually cannot
%         % use the 1-2 highest cohs because not enough error trials.
%         
%         keyboard
%         
%         for n = 1
%         
%             
%         end

        end % mods

        
        sessID(C) = n;
        C = C+1;
        
    end
    
%     pause
end



%% save, to save time later

clear newDataStruct raster_stimOn raster_RT raster_tun
save('/Users/chris/Documents/MATLAB/temp.mat');



%% reload

% clear
% load(['/Users/chris/Documents/MATLAB/temp.mat']);



%% plot average choice/conf-conditioned PSTH


for m = 1:length(mods)

    % aligned stimOn
clear psth

I = nansum(psthNorm_stimOn_one{m},2)>10; % pick some selection criteria; for now it's all cells (w nonzero spike count) 

psth(1,:) = nanmean(psthNorm_stimOn_one{m}(I,:));
psth(2,:) = nanmean(psthNorm_stimOn_two{m}(I,:));
psth(3,:) = nanmean(psthNorm_stimOn_three{m}(I,:));
psth(4,:) = nanmean(psthNorm_stimOn_four{m}(I,:));

figure(5000); set(gcf, 'Color', 'w', 'Position', [1 220 700 450], 'PaperPositionMode', 'auto');
tAxis = -extRaster(1) : psth_xlim1;
tAxis = tAxis + length(convKernel)/2+1; % offset for causal
clr = {'b-', 'b--', 'r-', 'r--'};            
YLim = [min(min(psth(:,1:length(tAxis))))*0.9 max(max(psth(:,1:length(tAxis))))*1.1];
XLim = [tAxis(1) psth_xlim1];
hold on;
for p = 1:size(psth,1)
    plot(tAxis, psth(p,1:length(tAxis)), clr{p},'LineWidth',2);
end
set(gca, 'XLim', XLim, 'YLim', YLim, 'TickDir','out'); box off;
xlabel('Time from motion on (ms)');
ylabel('Firing rate (normalized)');
changeAxesFontSize(gca, 22, 22);
h = legend('pref choice, high bet','pref choice, low bet','null choice, high bet','null choice, low bet','Location','Northwest'); legend('boxoff');
saveas(gcf, '3a', 'pdf')

    % aligned RT
clear psth
psth(1,:) = nanmean(psthNorm_RT_one(I,:));
psth(2,:) = nanmean(psthNorm_RT_two(I,:));
psth(3,:) = nanmean(psthNorm_RT_three(I,:));
psth(4,:) = nanmean(psthNorm_RT_four(I,:));

figure(5001); set(gcf, 'Color', 'w', 'Position', [350 220 700 450], 'PaperPositionMode', 'auto');
% subplot(1,2,2);
tAxis = -(rasterLen_exp+sum(extRaster)-extRaster(2)-1) : extRaster(2);
tAxis = tAxis + length(convKernel)/2+1; % offset for causal
XLim = [-400 20];
YLim = [min(min(psth(:,end-(psth_xlim2(2)-psth_xlim2(1)):end)))*0.9 max(max(psth(:,end-(psth_xlim2(2)-psth_xlim2(1)):end)))*1.1];
hold on;
for p = 1:size(psth,1)
    plot(tAxis, psth(p,:), clr{p},'LineWidth',2);
end
set(gca, 'XLim', XLim, 'YLim', YLim, 'TickDir','out', 'YAxisLocation', 'Right'); box off;
xlabel('Time from saccade (ms)');
changeAxesFontSize(gca, 22, 22);
% h = legend('pref choice, high bet','pref choice, low bet','null choice, high bet','null choice, low bet','Location','Northwest'); legend('boxoff');
saveas(gcf, '3b', 'pdf')



%% plot average choice/conf-conditioned PSTH, RESIDUALS

    % aligned dotsOn
clear psth

I = nansum(psthNorm_stimOn_one,2)>10; % pick some selection criteria; for now it's all cells (w nonzero spike count) 

psth(1,:) = nanmean(psthResid_stimOn_one{m}(I,:));
psth(2,:) = nanmean(psthResid_stimOn_two{m}(I,:));
psth(3,:) = nanmean(psthResid_stimOn_three{m}(I,:));
psth(4,:) = nanmean(psthResid_stimOn_four{m}(I,:));

figure(5100); set(gcf, 'Color', 'w', 'Position', [1 220 700 450], 'PaperPositionMode', 'auto');
tAxis = -extRaster(1) : psth_xlim1;
tAxis = tAxis + length(convKernel)/2+1; % offset for causal
clr = {'b-', 'b--', 'r-', 'r--'};            
XLim = [tAxis(1) psth_xlim1];
YLim = [-0.15 0.15];
hold on;
for p = 1:size(psth,1)
    plot(tAxis, psth(p,1:length(tAxis)), clr{p},'LineWidth',2);
end
set(gca, 'XLim', XLim, 'YLim', YLim, 'TickDir','out'); box off;
xlabel('Time from dots on (ms)');
ylabel('Residual firing rate (normalized)');
changeAxesFontSize(gca, 22, 22);
% h = legend('pref choice, high bet','pref choice, low bet','null choice, high bet','null choice, low bet','Location','Northeast'); legend('boxoff');
saveas(gcf, '4a', 'pdf')

    % aligned RT
clear psth
psth(1,:) = nanmean(psthResid_RT_one(I,:));
psth(2,:) = nanmean(psthResid_RT_two(I,:));
psth(3,:) = nanmean(psthResid_RT_three(I,:));
psth(4,:) = nanmean(psthResid_RT_four(I,:));

figure(5101); set(gcf, 'Color', 'w', 'Position', [350 220 700 450], 'PaperPositionMode', 'auto');
% subplot(1,2,2);
tAxis = -(rasterLen_exp+sum(extRaster)-extRaster(2)-1) : extRaster(2);
tAxis = tAxis + length(convKernel)/2+1; % offset for causal
XLim = [-400 20];
YLim = [-0.15 0.15];
hold on;
for p = 1:size(psth,1)
    plot(tAxis, psth(p,:), clr{p},'LineWidth',2);
end
set(gca, 'XLim', XLim, 'YLim', YLim, 'TickDir','out', 'YAxisLocation', 'Right'); box off;
xlabel('Time from saccade (ms)');
changeAxesFontSize(gca, 22, 22);
saveas(gcf, '4b', 'pdf')



%% CP and ConfP, vs each other and as a func of prefDir

nanmean(CP_high(:,m))
nanmean(CP_low(:,m))
nanmean(CP_all(:,m))
nanmean(ConfP_pref(:,m))
nanmean(ConfP_null(:,m))
nanmean(ConfP_all(:,m))


%% choose which to plot

CP = CP_high;
% or
% CP = CP_low;
% or
% CP = CP_all;

ConfP = ConfP_pref;
% or
% ConfP = ConfP_null;
% or
% ConfP = ConfP_all;

nans = isnan(CP) | isnan(ConfP);

CP(nans) = [];
ConfP(nans) = [];

CPse = std(CP)/sqrt(length(CP));
[~,P] = ttest(CP,0.5)

ConfPse = std(ConfP)/sqrt(length(ConfP));
[~,P] = ttest(ConfP,0.5)


% scatter whole-trial CP vs ConfP

figure; set(gcf, 'Color', 'w', 'Position', [400 200 450 450], 'PaperPositionMode', 'auto');
plot([0 1],[0 1],'k--',[0 1],[0.5 0.5],'k:',[0.5 0.5],[0 1],'k:');
hold on; plot(CP,ConfP,'ko','MarkerSize',11,'MarkerFaceColor',[1 1 1]); axis square;
% set(gca,'XLim',[0.17 0.83],'YLim',[0.17 0.83],'xtick',0.2:0.1:0.8,'ytick',0.2:0.1:0.8,'TickDir','out')
set(gca,'XLim',[0.1 0.9],'YLim',[0.1 0.9],'xtick',0.1:0.2:0.9,'ytick',0.1:0.2:0.9,'TickDir','out')
xlabel('Choice probability'); ylabel('Confidence probability');
changeAxesFontSize(gca, 22, 22);
saveas(gcf, '5a', 'pdf')

figure; set(gcf, 'Color', 'w', 'Position', [800 200 450 80], 'PaperPositionMode', 'auto');
[histCP,E] = histcounts(CP,0:0.05:1);
bar(E(2:end)-0.025,histCP); box('off');
set(gca,'XLim',[0.1 0.9],'xtick',0.1:0.2:0.9,'xticklabel',[],'TickDir','out');
changeAxesFontSize(gca, 14, 14);
saveas(gcf, '5b', 'pdf')

figure; set(gcf, 'Color', 'w', 'Position', [800 100 80 450], 'PaperPositionMode', 'auto');
[histCP,E] = histcounts(ConfP,0:0.05:1);
barh(E(2:end)-0.025,histCP); box('off');
set(gca,'YLim',[0.1 0.9],'ytick',0.1:0.2:0.9,'yticklabel',[],'TickDir','out');
changeAxesFontSize(gca, 14, 14);
saveas(gcf, '5c', 'pdf')

[r,p] = corrcoef(CP,ConfP)





end % mods






