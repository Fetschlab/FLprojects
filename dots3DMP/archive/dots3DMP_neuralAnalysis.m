% NEW: load/plot data from Miguel struct

clear
close all
load('/Users/chris/Documents/MATLAB/lucio_neuraldatacell_20220606.mat');
dbstop if error

%%
% dots3DMP_plotBehavFromNeuralStruct % need to fix, later


%%

extRaster = [50 150]; % extend raster this much before motion onset (1) and after motion offset (2)
latency = [300 -200]; % cutoffs for spike rate [post stimulus onset (1) and post stimulus offset (2)];
                      % latency(2) must be < extRaster(2)!
convKernel = fspecial('average', [1 40]); % N ms wide boxcar (acausal, centered)

% some manual plotting lims
psth_xlim1 = 400;
psth_xlim2 = [-300 100];

minRT = 0; 

excludes = [];

dots3DMP_neuralAnalysis_prelims

%%
%%%%%%%%%%%%%%%%
% BEGIN MAIN ANALYSIS LOOP: calc rasters, rates/counts, and psths for each unit
%%%%%%%%%%%%%%%%
C=1; % counter for valid units, across all sessions

for n = 1:length(dataCell)

%     if plotEverythingElse        
%         figure(n); set(gcf, 'Color', 'w', 'Position', [800 725 370*min([nU(n) 8]) 430], 'PaperPositionMode', 'auto');
%         figure(n*10); set(gcf, 'Color', 'w', 'Position', [900 425 370*min([nU(n) 8]) 430], 'PaperPositionMode', 'auto');
%         figure(n*100); set(gcf, 'Color', 'w', 'Position', [1000 125 370*min([nU(n) 8]) 430], 'PaperPositionMode', 'auto');
%     end
    
    for c = 1:nU(n)
        
        if ~isempty(excludes)
            if ismember([n c],excludes,'rows')
                continue
            end
        end
        
        disp(num2str([n c]));
        
    %*****************************
    % TUNING
    %*****************************
    
        mStart = dataCell{n}.Mapping.openEvents.motionStart';
        mEnd = dataCell{n}.Mapping.openEvents.motionEnd';
        dur = round((mEnd-mStart)*1000); % some variability (due to dropped frames???), but okay for now
        dir = dataCell{n}.Mapping.openEvents.direction;
        dirs = unique(dir);
        
        spikeTimes = dataCell{n}.Mapping.spikeTimes{c};
        startTrial = mStart-extRaster(1)/1000;
        endTrial = mEnd+extRaster(2)/1000;
        trialLen = endTrial - startTrial; % in seconds
        duration = dur;
        nTr = nTrials_map(n);
        spRate_map = nan(nTr,1);
        
        for t = 1:nTr
            if round(trialLen(t)*1000)<size(raster_map,3) 
                raster_map(C,t,1:round(trialLen(t)*1000)) = 0; % init bins to zero but keep nans after extRaster(2) cutoff
                spkInd = spikeTimes >= startTrial(t) & spikeTimes <= endTrial(t);
                spkTimes = spikeTimes(spkInd) - startTrial(t);                
                if isempty(spkTimes)
                    spRate_map(t) = 0;
                else
                    x = 0:0.001:ceil(spkTimes(end)*1000)/1000; % bin spikes only up to the last one; after that will be zeros, followed by NaNs, with a cutoff determined above (see comment)
                    raster_map(C,t,1:length(x)-1) = histcounts(spkTimes,x);
                    a = latency(1)+extRaster(1);
                    b = a + duration(t) - latency(1) + latency(2);
                    spRate_map(t) = sum(raster_map(C,t,a:b)) / ((b-a+1)/1000);
                end
            end % else, trial exceeds the max length specified in prealloc, likely invalid, loss of eye signal etc (will stay nans)
            
            a = latency(1)+extRaster(1);
            b = a + duration(t) - latency(1) + latency(2);
        end

        if sum(isnan(spRate_map(length(dir)-5:length(dir)))) > 3; keyboard; end % possible mismatch in ntrials given by dir, vs. spRate

        % fit the tuning curves to get pref dir
        FRmean = calc_mean(spRate_map(1:length(dir)),dir,dirs);
        FRmean(end+1) = FRmean(1);
        dirs(end+1) = 360;
        peak = dirs(FRmean==max(FRmean)); peak = peak(1);
            % vonMises params: b = [ampl, kappa, theta, offset]
        guess = [max(FRmean)-min(FRmean), 1.5, peak, min(FRmean)];
        [beta, fval] = fminsearch(@(x) tuning_vonMises_err(x,dir,spRate_map(1:length(dir))), guess);
        
        if plotTuningFits
               % temp, to check fits by eye: (they're all good!)
            gdirs = dirs(1):dirs(end);
            figure(1200); set(gcf, 'Color', 'w', 'Position', [1181 555 500 400], 'PaperPositionMode', 'auto'); clf; 
            plot(dirs,FRmean,'bo-',gdirs,tuning_vonMises(beta,gdirs),'r-');
            set(gca,'xtick',0:45:360,'xlim',[0 360]);
            title(num2str([n c]));
            pause
        end
        
  %**************
  %repeat for Remap
  %**************
  
        if ismember([n c],useRemap,'rows')                        

            if size(dataCell{n}.Remap.spikeTimes,2) ~= size(dataCell{n}.Mapping.spikeTimes,2)
                keyboard % channel mismatch?
            end
            
            mStart = dataCell{n}.Remap.openEvents.motionStart';
            mEnd = dataCell{n}.Remap.openEvents.motionEnd';
            dur = round((mEnd-mStart)*1000); % some variability (due to dropped frames???), but okay for now
            dir = dataCell{n}.Remap.openEvents.direction;
            dirs = unique(dir);

            spikeTimes = dataCell{n}.Remap.spikeTimes{c};
            startTrial = mStart-extRaster(1)/1000;
            endTrial = mEnd+extRaster(2)/1000;
            trialLen = endTrial - startTrial; % in seconds
            duration = dur;
            nTr = length(mStart);
            spRate_map = nan(nTr,1);

            for t = 1:nTr
                if round(trialLen(t)*1000)<size(raster_map,3) 
                    raster_map(C,t,1:round(trialLen(t)*1000)) = 0; % init bins to zero but keep nans after extRaster(2) cutoff
                    spkInd = spikeTimes >= startTrial(t) & spikeTimes <= endTrial(t);
                    spkTimes = spikeTimes(spkInd) - startTrial(t);                
                    if isempty(spkTimes)
                        spRate_map(t) = 0;
                    else
                        x = 0:0.001:ceil(spkTimes(end)*1000)/1000; % bin spikes only up to the last one; after that will be zeros, followed by NaNs, with a cutoff determined above (see comment)
                        raster_map(C,t,1:length(x)-1) = histcounts(spkTimes,x);
                        a = latency(1)+extRaster(1);
                        b = a + duration(t) - latency(1) + latency(2);
                        spRate_map(t) = sum(raster_map(C,t,a:b)) / ((b-a+1)/1000);
                    end
                end % else, trial exceeds the max length specified in prealloc, likely invalid, loss of eye signal etc (will stay nans)

                a = latency(1)+extRaster(1);
                b = a + duration(t) - latency(1) + latency(2);
            end

            if sum(isnan(spRate_map(length(dir)-5:length(dir)))) > 3; keyboard; end % possible mismatch in ntrials given by dir, vs. spRate

            % fit the tuning curves to get pref dir
            FRmean = calc_mean(spRate_map(1:length(dir)),dir,dirs);
            FRmean(end+1) = FRmean(1); dirs(end+1) = 360;
            peak = dirs(FRmean==max(FRmean)); peak = peak(1);
                % vonMises params: b = [ampl, kappa, theta, offset]
            guess = [max(FRmean)-min(FRmean), 1.5, peak, min(FRmean)];
            [beta, fval] = fminsearch(@(x) tuning_vonMises_err(x,dir,spRate_map(1:length(dir))), guess);

            if plotTuningFits
                   % temp, to check fits by eye: (they're all good!)
                gdirs = dirs(1):dirs(end);
                figure(1201); set(gcf, 'Color', 'w', 'Position', [1181 81 500 400], 'PaperPositionMode', 'auto'); clf; 
                plot(dirs,FRmean,'bo-',gdirs,tuning_vonMises(beta,gdirs),'r-');
                set(gca,'xtick',0:45:360,'xlim',[0 360]);
                title(num2str([n c]));
                pause
            end
                
        end

  %**************
  %end repeat for remap
  %**************
        
        % so beta is either from mapping or was overwritten by remap, and:
        prefDir_fit(C) = mod(beta(3),360); % sometimes it's negative, who knows why, just wrap it back around 360
        
        if plotEverythingElse && c<9
            figure(n);
            % simple tuning curve
            subplot(2,nU(n),c);
            dirPlot = dir; dirsPlot = dirs;
            [FRmean, FRse, ~] = calc_mean(spRate_map(1:length(dirPlot)),dirPlot,dirsPlot);
            errorbar(dirsPlot,FRmean,FRse,'ko-','MarkerFaceColor','k');
            set(gca,'xtick',dirsPlot,'xticklabel',{'0','','90','','180','','270',''},'tickdir','out')
            xlabel('direction (deg)'); ylabel('firing rate (sp/s)');
        %         title(sprintf('sess %d, unit %d',n,dataCell{n}.Mapping.Unit(c)));
            title(sprintf('sess %d, unit %d',n,c));
            changeAxesFontSize(gca, 20, 20);

            % PSTH for pref/null, aligned to dots onset
            subplot(2,nU(n),c+nU(n));
            pref = dirsPlot(FRmean==max(FRmean)); null = dirsPlot(FRmean==min(FRmean));
            psth = calc_mean(squeeze(raster_map(C,:,:)), dirPlot, [pref;null])*1e3;
            psth = smoothRaster(psth, convKernel);
            tAxis = -extRaster(1) : size(squeeze(raster_map(C,:,:)),2) - extRaster(1) - 1;
            tAxis = tAxis + length(convKernel)/2+1; % offset for causal
            clr = {'g','r'};
            YLim = [min(min(psth))*0.9 max(max(psth))*1.1];
            XLim = [tAxis(1) 350];
            hold on;
            for p = 1:size(psth,1)
                plot(tAxis, psth(p,1:length(tAxis)),'Color', clr{p},'LineWidth',2);
            end
            set(gca, 'XLim', XLim, 'YLim', YLim, 'xtick', 0:100:400, 'xticklabel', {'0','','200','','400'}, 'TickDir','out'); box off;
            xlabel('Time from dots on (ms)');
            ylabel('Firing rate (spikes/s)');
            changeAxesFontSize(gca, 22, 22);
%             pause;
        end


    %*****************************
    % TASK (exp)
    %*****************************

        mStart = dataCell{n}.Exp.openEvents.motionStart';
        nTrials = length(mStart);
        mEnd = dataCell{n}.Exp.openEvents.motionEnd';
        RT = round((mEnd-mStart)*1000);
        dur=RT;
% % %             invalidTr = dur > 2400; % exclude super long RTs which must have been invalid 
        dir = dataCell{n}.Exp.openEvents.direction;
        dirs = unique(dir);
        
        startTrial = mStart-extRaster(1)/1000;
        endTrial = mEnd+extRaster(2)/1000;
        trialLen = endTrial - startTrial; % in seconds
        nTr = nTrials_exp(n);
        spRate_exp = nan(nTr,1);
        spCount_exp = nan(nTr,1);
        
        for t = 1:nTr
            % raster, align dotsOn (& spike counts/rates)
            if round(trialLen(t)*1000)<size(raster_dotsOn,3)
                raster_dotsOn(C,t,1:round(trialLen(t)*1000)) = 0; % init bins to zero but keep nans after RT+extRaster(2)
                spkInd = dataCell{n}.Exp.spikeTimes{c} >= startTrial(t) & dataCell{n}.Exp.spikeTimes{c} <= endTrial(t);
                spkTimes = dataCell{n}.Exp.spikeTimes{c}(spkInd) - startTrial(t);                
                if isempty(spkTimes) % is there a possibility that no spikes means invalid trial that should be NaN? need to consider this w Miguel ********
                    spRate_exp(t) = 0;
                    spCount_exp(t) = 0;
                else
                    x = 0:0.001:ceil(spkTimes(end)*1000)/1000;
                    raster_dotsOn(C,t,1:length(x)-1) = histcounts(spkTimes,x);
                    a = latency(1)+extRaster(1);
                    b = a + dur(t) - latency(1) + latency(2);
                    spRate_exp(t) = sum(raster_dotsOn(C,t,a:b)) / ((b-a+1)/1000);
                    spCount_exp(t) = sum(raster_dotsOn(C,t,a:b));
                end
            end % else, trial exceeds the max length specified in prealloc (will stay nans)
                        
            % raster, align RT (dotsOff)
            if round(trialLen(t)*1000)<size(raster_RT,3)
                raster_RT(C,t,end-round(trialLen(t)*1000):end) = 0; % init bins to zero but keep nans before dotsOn-extRaster(1)
                spkTimes = flipud(endTrial(t)-dataCell{n}.Exp.spikeTimes{c}(spkInd)); % flipud instead, makes hist easier
                if ~isempty(spkTimes) && spkTimes(end)>0.001
                    x = 0:0.001:ceil(spkTimes(end)*1000)/1000;
                    raster_RT(C,t,end-length(x)+2:end) = fliplr(histcounts(spkTimes,x));
                end            
            end            
        end
        
        % identify pref choice targ [for referencing CPs below]
        
        % two options here: fit of mapping data, or just compare spike
        % rates on R and L trials during task (try both, but I think
        % mapping data makes the most sense in terms of what the readout
        % 'knows'; task data can be biased by, well, biases in RT or conf)        
        pref = double(prefDir_fit(C)<90 || prefDir_fit(C)>270);
        % OR
%         I = dataCell{n}.Exp.openEvents.direction==0 & dataCell{n}.Exp.openEvents.coherence > 0.25;   % high coh but some sess might exclude 0.512, 
%         J = dataCell{n}.Exp.openEvents.direction==180 & dataCell{n}.Exp.openEvents.coherence > 0.25; % plus we get more s/n this way
%         pref = double(nanmean(spRate_exp(I)) > nanmean(spRate_exp(J)));

%************************************************************************
% CP, confP, and choice/conf-conditioned PSTHs
%************************************************************************
        
        % First, conventional CP from full-trial spike rates, referenced to
        % preferred direction (as defined by the task trials themselves),
        % rather than Haefner-style which is to reference to a fixed choice.
        % We use rates instead of counts because in RT task, using counts
        % would overweight long-RT trials. Also, the choice-conditioned
        % count distributions may not be equal under the null hypothesis,
        % if RT isn't perfectly symmetric for the two choices.
        
        % weak motion only, or else there's a confound (later we'll pool across cohs)
            %   NEW: also enforce a min RT
        J = abs(dataCell{n}.Exp.openEvents.coherence)<=maxCoh & RT' >= minRT;

        % there are 4 choices available, so CP and ConfP each come in 2 flavors
        I = J & dataCell{n}.Exp.openEvents.choice==pref & dataCell{n}.Exp.openEvents.pdw==1;
            X = spRate_exp(I); % pref-high        
        I = J & dataCell{n}.Exp.openEvents.choice~=pref & dataCell{n}.Exp.openEvents.pdw==1;
            Y = spRate_exp(I); % null-high
        CP_high(C) = rocN(X(~isnan(X)),Y(~isnan(Y)),100);
                
        I = J & dataCell{n}.Exp.openEvents.choice==pref & dataCell{n}.Exp.openEvents.pdw==0;
            X = spRate_exp(I); % pref-low        
        I = J & dataCell{n}.Exp.openEvents.choice~=pref & dataCell{n}.Exp.openEvents.pdw==0;
            Y = spRate_exp(I); % null-low
        CP_low(C) = rocN(X(~isnan(X)),Y(~isnan(Y)),100);
        
        I = J & dataCell{n}.Exp.openEvents.choice==pref;
            X = spRate_exp(I); % pref-all        
        I = J & dataCell{n}.Exp.openEvents.choice~=pref;
            Y = spRate_exp(I); % null-all
        CP_all(C) = rocN(X(~isnan(X)),Y(~isnan(Y)),100);
        
        I = J & dataCell{n}.Exp.openEvents.choice==pref & dataCell{n}.Exp.openEvents.pdw==1;
            X = spRate_exp(I); % pref-high        
        I = J & dataCell{n}.Exp.openEvents.choice==pref & dataCell{n}.Exp.openEvents.pdw==0;
            Y = spRate_exp(I); % pref-low
        ConfP_pref(C) = rocN(X(~isnan(X)),Y(~isnan(Y)),100);
        
        I = J & dataCell{n}.Exp.openEvents.choice~=pref & dataCell{n}.Exp.openEvents.pdw==1;
            X = spRate_exp(I); % null-high        
        I = J & dataCell{n}.Exp.openEvents.choice~=pref & dataCell{n}.Exp.openEvents.pdw==0;
            Y = spRate_exp(I); % null-low
        ConfP_null(C) = rocN(X(~isnan(X)),Y(~isnan(Y)),100);
        
        I = J & dataCell{n}.Exp.openEvents.pdw==1;
            X = spRate_exp(I); % all-high        
        I = J & dataCell{n}.Exp.openEvents.pdw==0;
            Y = spRate_exp(I); % all-low
        ConfP_all(C) = rocN(X(~isnan(X)),Y(~isnan(Y)),100);
                
% % grab an example pair of dists, eg for talk figure
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
%     export_fig('ConfP_example_dists', '-eps');
% end
        



%         % Next, repeat the above but separately for each coh (within each
%         % cell), then take a weighted average across most cohs. This was an
%         % original strategy of Britten et al: even with a strong sensory
%         % signal, the choice-conditioned distributions (within a coh) are
%         % overlapping under the null hypothesis. Of course, usually cannot
%         % use the 1-2 highest cohs because not enough error trials.
%         
%         keyboard
%         dataCell{n}.Exp.openEvents.choice
%         
%         for n = 1
%         
%             
%         end



        % PSTHs
        
        % aligned dots-on
        clear psth;
        I = J & dataCell{n}.Exp.openEvents.choice==pref & dataCell{n}.Exp.openEvents.pdw==1; % pref-high
        psth(1,:) = smoothRaster(nanmean(squeeze(raster_dotsOn(C,I,:)))*1e3, convKernel);
        I = J & dataCell{n}.Exp.openEvents.choice==pref & dataCell{n}.Exp.openEvents.pdw==0; % pref-low
        psth(2,:) = smoothRaster(nanmean(squeeze(raster_dotsOn(C,I,:)))*1e3, convKernel);
        I = J & dataCell{n}.Exp.openEvents.choice~=pref & dataCell{n}.Exp.openEvents.pdw==1; % null-high
        psth(3,:) = smoothRaster(nanmean(squeeze(raster_dotsOn(C,I,:)))*1e3, convKernel);
        I = J & dataCell{n}.Exp.openEvents.choice~=pref & dataCell{n}.Exp.openEvents.pdw==0; % null-low
        psth(4,:) = smoothRaster(nanmean(squeeze(raster_dotsOn(C,I,:)))*1e3, convKernel);
        rMax = max(max(psth(:,1:300)));
        psthNorm_dotsOn_one(C,:) = psth(1,:)/rMax;
        psthNorm_dotsOn_two(C,:) = psth(2,:)/rMax;
        psthNorm_dotsOn_three(C,:) = psth(3,:)/rMax;
        psthNorm_dotsOn_four(C,:) = psth(4,:)/rMax;

        
        
        
        
        if plotEverythingElse
            figure(n*10); subplot(1,nU(n),c);
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
            xlabel('Time from dots on (ms)');
            ylabel('Firing rate (spikes/s)');
            changeAxesFontSize(gca, 20, 20);
            if c==nU(n)
                legend('pref high','pref low','null high','null low','Location','Northeast'); legend('boxoff');
            end
        end

        % aligned RT
        clear psth
        I = J & dataCell{n}.Exp.openEvents.choice==pref & dataCell{n}.Exp.openEvents.pdw==1; % pref-high
        psth(1,:) = smoothRaster(nanmean(squeeze(raster_RT(C,I,:)))*1e3, convKernel);
        I = J & dataCell{n}.Exp.openEvents.choice==pref & dataCell{n}.Exp.openEvents.pdw==0; % pref-low
        psth(2,:) = smoothRaster(nanmean(squeeze(raster_RT(C,I,:)))*1e3, convKernel);
        I = J & dataCell{n}.Exp.openEvents.choice~=pref & dataCell{n}.Exp.openEvents.pdw==1; % null-high
        psth(3,:) = smoothRaster(nanmean(squeeze(raster_RT(C,I,:)))*1e3, convKernel);
        I = J & dataCell{n}.Exp.openEvents.choice~=pref & dataCell{n}.Exp.openEvents.pdw==0; % null-low
        psth(4,:) = smoothRaster(nanmean(squeeze(raster_RT(C,I,:)))*1e3, convKernel);
        psthNorm_RT_one(C,:) = psth(1,:)/rMax;
        psthNorm_RT_two(C,:) = psth(2,:)/rMax;
        psthNorm_RT_three(C,:) = psth(3,:)/rMax;
        psthNorm_RT_four(C,:) = psth(4,:)/rMax;

        if plotEverythingElse
            figure(n*100); subplot(1,nU(n),c);
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
        
        sessID(C) = n;
        C = C+1;
        
    end
    
    if plotEverythingElse; pause; end
end



%% save

clear dataCell raster_dotsOn raster_RT raster_map
% save neuralAnalysis.mat
save neuralAnalysis_minRT500.mat


%% reload

clear
% load neuralAnalysis.mat
load neuralAnalysis_minRT500.mat


% randos

figure;hist(prefDir_fit,18);
xlabel('Preferred direction (deg)');
ylabel('Number of units');
set(gca,'xlim',[0 360],'xtick',0:60:360);
changeAxesFontSize(gca,16,16);



%% plot average choice/conf-conditioned PSTH

    % aligned dotsOn
clear psth
psth(1,:) = nanmean(psthNorm_dotsOn_one);
psth(2,:) = nanmean(psthNorm_dotsOn_two);
psth(3,:) = nanmean(psthNorm_dotsOn_three);
psth(4,:) = nanmean(psthNorm_dotsOn_four);

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
xlabel('Time from dots on (ms)');
ylabel('Firing rate (normalized)');
changeAxesFontSize(gca, 22, 22);
h = legend('pref choice, high bet','pref choice, low bet','null choice, high bet','null choice, low bet','Location','Northeast'); legend('boxoff');
export_fig('psthNorm_dotsOn_coh0', '-eps');
export_fig('psthNorm_dotsOn_coh32', '-eps');

    % aligned RT
clear psth
psth(1,:) = nanmean(psthNorm_RT_one);
psth(2,:) = nanmean(psthNorm_RT_two);
psth(3,:) = nanmean(psthNorm_RT_three);
psth(4,:) = nanmean(psthNorm_RT_four);

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
export_fig('psthNorm_RT_coh0', '-eps');
export_fig('psthNorm_RT_coh0_extended', '-eps');
export_fig('psthNorm_RT_coh32', '-eps');
export_fig('psthNorm_RT_coh32_extended', '-eps');


%% CP and ConfP, vs each other and as a func of prefDir

nanmean(CP_high(:))
nanmean(CP_low(:))
nanmean(CP_all(:))
nanmean(ConfP_pref(:))
nanmean(ConfP_null(:))
nanmean(ConfP_all(:))


%% choose

% CP = CP_high;
% or
CP = CP_low;
% or
% CP = CP_all;

% ConfP = ConfP_pref;
% or
ConfP = ConfP_null;
% or
% ConfP = ConfP_all;

nans = isnan(CP) | isnan(ConfP);

CP(nans) = [];
ConfP(nans) = [];

CPse = std(CP)/sqrt(length(CP));
[~,P] = ttest(CP,0.5)

ConfPse = std(ConfP)/sqrt(length(ConfP));
[~,P] = ttest(ConfP,0.5)



%%
figure; set(gcf, 'Color', 'w', 'Position', [400 200 450 450], 'PaperPositionMode', 'auto');
plot([0 1],[0 1],'k--',[0 1],[0.5 0.5],'k:',[0.5 0.5],[0 1],'k:');
hold on; plot(CP,ConfP,'ko','MarkerSize',11,'MarkerFaceColor',[1 1 1]); axis square;
% set(gca,'XLim',[0.17 0.83],'YLim',[0.17 0.83],'xtick',0.2:0.1:0.8,'ytick',0.2:0.1:0.8,'TickDir','out')
set(gca,'XLim',[0.1 0.9],'YLim',[0.1 0.9],'xtick',0.1:0.2:0.9,'ytick',0.1:0.2:0.9,'TickDir','out')
xlabel('Choice probability'); ylabel('Confidence probability');
changeAxesFontSize(gca, 22, 22);
export_fig('CPvsConfP_coh0', '-eps');

figure; set(gcf, 'Color', 'w', 'Position', [800 200 450 80], 'PaperPositionMode', 'auto');
[histCP,E] = histcounts(CP,0:0.05:1);
bar(E(2:end)-0.025,histCP); box('off');
set(gca,'XLim',[0.1 0.9],'xtick',0.1:0.2:0.9,'xticklabel',[],'TickDir','out');
changeAxesFontSize(gca, 14, 14);
export_fig('CPbar_coh0', '-eps');

figure; set(gcf, 'Color', 'w', 'Position', [800 100 80 450], 'PaperPositionMode', 'auto');
[histCP,E] = histcounts(ConfP,0:0.05:1);
barh(E(2:end)-0.025,histCP); box('off');
set(gca,'YLim',[0.1 0.9],'ytick',0.1:0.2:0.9,'yticklabel',[],'TickDir','out');
changeAxesFontSize(gca, 14, 14);
export_fig('ConfPbar_coh0', '-eps');

[r,p] = corrcoef(CP,ConfP)


%%
if any(isnan(prefDir_fit)); keyboard; end
prefDir = prefDir_fit(~nans);

figure; set(gcf, 'Color', 'w', 'Position', [400 200 450 450], 'PaperPositionMode', 'auto');
plot(prefDir,CP,'ko','MarkerSize',11,'MarkerFaceColor',[1 1 1]);
set(gca,'xlim',[0 360],'xtick',0:45:360,'tickDir','out')
xlabel('prefDir'); ylabel('CP');

figure; set(gcf, 'Color', 'w', 'Position', [400 200 450 450], 'PaperPositionMode', 'auto');
plot(prefDir,ConfP,'ko','MarkerSize',11,'MarkerFaceColor',[1 1 1]);
set(gca,'xlim',[0 360],'xtick',0:45:360,'tickDir','out')
xlabel('prefDir'); ylabel('ConfP');

I = cosd(prefDir)>sqrt(2)/2 | cosd(prefDir)<-sqrt(2)/2; % within +/- 45 deg of 0 or 180
CP_within45 = mean(CP(I))
CP_outside45 = mean(CP(~I))
% whaaaat!


%%
% running mean

dbstop if error

windowWidth = 60;

figure; set(gcf, 'Color', 'w', 'Position', [400 200 450 450], 'PaperPositionMode', 'auto');

[x_,y_,~] = running_mean_circular(prefDir, CP, windowWidth, 'window_width');
% [x_,y_,~] = running_mean_circular(prefDir, CP, windowWidth, 'data_point');
h(1) = plot(x_,y_,'k-','LineWidth',2); hold on;

[x_,y_,~] = running_mean_circular(prefDir, ConfP, windowWidth, 'window_width');
% [x_,y_,~] = running_mean_circular(prefDir, ConfP, windowWidth, 'data_point');
h(2) = plot(x_,y_,'r-','LineWidth',2);



% SWITCH TO POLAR PLOTS?



plot([0 360],[0.5 0.5],'k--');
set(gca,'xlim',[0 360],'xtick',0:45:360,'ylim',[0.35 0.65], 'ytick',0.4:0.05:0.6,'tickDir','out')
xlabel('Preferred direction (deg)'); ylabel('Choice/Confidence Probability');
legend(h(1:2),'CP','ConfP','location','northwest'); legend('boxoff');
changeAxesFontSize(gca, 16, 16);

% export_fig('CP-v-prefdir_runningMean', '-eps');


% % sanity check for running_mean_circular
% x = unifrnd(0,359.99,10000,1);
% y = normrnd(1,0.2,size(x));
% [x_,y_,~] = running_mean_circular(x, y, 60, 'window_width');
% h(1) = plot(x_,y_,'k-','LineWidth',2); hold on;
% ylim([0.5 1.5]);


%% repeat CP using residuals (grand-CP)

% [actually make that a toggle above]




%% need to validate CPs by checking RT dists of the two groups being compared,
%% and trying a cap on the time window,
%% and by using spike count instead of rate (just curious)
%% and by inspecting a sample of trials for anomalies in the FR dists (bimodality?)
%% and by trying different parameters of ROCn







%% convert this to readout weights using Haefner equation?
% [waiting for corr info]





%% split CP by time; first just by short/long RT,
% then within long-RT only, by early vs late

























