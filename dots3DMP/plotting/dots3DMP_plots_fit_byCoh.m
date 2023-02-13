function dots3DMP_plots_fit_byCoh(data,fit,conftask,RTtask)

% somewhat unwieldy amalgam of parseData and dots3DMP_plots, to show model 
% fits versus data in dots3DMP experiment

% calls parseData on struct data and plots just the data points, then does
% the analogous step for struct fit, and plots the result as lines (interp)
%   [CF removed fitGauss as it should be obsolete with non-MC]

if nargin<3, conftask=1; end
if nargin<4, RTtask = 0; end

fsz = 14; % fontsize

mods   = unique(data.modality);
cohs   = unique(data.coherence);
deltas = unique(data.delta);
hdgs   = unique(data.heading);

if all(mods==1), cohs=1; end

xLab = sprintf('heading angle (%s)',char(176));
modlabels = {'Ves','Vis','Comb'};

if conftask==1 
    yLab = 'SaccEP'; confYlims = [0.2 0.9]; 
    if all(mods==1), RTylims = [0.9 1.5]; 
    else,            RTylims = [0.9 1.9]; 
    end
    xt = -10:5:10;
     if length(mods)>1, cohlabs = {'Low Coh','High Coh'}; end
elseif conftask==2
    yLab = 'P(High Bet)'; confYlims = [0.4 1.0]; 
    if all(mods==1), RTylims = [0.5 0.72]; 
    else,            RTylims = [0.5 0.9]; 
    end
    xt = -12:6:12;
    if length(mods)>1, cohlabs = {['coh = ' num2str(cohs(1))],['coh = ' num2str(cohs(2))]}; end
end 

%% actual data
% plot as individual points with error bars

useAbsHdg = 0;
RTCorrOnly = 0; % when fitting correct RTs only, only compare to those in the data

parsedData = dots3DMP_parseData(data,mods,cohs,deltas,hdgs,conftask,RTtask,useAbsHdg,RTCorrOnly); 

if useAbsHdg
    hdgs = unique(abs(hdgs));
end

spRows = 1 + double(conftask>0) + double(RTtask);

if conftask==1, confYL = 'Sacc EP';
elseif conftask==2, confYL = 'P(High Bet)';
end

% first, single-cues + zero-delta comb
D = find(deltas==0);
         %ves %vis %comb
% clr{1} = {'ko','mo','co'};
clr{1} = {'ko','ro','bo'};
clr{2} = {'ko','ro','bo'};
clr{3} = {'ko','yo','go'};
figure(101);
set(gcf,'Color',[1 1 1],'Position',[300 1000 230+300*(length(cohs)-1) 200+150*(conftask>0)+150*RTtask],'PaperPositionMode','auto'); clf;
for c = 1:length(cohs)
    subplot(spRows,length(cohs),c); hold on;
    for m = 1:length(mods)     % m c d h
        h(m) = errorbar(hdgs, squeeze(parsedData.pRight(m,c,D,:)), squeeze(parsedData.pRightSE(m,c,D,:)), [clr{c}{m}],'linewidth',1.5); 
    end
    
    if conftask>0
        subplot(spRows,length(cohs),c+length(cohs)); hold on
        for m = 1:length(mods)
            h(m) = errorbar(hdgs, squeeze(parsedData.confMean(m,c,D,:)), squeeze(parsedData.confSE(m,c,D,:)), [clr{c}{m}],'linewidth',1.5);
        end

    end
    
    if RTtask
        subplot(spRows,length(cohs),c+length(cohs)*(2-(conftask==0))); hold on
        for m = 1:length(mods)
            h(m) = errorbar(hdgs, squeeze(parsedData.RTmean(m,c,D,:)), squeeze(parsedData.RTse(m,c,D,:)), [clr{c}{m}],'linewidth',1.5); 
        end
    end
    
end

% now comb separated by delta
if length(deltas)>1

clr{1} = {'bs','cs','gs'};
clr{2} = {'b^','c^','g^'};
clr{3} = {'bo','co','go'};

clear L;
figure(108);
set(gcf,'Color',[1 1 1],'Position',[50 20 230+300*(length(cohs)-1) 200+150*(conftask>0)+150*RTtask],'PaperPositionMode','auto'); clf;
for c = 1:length(cohs)
    subplot(spRows,length(cohs),c); hold on
    for d = 1:length(deltas)     % m c d h
        h(d) = errorbar(hdgs, squeeze(parsedData.pRight(3,c,d,:)), squeeze(parsedData.pRightSE(3,c,d,:)), [clr{c}{d}],'linewidth',1.5);
%         L{d} = sprintf('\x0394=%d',deltas(d));
    end

    if conftask>0
        subplot(spRows,length(cohs),c+length(cohs)); hold on
        for d = 1:length(deltas)
            h(d) = errorbar(hdgs, squeeze(parsedData.confMean(3,c,d,:)), squeeze(parsedData.confSE(3,c,d,:)), [clr{c}{d}],'linewidth',1.5);
        end
    end
    
    if RTtask
        subplot(spRows,length(cohs),c+length(cohs)*(2-(conftask==0))); hold on
        for d = 1:length(deltas)
            h(d) = errorbar(hdgs, squeeze(parsedData.RTmean(3,c,d,:)), squeeze(parsedData.RTse(3,c,d,:)), [clr{c}{d}],'linewidth',1.5); 
        end
    end
    
end

end


%% fit
% plot curve only, overlaid on actual data points

mods   = unique(fit.modality);
cohs   = unique(fit.coherence);
deltas = unique(fit.delta);
hdgs   = unique(fit.heading);
D = find(deltas==0);

%%% parsedData = dots3DMP_parseData(fit,mods,cohs,deltas,hdgs,conftask,RTtask); 
% ^ because fit may have continuous pRight/pHigh instead of binary
% choice/PDW, use this abridged version of parseData instead:
n = nan(length(mods),length(cohs),length(deltas),length(hdgs));
pRight = n;
pRightSE = n;
RTmean = n; RTse = n;
confMean = n; confSE = n;
for m = 1:length(mods)
for c = 1:length(cohs)
for d = 1:length(deltas)     
    for h = 1:length(hdgs)
        J = fit.modality==mods(m) & fit.coherence==cohs(c) & fit.heading==hdgs(h) & fit.delta==deltas(d);
        n(m,c,d,h) = nansum(J);
        pRight(m,c,d,h) = nanmean(fit.pRight(J));
        pRightSE(m,c,d,h) = nanstd(fit.pRight(J))/sqrt(n(m,c,d,h));
        if RTtask
            RTmean(m,c,d,h) = nanmean(fit.RT(J));
            RTse(m,c,d,h) = nanstd(fit.RT(J))/sqrt(n(m,c,d,h));
        else
            RTmean(m,c,d,h) = NaN;
            RTse(m,c,d,h) = NaN;
        end        
        if conftask
            confMean(m,c,d,h) = nanmean(fit.conf(J));
            confSE(m,c,d,h) = nanstd(fit.conf(J))/sqrt(n(m,c,d,h));        
        else
            confMean(m,c,d,h) = NaN;
            confSE(m,c,d,h) = NaN;
        end
    end
end
end
end

% copy vestib-only data to all coherences, to aid plotting
for c=1:length(cohs)
    n(1,c,:,:) = n(1,1,:,:);
    pRight(1,c,:,:) = pRight(1,1,:,:);
    pRightSE(1,c,:,:) = pRightSE(1,1,:,:);
    confMean(1,c,:,:) = confMean(1,1,:,:);
    confSE(1,c,:,:) = confSE(1,1,:,:);
    RTmean(1,c,:,:) = RTmean(1,1,:,:);
    RTse(1,c,:,:) = RTse(1,1,:,:);
end

% make analogous 'parsed' struct
parsedFit = struct();
parsedFit.pRight = pRight;
if conftask
    parsedFit.confMean = confMean;
end
if RTtask
    parsedFit.RTmean = RTmean;
    end


% plot it!

% first, single-cues + zero-delta comb

%ves %vis %comb
% clr{1} = {'k-','m-','c-'};
clr{1} = {'k-','r-','b-'};
clr{2} = {'k-','r-','b-'};
clr{3} = {'k-','y-','g-'};

figure(101);
for c = 1:length(cohs)
    subplot(spRows,length(cohs),c);
    for m = 1:length(mods)     % m c d h
        h(m) = plot(hdgs, squeeze(parsedFit.pRight(m,c,D,:)), [clr{c}{m}]); hold on;
        ylim([0 1]);
    end
    if conftask
        subplot(spRows,length(cohs),c+length(cohs));
        for m = 1:length(mods)
            h(m) = plot(hdgs, squeeze(parsedFit.confMean(m,c,D,:)), [clr{c}{m}]);
            ylim([0 1]); hold on;
        end
    end
    if RTtask
        subplot(spRows,length(cohs),c+length(cohs)*2);
        for m = 1:length(mods)
            h(m) = plot(hdgs, squeeze(parsedFit.RTmean(m,c,D,:)), [clr{c}{m}]); hold on;
        end
    end
end

% now comb separated by delta

clr{1} = {'b-','c-','g-'};
clr{2} = {'b-','c-','g-'};
clr{3} = {'b-','c-','g-'};

if length(deltas)>1
    figure(108);
    for c = 1:length(cohs)
        subplot(spRows,length(cohs),c); hold on
        for d = 1:length(deltas)     % m c d h
            h(d) = plot(hdgs, squeeze(parsedFit.pRight(3,c,d,:)), [clr{c}{d}]); 
        end
        if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
        ylim([0 1]);
        if conftask>0
            subplot(spRows,length(cohs),c+length(cohs)); hold on
            for d = 1:length(deltas)
                h(d) = plot(hdgs, squeeze(parsedFit.confMean(3,c,d,:)), [clr{c}{d}]);
            end
        end
        if RTtask
            subplot(spRows,length(cohs),c+length(cohs)*(2-(conftask==0))); hold on;
            for d = 1:length(deltas)
                h(d) = plot(hdgs, squeeze(parsedFit.RTmean(3,c,d,:)), [clr{c}{d}]);
            end
        end
    end        
end


