function dots3DMP_plots_fit_byConf(data,fit,conftask,RTtask)
% similar to dots3DMP_plots_fit_byCoh, but for different confidence levels

% calls parseData on struct data and plots just the data points, then does
% the analogous step for struct fit, and plots the result as lines (interp)
%   [CF removed fitGauss as it should be obsolete with non-MC]

if nargin<3, conftask=2; end
if nargin<4, RTtask = 0; end

fsz = 14; % fontsize

mods   = unique(data.modality);
cohs   = unique(data.coherence);
deltas = unique(data.delta);
hdgs   = unique(data.heading);

if ~isfield(data,'oneTargConf')
    data.oneTargConf = false(size(data.heading));
end

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
RTCorrOnly = 1; % when fitting correct RTs only, only compare to those in the data

parsedData = dots3DMP_parseData_byConf(data,mods,cohs,deltas,hdgs,conftask,RTtask);

if useAbsHdg
    hdgs = unique(abs(hdgs));
end

spRows = 1 + double(RTtask);

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
        h(m,1) = plot(hdgs, squeeze(parsedData.pRight(m,c,D,:,1)), clr{c}{m},'linewidth',1.5,'linestyle','none');
        h(m,2) = plot(hdgs, squeeze(parsedData.pRight(m,c,D,:,2)), clr{c}{m}(1),'linewidth',1.5,'linestyle','none','marker','x');
    end
    
    if RTtask
        subplot(spRows,length(cohs),c+length(cohs)); hold on
        for m = 1:length(mods)
%             h(m) = errorbar(hdgs, squeeze(parsedData.RTmean(m,c,D,:)), squeeze(parsedData.RTse(m,c,D,:)), [clr{c}{m}],'linewidth',1.5); 


            plot(hdgs, squeeze(parsedData.RTmean(m,c,D,:,1)), clr{c}{m}(1),'linestyle','none','linewidth',1.5,'marker',clr{c}{m}(2));
            plot(hdgs, squeeze(parsedData.RTmean(m,c,D,:,2)), clr{c}{m}(1),'linestyle','none','linewidth',1.5,'marker','x');
        
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

try 
mods   = unique(fit.modality);
cohs   = unique(fit.coherence);
deltas = unique(fit.delta);
% deltas = 0;
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
        
        K = J & fit.PDW==1;
        n(m,c,d,h,1) = nansum(K);
        pRightHigh(m,c,d,h) = nanmean(fit.pRight(K));
        pRightHighSE(m,c,d,h) = nanstd(fit.pRight(K))/sqrt(n(m,c,d,h,1));

        K = J & fit.PDW==0;
        n(m,c,d,h,2) = nansum(K);
        pRightLow(m,c,d,h) = nanmean(fit.pRight(K));
        pRightLowSE(m,c,d,h) = nanstd(fit.pRight(K))/sqrt(n(m,c,d,h,2));
        
        if RTtask
            K = J & fit.PDW==1;
            RTmeanHigh(m,c,d,h) = nanmean(fit.RT(K));
            RTHighse(m,c,d,h) = nanstd(fit.RT(J))/sqrt(n(m,c,d,h,1));

            K = J & fit.PDW==0;
            RTmeanLow(m,c,d,h) = nanmean(fit.RT(K));
            RTLowse(m,c,d,h) = nanstd(fit.RT(J))/sqrt(n(m,c,d,h,2));
        else
            RTmeanHigh(m,c,d,h) = NaN;
            RTHighse(m,c,d,h) = NaN;
        end        
    end
end
end
end

% copy vestib-only data to all coherences, to aid plotting
for c=1:length(cohs)
    n(1,c,:,:,:) = n(1,1,:,:,:);
    pRightHigh(1,c,:,:) = pRightHigh(1,1,:,:);
    pRightHighSE(1,c,:,:) = pRightHighSE(1,1,:,:);
    RTmeanHigh(1,c,:,:) = RTmeanHigh(1,1,:,:);
    RTHighse(1,c,:,:) = RTHighse(1,1,:,:);

    pRightLow(1,c,:,:) = pRightLow(1,1,:,:);
    pRightLowSE(1,c,:,:) = pRightLowSE(1,1,:,:);
    RTmeanLow(1,c,:,:) = RTmeanLow(1,1,:,:);
    RTLowse(1,c,:,:) = RTLowse(1,1,:,:);
end

% make analogous 'parsed' struct
parsedFit = struct();
parsedFit.pRightHigh = pRightHigh;
parsedFit.pRightLow = pRightLow;
if RTtask
    parsedFit.RTmeanHigh = RTmeanHigh;
    parsedFit.RTmeanLow = RTmeanLow;
end

catch
    disp('using passed parsedFit')
    parsedFit = fit;
end


% plot it!

% first, single-cues + zero-delta comb

%ves %vis %comb
% clr{1} = {'k-','m-','c-'};
clr{1} = {'k-','r-','b-'};
clr{2} = {'k-','r-','b-'};
clr{3} = {'k-','y-','g-'};

hdgs = linspace(-12,12,33);
figure(101);
for c = 1:length(cohs)
    subplot(spRows,length(cohs),c);
    for m = 1:length(mods)     % m c d h
        h(m,1) = plot(hdgs, squeeze(parsedFit.pRightHigh(m,c,D,:)), [clr{c}{m}]); hold on;
        h(m,2) = plot(hdgs, squeeze(parsedFit.pRightLow(m,c,D,:)), 'color',[clr{c}{m}(1)] ,'Linestyle','--'); hold on;

        ylim([0 1]);
    end
    if RTtask
        subplot(spRows,length(cohs),c+length(cohs));
        for m = 1:length(mods)
            h(m) = plot(hdgs, squeeze(parsedFit.RTmeanHigh(m,c,D,:)), [clr{c}{m}]); hold on;
            h(m) = plot(hdgs, squeeze(parsedFit.RTmeanLow(m,c,D,:)), 'color',[clr{c}{m}(1)],'linestyle','--'); hold on;
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


