function fh = dots3DMP_plots(parsedData,mods,cohs,deltas,hdgs,conftask,RTtask,D)

if nargin < 8 || isempty(D)
    % first, for all trials irrespective of delta
    D = length(deltas)+1; % (the extra column we made for pooling across deltas)
    % OR select just delta=0:
%     D = find(deltas==0);
end

         %ves %vis %comb
% clr{1} = {'ko','mo','co'};
clr{1} = {'ko','ro','bo'};
clr{2} = {'ko','ro','bo'};
clr{3} = {'ko','yo','go'};
lh = {'ves','vis','comb'};
fh(1) = figure(101+D);
% set(gcf,'Color',[1 1 1],'Position',[300 1000 450+300*(length(cohs)-2) 200+150*(conftask>0)+150*RTtask],'PaperPositionMode','auto'); clf;
set(gcf,'Color',[1 1 1],'Position',[200 80 700 900],'PaperPositionMode','auto'); clf;
for c = 1:length(cohs)
    
    % CHOICE
    subplot(1+double(conftask>0)+double(RTtask),length(cohs),c); box off; hold on;
    for m = 1:length(mods)     % m c d h
        if parsedData.plotLogistic(m,c,D)
            h(m) = plot(parsedData.xVals,squeeze(parsedData.yVals(m,c,D,:)),[clr{c}{mods(m)}(1) '-'],'linewidth',1.5); hold on;
            errorbar(hdgs, squeeze(parsedData.pRight(m,c,D,:)), squeeze(parsedData.pRightSE(m,c,D,:)), clr{c}{mods(m)},'linewidth',1.5);
        else
            h(m) = errorbar(hdgs, squeeze(parsedData.pRight(m,c,D,:)), squeeze(parsedData.pRightSE(m,c,D,:)), [clr{c}{mods(m)} '-'],'linewidth',1.5); hold on;
        end
        ylim([0 1]);
        if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
    end
    legend(h,lh(mods),'Location','northwest');
    xlabel('heading angle (deg)'); ylabel('P(right)');
    try changeAxesFontSize(gca,20,20); catch; end
    
    % CONFIDENCE
    if conftask
        subplot(1+double(conftask>0)+double(RTtask),length(cohs),c+length(cohs)); box off; hold on;
        for m = 1:length(mods)
            h(m) = errorbar(hdgs, squeeze(parsedData.confMean(m,c,D,:)), squeeze(parsedData.confSE(m,c,D,:)), [clr{c}{mods(m)} '-'],'linewidth',1.5);
            ylim([0.4 0.8]); hold on;
        end
        xlabel('heading angle (deg)');
         if conftask==1, ylim([0.4 0.8]); ylabel('SEP (''confidence'', %)');
        elseif conftask==2, ylim([0.4 0.1]); ylabel('P(high bet)');
        end
        try changeAxesFontSize(gca,20,20); catch; end
    end
    
    % RT
    if RTtask
        subplot(1+double(conftask>0)+double(RTtask),length(cohs),c+length(cohs)*2); box off; hold on;
        for m = 1:length(mods)
            h(m) = errorbar(hdgs, squeeze(parsedData.RTmean(m,c,D,:)), squeeze(parsedData.RTse(m,c,D,:)), [clr{c}{mods(m)} '-'],'linewidth',1.5); hold on;
        end
        xlabel('heading angle (deg)'); ylabel('RT (s)');
        yRng = [min(parsedData.RTmean(:)) max(parsedData.RTmean(:))];
        ylim(yRng.*[0.9 1.1])
    end
    try changeAxesFontSize(gca,20,20); catch; end
end


%% now separate by delta

if length(deltas)>1

clr{1} = {'bs','cs','gs'};
clr{2} = {'b^','c^','g^'};
clr{3} = {'bo','co','go'};

clear L;
fh(2) = figure(108);
% set(gcf,'Color',[1 1 1],'Position',[50 20 450+300*(length(cohs)-2) 200+150*(conftask>0)+150*RTtask],'PaperPositionMode','auto'); clf;
set(gcf,'Color',[1 1 1],'Position',[900 80 700 900],'PaperPositionMode','auto'); clf;
for c = 1:length(cohs)
    subplot(1+double(conftask>0)+RTtask,length(cohs),c); box off; hold on;
    for d = 1:length(deltas)     % m c d h
        if parsedData.plotLogistic(3,c,d)
            h(d) = plot(parsedData.xVals,squeeze(parsedData.yVals(3,c,d,:)),[clr{c}{d}(1) '-'], 'linewidth', 1.5); hold on;
            errorbar(hdgs, squeeze(parsedData.pRight(3,c,d,:)), squeeze(parsedData.pRightSE(3,c,d,:)), clr{c}{d});
        else
            h(d) = errorbar(hdgs, squeeze(parsedData.pRight(3,c,d,:)), squeeze(parsedData.pRightSE(3,c,d,:)), [clr{c}{d} '-'], 'linewidth', 1.5); hold on
        end
        L{d} = sprintf('\x0394=%d',deltas(d));
        ylim([0 1]);
        if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
    end
    legend(h,L,'location','northwest');
    xlabel('heading angle (deg)'); ylabel('P(right)');
    try changeAxesFontSize(gca,20,20); catch; end
    
    if conftask
        subplot(1+double(conftask>0)+double(RTtask),length(cohs),c+length(cohs)); box off; hold on;
        for d = 1:length(deltas)
            h(d) = errorbar(hdgs, squeeze(parsedData.confMean(3,c,d,:)), squeeze(parsedData.confSE(3,c,d,:)), [clr{c}{d} '-'], 'linewidth', 1.5);
            hold on;
        end
        
        xlabel('heading angle (deg)'); 
        if conftask==1, ylim([0.5 0.8]); ylabel('SEP (''confidence'', %)');
        elseif conftask==2, ylim([0.4 0.1]); ylabel('P(high bet)');
        end
        try changeAxesFontSize(gca,20,20); catch; end
    end 
    
    if RTtask
        subplot(1+double(conftask>0)+double(RTtask),length(cohs),c+length(cohs)*2); box off; hold on;
        for d = 1:length(deltas)
            h(d) = errorbar(hdgs, squeeze(parsedData.RTmean(3,c,d,:)), squeeze(parsedData.RTse(3,c,d,:)), [clr{c}{d} '-'], 'linewidth', 1.5); hold on;
        end
        xlabel('heading angle (deg)'); ylabel('RT (s)');
        yRng = [min(parsedData.RTmean(:)) max(parsedData.RTmean(:))];
        ylim(yRng.*[0.9 1.1])
    end
    try changeAxesFontSize(gca,20,20); catch; end
end

end