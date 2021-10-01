 sigmaPMF(3,1,3) = sigmaPMF(3,1,3)*.7;
 sigmaPMF(3,2,1) = sigmaPMF(3,2,1)*.6;
 muPMF(3,2,1) = muPMF(3,2,1)-0.2;

%% first, for all trials irrespective of delta
D = length(deltas)+1; % (the extra column we made for pooling across deltas)
% OR select just delta=0:
% D = find(deltas==0);

modlabels = {'Ves','Vis','Comb'};

         %ves %vis %comb
clr{1} = {'ko','mo','co'};
clr{2} = {'ko','ro','bo'};
clr{3} = {'ko','yo','go'};
figure(201+D);
set(gcf,'Color',[1 1 1],'Position',[300 500 950+300*(length(cohs)-2) 800],'PaperPositionMode','auto'); clf;
for c = 1:length(cohs)
    % choice
    subplot(1+double(conftask>0)+double(RTtask),length(cohs),c); box off; hold on;
    for m = 1:length(mods)     % m c d h
        beta = [muPMF(m,c,D) sigmaPMF(m,c,D)];
        h(m) = plot(xVals, cgauss(beta,xVals), [clr{c}{m}(1) '-'],'linewidth',1.5); hold on;
        errorbar(hdgs, squeeze(pRight(m,c,D,:)), squeeze(pRightSE(m,c,D,:)), clr{c}{m},'linewidth',1.5);
        ylim([0 1]); if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
        text(hdgs(1)+0.5,1.0-m*0.07,sprintf('%s: mu = %.2f, s = %.2f',modlabels{m},beta(1),beta(2)),'color',clr{c}{m}(1))
    end
%     legend(h,'vestib','visual','comb','Location','northwest');
    
    xlabel('heading angle (deg)'); ylabel('P(Right)');
    try changeAxesFontSize(gca,15,15); catch; end

    % conf
    if conftask
        subplot(2+double(RTtask),length(cohs),c+length(cohs)); box off; hold on;
        for m = 1:length(mods)
            beta = [amplConf(m,c,D) muConf(m,c,D) sigmaConf(m,c,D) baselineConf(m,c,D)];
            h(m) = plot(xVals, flippedGauss(beta,xVals), [clr{c}{m}(1) '-'],'linewidth',1.5); hold on;
            errorbar(hdgs, squeeze(confMean(m,c,D,:)), squeeze(confSE(m,c,D,:)), clr{c}{m},'linewidth',1.5);
            ylim([0 max(max(confMean(:)),1)]); if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
        end
        %     legend(h,'vestib','visual','comb','location','northwest');
        xlabel('heading angle (deg)');
        if conftask==1, ylabel('SEP (''confidence'', %)');
        elseif conftask==2, ylabel('proportion high bet');
        end
        try changeAxesFontSize(gca,15,15); catch; end
    end
    
    % RT
    if RTtask
        subplot(2+double(conftask>0),length(cohs),c+length(cohs)*2); box off; hold on;
        for m = 1:length(mods)        
            beta = [amplRT(m,c,D) muRT(m,c,D) sigmaRT(m,c,D) baselineRT(m,c,D)];
            h(m) = plot(xVals, gauss(beta,xVals), [clr{c}{m}(1) '-'],'linewidth',1.5); hold on;       
            errorbar(hdgs, squeeze(RTmean(m,c,D,:)), squeeze(RTse(m,c,D,:)), clr{c}{m},'linewidth',1.5);
            if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
        end
        if strcmp(subject,'lucio')
            ylim([0.5 0.9])
        else
            ylim([0.6 1.8]);
        end
        ylabel('RT (s)')
        try changeAxesFontSize(gca,15,15); catch; end


    end
end


%% now separate by delta

if length(deltas)>1
    
clr{1} = {'bs','cs','gs'};
clr{2} = {'b^','c^','g^'};
clr{3} = {'bo','co','go'};

clear L;
figure(208);
set(gcf,'Color',[1 1 1],'Position',[50 20 950+300*(length(cohs)-2) 800],'PaperPositionMode','auto'); clf;
for c = 1:length(cohs)
    % choice
    subplot(1+double(conftask>0)+double(RTtask),length(cohs),c); box off; hold on;
    for d = 1:length(deltas)     % m c d h
        beta = [muPMF(3,c,d) sigmaPMF(3,c,d)];
        h(d) = plot(xVals, cgauss(beta,xVals), [clr{c}{d}(1) '-'],'linewidth',1.5); hold on;
        errorbar(hdgs, squeeze(pRight(3,c,d,:)), squeeze(pRightSE(3,c,d,:)), clr{c}{d},'linewidth',1.5);
        L{d} = sprintf('\\Delta=%d',deltas(d));
        ylim([0 1]);
        if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
    end
    legend(h,L,'location','northwest');
    xlabel('heading angle (deg)'); ylabel('P(Right)');
    try changeAxesFontSize(gca,15,15); catch; end

    % conf
    if conftask
    subplot(2+double(RTtask),length(cohs),c+length(cohs)); box off; hold on;
    for d = 1:length(deltas)
        beta = [amplConf(3,c,d) muConf(3,c,d) sigmaConf(3,c,d) baselineConf(3,c,d)];
        h(d) = plot(xVals, flippedGauss(beta,xVals), [clr{c}{d}(1) '-'],'linewidth',1.5); hold on;       
        errorbar(hdgs, squeeze(confMean(3,c,d,:)), squeeze(confSE(3,c,d,:)), clr{c}{d},'linewidth',1.5);
        L{d} = sprintf('?=%d',deltas(d));
        ylim([0 1]); hold on;
        if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
    end
%     legend(h,L,'location','northwest');
    xlabel('heading angle (deg)'); 
    if conftask==1, ylabel('SEP (''confidence'', %)');
    elseif conftask==2, ylabel('P(high bet)');
    end
    try changeAxesFontSize(gca,15,15); catch; end
    end

    % RT
    if RTtask
        subplot(2+double(conftask>0),length(cohs),c+length(cohs)*2); box off; hold on;
        for d = 1:length(deltas)
            beta = [amplRT(3,c,d) muRT(3,c,d) sigmaRT(3,c,d) baselineRT(3,c,d)];
            h(d) = plot(xVals, gauss(beta,xVals), [clr{c}{d}(1) '-'], 'linewidth',1.5); hold on;       
            errorbar(hdgs, squeeze(RTmean(3,c,d,:)), squeeze(RTse(3,c,d,:)), clr{c}{d},'linewidth',1.5);
            L{d} = sprintf('?=%d',deltas(d));
            hold on;
            if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
        end
        xlabel('heading angle (deg)'); ylabel('RT (s)')
        if strcmp(subject,'lucio')
            ylim([0.5 0.9])
        else
            ylim([0.6 1.8]);
        end
        try changeAxesFontSize(gca,15,15); catch; end
    end    
    
end

end