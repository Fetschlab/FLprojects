function dots3DMP_plots_cgauss_byCoh(gfit,parsedData,mods,cohs,deltas,hdgs,conftask,RTtask)
% SJ 07-2021 converted to function for cleaner workspace

% first, for all trials irrespective of delta
% D = length(deltas)+1; % (the extra column we made for pooling across deltas)
% OR select just delta=0:
% D = find(deltas==0);
D = gfit.D;

% modlabels = {'Ves','Vis','Comb'};

         %ves %vis %comb
clr{1} = {'ko','mo','co'};
clr{2} = {'ko','ro','bo'};
clr{3} = {'ko','yo','go'};
clr{4} = clr{1};
figure(101+D);
<<<<<<< HEAD:dots3DMP/dots3DMP_plots_cgauss_byCoh.m
set(gcf,'Color',[1 1 1],'Position',[300 1000 450+300*(length(cohs)-2) 400+150*(conftask>0)+150*RTtask],'PaperPositionMode','auto'); clf;
=======
set(gcf,'Color',[1 1 1],'Position',[300 1000 450+300*(length(cohs)-2) 200+150*(conftask>0)+150*RTtask],'PaperPositionMode','auto'); clf;
% set(gcf,'Color',[1 1 1],'Position',[300 600 600 800],'PaperPositionMode','auto'); clf;
>>>>>>> 442f2b6e8eab5ad9e235836908aa3f5d8c3a02c4:dots3DMP/dots3DMP_plots_cgauss_func.m
for c = 1:length(cohs)
    % choice
    subplot(1+double(conftask>0)+double(RTtask),length(cohs),c); box off; hold on;
    for m = 1:length(mods)     % m c d h
        beta = [gfit.choice.mu(m,c,D) gfit.choice.sigma(m,c,D)];
        h(m) = plot(parsedData.xVals, gfit.choice.func(beta,parsedData.xVals), [clr{c}{m}(1) '-'],'linewidth',1.5); hold on;
        errorbar(hdgs, squeeze(parsedData.pRight(m,c,D,:)), squeeze(parsedData.pRightSE(m,c,D,:)), clr{c}{m},'linewidth',1.5);
        ylim([0 1]); if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
%         text(hdgs(1)+0.5,1.0-m*0.07,sprintf('%s: mu = %.2f, s = %.2f',modlabels{m},beta(1),beta(2)),'color',clr{c}{m}(1))
    end
%     legend(h,'vestib','visual','comb','Location','northwest');
    
    xlabel('heading angle (deg)');
    ylabel('P(right)');
    try changeAxesFontSize(gca,15,15); catch; end
    
        % conf
    if conftask
        subplot(2+double(RTtask),length(cohs),c+length(cohs)); box off; hold on;
        for m = 1:length(mods)
            beta = [gfit.conf.ampl(m,c,D) gfit.conf.mu(m,c,D) gfit.conf.sigma(m,c,D) gfit.conf.bsln(m,c,D)];
            h(m) = plot(parsedData.xVals, gfit.conf.func(beta,parsedData.xVals), [clr{c}{m}(1) '-'],'linewidth',1.5); hold on;
            errorbar(hdgs, squeeze(parsedData.confMean(m,c,D,:)), squeeze(parsedData.confSE(m,c,D,:)), clr{c}{m},'linewidth',1.5);
            ylim([0 max(max(parsedData.confMean(:)),1)]); if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
        end
        %     legend(h,'vestib','visual','comb','location','northwest');
        xlabel('heading angle (deg)');
        if conftask==1, ylabel('SEP (''confidence'', %)');
        elseif conftask==2, ylabel('P(high bet)');
        end
        try changeAxesFontSize(gca,15,15); catch; end
    end
    
    % RT
    if RTtask
        subplot(2+double(conftask>0),length(cohs),c+length(cohs)*2); box off; hold on;
        for m = 1:length(mods)        
            beta = [gfit.RT.ampl(m,c,D) gfit.RT.mu(m,c,D) gfit.RT.sigma(m,c,D) gfit.RT.bsln(m,c,D)];
            h(m) = plot(parsedData.xVals, gfit.RT.func(beta,parsedData.xVals), [clr{c}{m}(1) '-'],'linewidth',1.5); hold on;       
            errorbar(hdgs, squeeze(parsedData.RTmean(m,c,D,:)), squeeze(parsedData.RTse(m,c,D,:)), clr{c}{m},'linewidth',1.5);
            if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
        end
        yRng = [min(parsedData.RTmean(:)) max(parsedData.RTmean(:))];
        ylim(yRng.*[0.9 1.1])
        xlabel('heading angle (deg)'); ylabel('RT (s)')
        try changeAxesFontSize(gca,15,15); catch; end
    end
end


%% now separate by delta

if length(deltas)>1
    
clr{1} = {'bs','cs','gs'};
clr{2} = {'b^','c^','g^'};
clr{3} = {'bo','co','go'};
clr{4} = {'bd','cd','gd'};

clear L;
figure(208);
set(gcf,'Color',[1 1 1],'Position',[50 20 450+300*(length(cohs)-2) 200+150*(conftask>0)+150*RTtask],'PaperPositionMode','auto'); clf;
% set(gcf,'Color',[1 1 1],'Position',[900 600 600 800],'PaperPositionMode','auto'); clf;
for c = 1:length(cohs)
    % choice
    subplot(1+double(conftask>0)+double(RTtask),length(cohs),c); box off; hold on;
    for d = 1:length(deltas)     % m c d h
        beta = [gfit.choice.mu(3,c,d) gfit.choice.sigma(3,c,d)];
        h(d) = plot(parsedData.xVals, gfit.choice.func(beta,parsedData.xVals), [clr{c}{d}(1) '-'],'linewidth',1.5); hold on;
        errorbar(hdgs, squeeze(parsedData.pRight(3,c,d,:)), squeeze(parsedData.pRightSE(3,c,d,:)), clr{c}{d},'linewidth',1.5);
        L{d} = sprintf('\\Delta=%d',deltas(d));
        ylim([0 1]);
        if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
    end
    lh=legend(h,L,'location','southeast'); set(lh,'box','off');
    xlabel('heading angle (deg)'); ylabel('P(Right)');
    try changeAxesFontSize(gca,15,15); catch; end

    % conf
    if conftask
        subplot(2+double(RTtask),length(cohs),c+length(cohs)); box off; hold on;
        for d = 1:length(deltas)
            beta = [gfit.conf.ampl(3,c,d) gfit.conf.mu(3,c,d) gfit.conf.sigma(3,c,d) gfit.conf.bsln(3,c,d)];
            h(d) = plot(parsedData.xVals, gfit.conf.func(beta,parsedData.xVals), [clr{c}{d}(1) '-'],'linewidth',1.5); hold on;       
            errorbar(hdgs, squeeze(parsedData.confMean(3,c,d,:)), squeeze(parsedData.confSE(3,c,d,:)), clr{c}{d},'linewidth',1.5);
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
            beta = [gfit.RT.ampl(3,c,d) gfit.RT.mu(3,c,d) gfit.RT.sigma(3,c,d) gfit.RT.bsln(3,c,d)];
            h(d) = plot(parsedData.xVals, gfit.RT.func(beta,parsedData.xVals), [clr{c}{d}(1) '-'],'linewidth',1.5); hold on;
            errorbar(hdgs, squeeze(parsedData.RTmean(3,c,d,:)), squeeze(parsedData.RTse(3,c,d,:)), clr{c}{d},'linewidth',1.5);
            L{d} = sprintf('?=%d',deltas(d));
            hold on;
            if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
        end
        xlabel('heading angle (deg)'); ylabel('RT (s)')
        yRng = [min(parsedData.RTmean(:)) max(parsedData.RTmean(:))];
        ylim(yRng.*[0.9 1.1])
        try changeAxesFontSize(gca,15,15); catch; end
    end    

end

end