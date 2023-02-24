function dots3DMP_plots_cgauss_byModDelta(gfit,parsedData,mods,cohs,deltas,hdgs,conftask,RTtask)
% SJ 07-2021 converted to function for cleaner workspace

% first, for all trials irrespective of delta
D = length(deltas)+1; % (the extra column we made for pooling across deltas)
% OR select just delta=0:
% D = find(deltas==0);

modlabels = {'Ves','Vis','Comb'};
cohlabels = {'Low','High'};

clinstl = {':','-'};

         %ves %vis %comb
clr{1} = {'ko','mo','co'};
clr{2} = {'ko','ro','bo'};
clr{3} = {'ko','yo','go'};
% set(gcf,'Color',[1 1 1],'Position',[300 1000 450+300*(length(cohs)-2) 400+150*(conftask>0)+150*RTtask],'PaperPositionMode','auto'); clf;
for m = 1:length(mods)
    
    % choice
    fh=figure(101);
    set(fh,'color','w','position',[100 100 700 500])
    subplot(2,length(mods),m); box off; hold on;
    for c = 1:length(cohs) % m c d h
        if mods(m)==1 && c>1, continue, end 
        beta = [gfit.choice.mu(m,c,D) gfit.choice.sigma(m,c,D)];
        h(m) = plot(parsedData.xVals, gfit.choice.func(beta,parsedData.xVals), [clr{c}{m}(1) clinstl{c}],'linewidth',1.5); hold on;
        errorbar(hdgs, squeeze(parsedData.pRight(m,c,D,:)), squeeze(parsedData.pRightSE(m,c,D,:)), clr{c}{m},'linewidth',1.5);
        ylim([0 1]); 
%         text(hdgs(1)+0.5,1.0-c*0.07,sprintf('%s: mu = %.2f, s = %.2f',cohlabels{c},beta(1),beta(2)),'color',clr{c}{m}(1))
    end
    title(modlabels{m});    
    xlabel('heading angle (deg)');
    ylabel('P(right)');
    try changeAxesFontSize(gca,15,15); catch; end
    
        % conf
    if conftask
        fh=figure(102);
        set(fh,'color','w','position',[100 100 700 500])
        subplot(2,length(mods),m); box off; hold on;
        for c = 1:length(cohs)
            if mods(m)==1 && c>1, continue, end
            beta = [gfit.conf.ampl(m,c,D) gfit.conf.mu(m,c,D) gfit.conf.sigma(m,c,D) gfit.conf.bsln(m,c,D)];
            h(m) = plot(parsedData.xVals, gfit.conf.func(beta,parsedData.xVals), [clr{c}{m}(1) clinstl{c}],'linewidth',1.5); hold on;
            errorbar(hdgs, squeeze(parsedData.confMean(m,c,D,:)), squeeze(parsedData.confSE(m,c,D,:)), clr{c}{m},'linewidth',1.5);
            ylim([0 max(max(parsedData.confMean(:)),1)]); title(modlabels{m});
%             text(hdgs(1)+0.5,0.5-c*0.07,sprintf('%s: mu = %.2f, s = %.2f',cohlabels{c},beta(1),beta(2)),'color',clr{c}{m}(1))
        end
        title(modlabels{m});
        %     legend(h,'vestib','visual','comb','location','northwest');
        xlabel('heading angle (deg)');
        if conftask==1, ylabel('SEP (''confidence'', %)');
        elseif conftask==2, ylabel('P(high bet)');
        end
        try changeAxesFontSize(gca,15,15); catch; end
    end
    
    % RT
    if RTtask
        fh=figure(103);
        set(fh,'color','w','position',[100 100 700 500])
        subplot(2,length(mods),m); box off; hold on;
        for c = 1:length(cohs) 
            if mods(m)==1 && c>1, continue, end
            beta = [gfit.RT.ampl(m,c,D) gfit.RT.mu(m,c,D) gfit.RT.sigma(m,c,D) gfit.RT.bsln(m,c,D)];
            h(m) = plot(parsedData.xVals, gfit.RT.func(beta,parsedData.xVals), [clr{c}{m}(1) clinstl{c}],'linewidth',1.5); hold on;       
            errorbar(hdgs, squeeze(parsedData.RTmean(m,c,D,:)), squeeze(parsedData.RTse(m,c,D,:)), clr{c}{m},'linewidth',1.5);
%             text(hdgs(1)+0.5,0.4+c*0.07,sprintf('%s: mu = %.2f, s = %.2f',cohlabels{c},beta(1),beta(2)),'color',clr{c}{m}(1))
        end
        title(modlabels{m});
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

clear L;
for d=1:length(deltas)
   
    % choice
    figure(101);
    subplot(2,length(mods),length(mods)+d); box off; hold on; % assumes 3 deltas and 3 modalities!
    for c = 1:length(cohs)
        beta = [gfit.choice.mu(3,c,d) gfit.choice.sigma(3,c,d)];
        h(d) = plot(parsedData.xVals, gfit.choice.func(beta,parsedData.xVals), [clr{c}{d}(1) clinstl{c}],'linewidth',1.5); hold on;
        errorbar(hdgs, squeeze(parsedData.pRight(3,c,d,:)), squeeze(parsedData.pRightSE(3,c,d,:)), clr{c}{d},'linewidth',1.5);
        ylim([0 1]);
%         text(hdgs(1)+0.5,1.0-c*0.07,sprintf('%s: mu = %.2f, s = %.2f',cohlabels{c},beta(1),beta(2)),'color',clr{c}{d}(1))
    end
    title(sprintf('\\Delta=%d',deltas(d)))
    xlabel('heading angle (deg)'); ylabel('P(Right)');
    try changeAxesFontSize(gca,15,15); catch; end

    % conf
    if conftask
    
        figure(102)
        subplot(2,length(mods),length(mods)+d); box off; hold on; % assumes 3 deltas and 3 modalities!
        for c = 1:length(cohs)
            beta = [gfit.conf.ampl(3,c,d) gfit.conf.mu(3,c,d) gfit.conf.sigma(3,c,d) gfit.conf.bsln(3,c,d)];
            h(d) = plot(parsedData.xVals, gfit.conf.func(beta,parsedData.xVals), [clr{c}{d}(1) clinstl{c}],'linewidth',1.5); hold on;
            errorbar(hdgs, squeeze(parsedData.confMean(3,c,d,:)), squeeze(parsedData.confSE(3,c,d,:)), clr{c}{d},'linewidth',1.5);
            ylim([0 1]); hold on;
%             text(hdgs(1)+0.5,1.0-c*0.07,sprintf('%s: mu = %.2f, s = %.2f',cohlabels{c},beta(1),beta(2)),'color',clr{c}{d}(1))
        end
       
    title(sprintf('\\Delta=%d',deltas(d)))
    xlabel('heading angle (deg)'); 
    if conftask==1, ylabel('SEP (''confidence'', %)');
    elseif conftask==2, ylabel('P(high bet)');
    end
    try changeAxesFontSize(gca,15,15); catch; end
    end

    % RT
    if RTtask
        figure(103)
        subplot(2,length(mods),length(mods)+d); box off; hold on; % assumes 3 deltas and 3 modalities!
        for c = 1:length(cohs)
            beta = [gfit.RT.ampl(3,c,d) gfit.RT.mu(3,c,d) gfit.RT.sigma(3,c,d) gfit.RT.bsln(3,c,d)];
            h(d) = plot(parsedData.xVals, gfit.RT.func(beta,parsedData.xVals), [clr{c}{d}(1) clinstl{c}],'linewidth',1.5); hold on;
            errorbar(hdgs, squeeze(parsedData.RTmean(3,c,d,:)), squeeze(parsedData.RTse(3,c,d,:)), clr{c}{d},'linewidth',1.5);
%             text(hdgs(1)+0.5,1.0-c*0.07,sprintf('%s: mu = %.2f, s = %.2f',cohlabels{c},beta(1),beta(2)),'color',clr{c}{d}(1))
        end
        title(sprintf('\\Delta=%d',deltas(d)))
        xlabel('heading angle (deg)'); ylabel('RT (s)')
        yRng = [min(parsedData.RTmean(:)) max(parsedData.RTmean(:))];
        ylim(yRng.*[0.9 1.1])
        try changeAxesFontSize(gca,15,15); catch; end
    end    

end

end