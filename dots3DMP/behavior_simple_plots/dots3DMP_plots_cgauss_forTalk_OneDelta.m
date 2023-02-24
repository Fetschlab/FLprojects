exfig=1;

%% now separate by delta

% D = 1:length(deltas);
% OR omit the zero delta, for some plots
D = find(deltas~=0)';

    
% clr{1} = {[255 140 0]./255, [175 238 238]./255, [221 160 221]./255}; % orange, cyan, lavender
clr{1} = {[1 0 0], [0 0 1]}; % red, blue
clr{2} = clr{1};
clr{3} = clr{1};

for d = D
    % choice
    clear L h; figure(300+d); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
    lind = 1;
    for c = 1:length(cohs)
        beta = [muPMF(3,c,d) sigmaPMF(3,c,d)];
        h(lind) = plot(xVals, cgauss(beta,xVals), '-', 'Color', clr{d}{c}, 'Linewidth', 3); hold on;
        errorbar(hdgs, squeeze(pRight(3,c,d,:)), squeeze(pRightSE(3,c,d,:)), 'o', 'Color', clr{d}{c}, 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
        L{lind} = sprintf('Coh=%0.1f',cohs(c)); lind = lind+1;
        xlabel('Heading angle (deg)'); ylabel('Proportion rightward choices'); ylim([0 1]);
        changeAxesFontSize(gca,20,20); set(gca,'box','off')
        legend(h,L,'location','northwest'); legend('boxoff');
    end
    if exfig; export_fig(['d=' num2str(deltas(d)) '_pmf'],'-eps'); end

    % conf
    figure(400+d); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
    for c = 1:length(cohs)
        beta = [amplConf(3,c,d) muConf(3,c,d) sigmaConf(3,c,d) baselineConf(3,c,d)];
        plot(xVals, flippedGauss(beta,xVals), '-', 'Color', clr{d}{c}, 'Linewidth', 3); hold on;
        errorbar(hdgs, squeeze(confMean(3,c,d,:)), squeeze(confSE(3,c,d,:)), 'o', 'Color', clr{d}{c}, 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
        xlabel('Heading angle (deg)'); ylabel('Confidence'); ylim([0 1]);
        changeAxesFontSize(gca,20,20); set(gca,'box','off')
    end
    if exfig; export_fig(['d=' num2str(deltas(d)) '_conf'],'-eps'); end
    
    % RT
    if ~isnan(RTmean(1,1,2))
        figure(500+d); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
        for c = 1:length(cohs)
            beta = [amplRT(3,c,d) muRT(3,c,d) sigmaRT(3,c,d) baselineRT(3,c,d)];
            plot(xVals, gauss(beta,xVals), '-', 'Color', clr{d}{c}, 'Linewidth', 3); hold on;
            errorbar(hdgs, squeeze(RTmean(3,c,d,:)), squeeze(RTse(3,c,d,:)), 'o', 'Color', clr{d}{c}, 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
            xlabel('Heading angle (deg)'); ylabel('Reaction time (s)');
            changeAxesFontSize(gca,20,20); set(gca,'box','off')
        end
        if exfig; export_fig(['d=' num2str(deltas(d)) '_rt'],'-eps'); end
    end    
end
