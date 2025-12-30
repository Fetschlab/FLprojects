function plotBehav(parsedData,cohs,conftask,RTtask,wFit)


xl = [-1.05*max(cohs) 1.05*max(cohs)];

if nargin<5, wFit=0; end

nplots = 1+RTtask+double(conftask>0);

if conftask==1
    conf = parsedData.confMean;
    confSE = parsedData.confSE;
    confCorr = parsedData.confMeanCorr;
    confSEcorr = parsedData.confSEcorr;
    confErr = parsedData.confMeanErr;
    confSEerr = parsedData.confSEerr;
    confAxisLabel = 'Confidence (mean SEP or rating)';
    legEntry{1} = 'high conf'; legEntry{2} = 'low conf';
end
if conftask==2
    conf = parsedData.pHigh;
    confSE = parsedData.pHighSE;
    confCorr = parsedData.pHighCorr;
    confSEcorr = parsedData.pHighSEcorr;
    confErr = parsedData.pHighErr;
    confSEerr = parsedData.pHighSEerr;
    confAxisLabel = 'Proportion high bet';
    legEntry{1} = 'high bet'; legEntry{2} = 'low bet';
%    confAxisLabel = 'proportion high conf rating'; % temp, for doubtConf
%    legEntry{1} = 'high conf'; legEntry{2} = 'low conf'; % temp, for doubtConf
end
    
figure(101);
set(gcf,'Color',[1 1 1],'Position',[300 500 450 200+200*nplots],'PaperPositionMode','auto'); clf;

% which cohs to label on the x axis:
if length(cohs)==12
    ticks = cohs([1 2 4 9 11 12]); 
elseif length(cohs)==11
    ticks = cohs([1 2 4 6 8 10 11]); 
elseif length(cohs)==10 % doubtconf
    ticks = [cohs(1:5)' 0 cohs(6:10)'];
%     xlabels = {'-0.7','-0.45','-0.25','','','0','','','0.25','0.45','0.7'};
    xlabels = {'-0.7','-0.45','','-0.15','','0','','0.15','','0.45','0.7'};
else
    ticks = cohs;
end

if wFit
    line1=[];
    line2=[];
else
    line1='-';
    line2='--';
end

subplot(nplots,1,1);
if wFit==0; plot(parsedData.xVals,parsedData.yVals1,'k-'); hold on; end
errorbar(cohs, parsedData.pRight, parsedData.pRightSE, 'ko'); hold on;
try
    set(gca,'xtick',ticks,'tickdir','out','xticklabel',xlabels);
catch
    set(gca,'xtick',ticks,'tickdir','out');
end
ylim([0 1]); xlim(xl);
xlabel('Motion strength (%coh)');
ylabel('Proportion rightward choices');
changeAxesFontSize(gca,14,14);

if RTtask
    subplot(nplots,1,2);
    errorbar(cohs, parsedData.RTmean, parsedData.RTse, ['bo' line1]);
    try
        set(gca,'xtick',ticks,'tickdir','out','xticklabel',xlabels);
    catch
        set(gca,'xtick',ticks,'tickdir','out');
    end
    xlim(xl);
    xlabel('Motion strength (%coh)'); ylabel('Reaction time (s)');
    changeAxesFontSize(gca,14,14);
    if conftask
        subplot(nplots,1,3);
        errorbar(cohs, conf, confSE, ['ro' line1]);
        try
            set(gca,'xtick',ticks,'tickdir','out','xticklabel',xlabels);
        catch
            set(gca,'xtick',ticks,'tickdir','out');
        end
        ylim([0 1]); xlim(xl);
        xlabel('Motion strength (%coh)'); ylabel(confAxisLabel);
        changeAxesFontSize(gca,14,14);
    end    
else 
    if conftask
        subplot(nplots,1,2);
        errorbar(cohs, conf, confSE, ['ro' line1]);
        try
            set(gca,'xtick',ticks,'tickdir','out','xticklabel',xlabels);
        catch
            set(gca,'xtick',ticks,'tickdir','out');
        end
        ylim([0 1]); xlim(xl);
        xlabel('Motion strength (%coh)'); ylabel(confAxisLabel);
        changeAxesFontSize(gca,14,14);
    end
end


%% separate choice+RT by high/low bet, and p(high) (or conf rating) by corr/err

if conftask

figure(102);
set(gcf,'Color',[1 1 1],'Position',[600 500 450 200+200*nplots],'PaperPositionMode','auto'); clf;

subplot(nplots,1,1);
if wFit==0; plot(parsedData.xVals,parsedData.yVals2,'k-',parsedData.xVals,parsedData.yVals3,'k--'); hold on; end
errorbar(cohs, parsedData.pRightHigh, parsedData.pRightSEhigh, 'ko', 'MarkerFaceColor', 'k'); hold on;
errorbar(cohs, parsedData.pRightLow, parsedData.pRightSElow, 'ko', 'MarkerFaceColor', 'w');
try
    set(gca,'xtick',ticks,'tickdir','out','xticklabel',xlabels);
catch
    set(gca,'xtick',ticks,'tickdir','out');
end
ylim([0 1]); xlim(xl);
xlabel('Motion strength (%coh)'); ylabel('Proportion rightward choices');
legend(legEntry{1},legEntry{2},'Location','Northwest'); legend('boxoff')
changeAxesFontSize(gca,14,14);

if RTtask
    subplot(nplots,1,2);
    errorbar(cohs, parsedData.RTmeanHigh, parsedData.RTseHigh, ['bo' line1], 'MarkerFaceColor', 'b'); hold on;
    errorbar(cohs, parsedData.RTmeanLow, parsedData.RTseLow, ['bo' line2], 'MarkerFaceColor', 'w');
    try
        set(gca,'xtick',ticks,'tickdir','out','xticklabel',xlabels);
    catch
        set(gca,'xtick',ticks,'tickdir','out');
    end
    xlim(xl);
    xlabel('Motion strength (%coh)'); ylabel('Reaction time (s)');
    legend(legEntry{1},legEntry{2},'Location','Northwest'); legend('boxoff')
    changeAxesFontSize(gca,14,14);
    if conftask
        subplot(nplots,1,3);
        errorbar(cohs, confCorr, confSEcorr, ['ro' line1], 'MarkerFaceColor', 'r'); hold on;
        errorbar(cohs, confErr, confSEerr, ['ro' line2], 'MarkerFaceColor', 'w');
        try
            set(gca,'xtick',ticks,'tickdir','out','xticklabel',xlabels);
        catch
            set(gca,'xtick',ticks,'tickdir','out');
        end
        ylim([0 1]); xlim(xl);
        xlabel('Motion strength (%coh)'); ylabel(confAxisLabel);
        legend('Corr','Err','Location','Southwest'); legend('boxoff')
        changeAxesFontSize(gca,14,14);
    end
else
    if conftask
        subplot(nplots,1,2);
        errorbar(cohs, confCorr, confSEcorr, ['ro' line1], 'MarkerFaceColor', 'r'); hold on;
        errorbar(cohs, confErr, confSEerr, ['ro' line2], 'MarkerFaceColor', 'w');
        try
            set(gca,'xtick',ticks,'tickdir','out','xticklabel',xlabels);
        catch
            set(gca,'xtick',ticks,'tickdir','out');
        end
        ylim([0 1]); xlim(xl);
        xlabel('Motion strength (%coh)'); ylabel(confAxisLabel);
        legend('Corr','Err','Location','North'); legend('boxoff')
        changeAxesFontSize(gca,14,14);
    end
end

end


