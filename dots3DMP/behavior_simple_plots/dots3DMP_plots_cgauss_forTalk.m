exfig=1;


%% first, for all trials irrespective of delta
D = length(deltas)+1; % (the extra column we made for pooling across deltas)
% OR select just delta=0:
% D = find(deltas==0);

         %ves %vis %comb
clr{1} = {[0 0 0], [255 0 255]./255, [0 255 255]./255};
clr{2} = {[0 0 0], [255 0 0]./255, [0 85 255]./255};
clr{3} = {[0 0 0], [155 0 0]./255, [0 0 225]./255};

modtxt = {'ves','vis','comb'};

for c = 1:length(cohs)
    % choice
    clear L h; figure(c); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
    for m = 1:length(mods)     % m c d h
        beta = [muPMF(m,c,D) sigmaPMF(m,c,D)];
        h(m) = plot(xVals, cgauss(beta,xVals), '-', 'Color', clr{c}{m}, 'Linewidth', 3); hold on;
        errorbar(hdgs, squeeze(pRight(m,c,D,:)), squeeze(pRightSE(m,c,D,:)), 'o', 'Color', clr{c}{m}, 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
        L{m} = modtxt{m};
        xlabel('Heading angle (deg)'); ylabel('Proportion rightward choices'); ylim([0 1]);
        changeAxesFontSize(gca,20,20); set(gca,'box','off')
        legend(h,L,'location','northwest'); legend('boxoff');
        if exfig; export_fig(['c=' num2str(cohs(c)) '_' modtxt{m} '_pmf'],'-eps'); end
    end

    % conf
    figure(100+c); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
    for m = 1:length(mods)        
        beta = [amplConf(m,c,D) muConf(m,c,D) sigmaConf(m,c,D) baselineConf(m,c,D)];
        h(d) = plot(xVals, flippedGauss(beta,xVals), '-', 'Color', clr{c}{m}, 'Linewidth', 3); hold on;
        errorbar(hdgs, squeeze(confMean(m,c,D,:)), squeeze(confSE(m,c,D,:)), 'o', 'Color', clr{c}{m}, 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
        xlabel('Heading angle (deg)'); ylabel('Confidence'); ylim([0 1]);
        changeAxesFontSize(gca,20,20); set(gca,'box','off')
        if exfig; export_fig(['c=' num2str(cohs(c)) '_' modtxt{m} '_conf'],'-eps'); end
    end

    % RT
    if ~isnan(RTmean(1,1,2))
        figure(200+c); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
        for m = 1:length(mods)        
            beta = [amplRT(m,c,D) muRT(m,c,D) sigmaRT(m,c,D) baselineRT(m,c,D)];
            h(d) = plot(xVals, gauss(beta,xVals), '-', 'Color', clr{c}{m}, 'Linewidth', 3); hold on;
            errorbar(hdgs, squeeze(RTmean(m,c,D,:)), squeeze(RTse(m,c,D,:)), 'o', 'Color', clr{c}{m}, 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
            xlabel('Heading angle (deg)'); ylabel('Reaction time (s)'); ylim([0.5 0.9]);
            changeAxesFontSize(gca,20,20); set(gca,'box','off')
            if exfig; export_fig(['c=' num2str(cohs(c)) '_' modtxt{m} '_rt'],'-eps'); end
        end
    end
end


% % alternate version: single-cues only
% modtxt = {'ves','vis: low coh','vis: high coh'};
% clr = {'k','m','b'};
% 
% % choice
% clear L h;
% figure(990); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
% for m = 1:3
%     if m==1
%         c = 1;
%         I = data.modality==mods(m);
%     else
%         c = m-1;
%         I = data.modality==mods(m) & data.coherence==cohs(c);
%     end
%     [beta,fval,exitflag,output] = fminsearch(@(x) cgauss_err(x,data.choice(I)==2,data.heading(I)), [0 3]);
%     h(m) = plot(xVals, cgauss(beta,xVals), '-', 'Color', clr{m}, 'Linewidth', 3); hold on;
%     errorbar(hdgs, squeeze(pRight(m,c,D,:)), squeeze(pRightSE(m,c,D,:)), 'o', 'Color', clr{m}, 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
%     L{m} = modtxt{m};
%     xlabel('Heading angle (deg)'); ylabel('Proportion rightward choices'); ylim([0 1]);
%     changeAxesFontSize(gca,20,20); set(gca,'box','off')
%     ll = legend(h,L,'location',[0.1849    0.7938    0.3972    0.2219]); legend('boxoff');
%     if exfig; export_fig(['altPMF' num2str(m)],'-eps'); end
% end
% 
% % conf
% figure(991); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
% for m = 1:3     
%     if m==1
%         c = 1;
%         I = data.modality==mods(m);
%     else
%         c = m-1;
%         I = data.modality==mods(m) & data.coherence==cohs(c);
%     end
%     [beta,fval,exitflag,output] = fminsearch(@(x) flippedGauss_err(x,data.conf(I),data.heading(I)), guess_fgauss);
%     h(d) = plot(xVals, flippedGauss(beta,xVals), '-', 'Color', clr{m}, 'Linewidth', 3); hold on;
%     errorbar(hdgs, squeeze(confMean(m,c,D,:)), squeeze(confSE(m,c,D,:)), 'o', 'Color', clr{m}, 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
%     xlabel('Heading angle (deg)'); ylabel('Confidence'); ylim([0 1]);
%     changeAxesFontSize(gca,20,20); set(gca,'box','off')
%     if exfig; export_fig(['altConf' num2str(m)],'-eps'); end;
% end





% % % % TEMP: pool all conditions for conf vs. RT, globally
% % % for h = 1:length(hdgs)
% % %     J = data.heading==hdgs(h);
% % %     nn(h) = sum(J);
% % %     pR(h) = sum(J & data.choice==2) / nn(h); % 2 is rightward
% % %     RTm(h) = mean(data.RT(J));
% % %     RTs(h) = std(data.RT(J))/sqrt(nn(h));
% % %     confM(h) = mean(data.conf(J));
% % %     confS(h) = std(data.conf(J))/sqrt(nn(h));
% % % end
% % % pRS = sqrt( (pR.*(1-pR)) ./ nn );
% % % 
% % % X = data.heading;
% % % y = data.choice==2; % 2 is rightward
% % % [BB, ~, ~] = glmfit(X, y, 'binomial');
% % % yV = glmval(BB,hdgs(1):0.1:hdgs(end),'logit');
% % %         
% % % figure(992); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
% % % h(m) = plot(xVals,yV,'k-','Linewidth', 3); hold on;
% % % % OR
% % % % [beta,~,~,~] = fminsearch(@(x) cgauss_err(x,data.choice==2,data.heading), [0 3]);
% % % % h(m) = plot(xVals, cgauss(beta,xVals), '-', 'Color', clr{c}{m}, 'Linewidth', 3); hold on;
% % % errorbar(hdgs, pR, pRS, 'o', 'Color', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
% % % xlabel('Heading angle (deg)'); ylabel('Proportion rightward choices'); ylim([0 1]);
% % % changeAxesFontSize(gca,20,20); set(gca,'box','off')
% % % export_fig(['choiceGlobal' num2str(m)],'-eps');
% % % 
% % % figure(993); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
% % % % [beta,~,~,~] = fminsearch(@(x) flippedGauss_err(x,data.conf,data.heading), guess_fgauss);
% % % % h(d) = plot(xVals, flippedGauss(beta,xVals), '-', 'Color', 'k', 'Linewidth', 3); hold on;
% % % errorbar(hdgs, confM, confS, 'o-', 'Color', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
% % % xlabel('Heading angle (deg)'); ylabel('Confidence'); ylim([0 1]);
% % % changeAxesFontSize(gca,20,20); set(gca,'box','off')
% % % export_fig(['confGlobal' num2str(m)],'-eps');
% % % 
% % % figure(994); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
% % % % [beta,fval,exitflag,output] = fminsearch(@(x) gauss_err(x,data.RT,data.heading), guess_gauss);
% % % % h(d) = plot(xVals, gauss(beta,xVals), '-', 'Color', 'k', 'Linewidth', 3); hold on;
% % % errorbar(hdgs, RTm, RTs, 'o-', 'Color', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
% % % xlabel('Heading angle (deg)'); ylabel('Reaction time (s)');
% % % changeAxesFontSize(gca,20,20); set(gca,'box','off')
% % % export_fig(['RTglobal' num2str(m)],'-eps');






%% now separate by delta

if length(deltas)>1

% D = 1:length(deltas);
% OR omit the zero delta, for some plots
D = find(deltas~=0)';

    
% clr{1} = {[255 140 0]./255, [175 238 238]./255, [221 160 221]./255}; % orange, cyan, lavender
clr{1} = {[1 0 0], [0 1 0], [0 0 1]}; % red, green, blue
clr{2} = clr{1};
clr{3} = clr{1};

for c = 1:length(cohs)
    % choice
    clear L h; figure(300+c); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
%     for d = 1:length(deltas)     % m c d h
    lind = 1;
    for d = D     % m c d h
        beta = [muPMF(3,c,d) sigmaPMF(3,c,d)];
        h(lind) = plot(xVals, cgauss(beta,xVals), '-', 'Color', clr{c}{d}, 'Linewidth', 3); hold on;
        errorbar(hdgs, squeeze(pRight(3,c,d,:)), squeeze(pRightSE(3,c,d,:)), 'o', 'Color', clr{c}{d}, 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
        L{lind} = sprintf('\\Delta=%d',deltas(d)); lind = lind+1;
        xlabel('Heading angle (deg)'); ylabel('Proportion rightward choices'); ylim([0 1]);
        changeAxesFontSize(gca,20,20); set(gca,'box','off')
        legend(h,L,'location','northwest'); legend('boxoff');
        if exfig; export_fig(['c=' num2str(cohs(c)) '_d=' num2str(deltas(d)) '_pmf'],'-eps'); end
    end

    % conf
    figure(400+c); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
    for d = D
        beta = [amplConf(3,c,d) muConf(3,c,d) sigmaConf(3,c,d) baselineConf(3,c,d)];
        plot(xVals, flippedGauss(beta,xVals), '-', 'Color', clr{c}{d}, 'Linewidth', 3); hold on;
        errorbar(hdgs, squeeze(confMean(3,c,d,:)), squeeze(confSE(3,c,d,:)), 'o', 'Color', clr{c}{d}, 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
        xlabel('Heading angle (deg)'); ylabel('Confidence'); ylim([0 1]);
        changeAxesFontSize(gca,20,20); set(gca,'box','off')
        if exfig; export_fig(['c=' num2str(cohs(c)) '_d=' num2str(deltas(d)) '_conf'],'-eps'); end
    end
    
    % RT
    if ~isnan(RTmean(1,1,2))
        figure(500+c); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
        for d = D
            beta = [amplRT(3,c,d) muRT(3,c,d) sigmaRT(3,c,d) baselineRT(3,c,d)];
            plot(xVals, gauss(beta,xVals), '-', 'Color', clr{c}{d}, 'Linewidth', 3); hold on;
            errorbar(hdgs, squeeze(RTmean(3,c,d,:)), squeeze(RTse(3,c,d,:)), 'o', 'Color', clr{c}{d}, 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
            xlabel('Heading angle (deg)'); ylabel('Reaction time (s)');
            changeAxesFontSize(gca,20,20); set(gca,'box','off')
            if exfig; export_fig(['c=' num2str(cohs(c)) '_d=' num2str(deltas(d)) '_rt'],'-eps'); end
        end
    end    
end

end