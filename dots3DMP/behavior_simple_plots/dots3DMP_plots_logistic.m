function dots3DMP_plots_logistic(betaLow,betaHigh,parsedData,mods,cohs,deltas,hdgs)

%% first, for all trials irrespective of delta
D = length(deltas)+1; % (the extra column we made for pooling across deltas)
% OR select just delta=0:
% D = find(deltas==0);

modlabels = {'Ves','Vis','Comb'};

yVals1 = glmval(B1,xVals,'logit');
% [betaHigh, devHigh, statsHigh] = glmfit(X, y, 'binomial','link',linkfun);


         %ves %vis %comb
clr{1} = {'ko','mo','co'};
clr{2} = {'ko','ro','bo'};
clr{3} = {'ko','yo','go'};
figure(101+D);
set(gcf,'Color',[1 1 1],'Position',[300 1000 450+300*(length(cohs)-1) 400],'PaperPositionMode','auto'); clf;
for c = 1:length(cohs)
    % choice
    subplot(2,length(cohs),c); box off; hold on;
    for m = 1:length(mods)     % m c d h
        beta = [gfit.muPMF(m,c,D) gfit.sigmaPMF(m,c,D)];
        h(m) = plot(parsedData.xVals, gfit.func.cgauss(beta,parsedData.xVals), [clr{c}{m}(1) '-'],'linewidth',1.5); hold on;
        errorbar(hdgs, squeeze(parsedData.pRight(m,c,D,:)), squeeze(parsedData.pRightSE(m,c,D,:)), clr{c}{m},'linewidth',1.5);
        ylim([0 1]); if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
        text(hdgs(1)+0.5,1.0-m*0.07,sprintf('%s: mu = %.2f, s = %.2f',modlabels{m},beta(1),beta(2)),'color',clr{c}{m}(1))
    end
%     legend(h,'vestib','visual','comb','Location','northwest');
    
    if length(deltas)==1,xlabel('heading angle (deg)'); end
    ylabel('proportion rightward choices');
    try, changeAxesFontSize(gca,15,15); end
end


%% now separate by delta

if length(deltas)>1
    
clr{1} = {'bs','cs','gs'};
clr{2} = {'b^','c^','g^'};
clr{3} = {'bo','co','go'};

clear L;
%figure(208);
%set(gcf,'Color',[1 1 1],'Position',[50 20 950+300*(length(cohs)-2) 800],'PaperPositionMode','auto'); clf;
for c = 1:length(cohs)
    % choice
    subplot(2,length(cohs),c+length(cohs)); box off; hold on;
    for d = 1:length(deltas)     % m c d h
        beta = [gfit.muPMF(3,c,d) gfit.sigmaPMF(3,c,d)];
        h(d) = plot(parsedData.xVals, gfit.func.cgauss(beta,parsedData.xVals), [clr{c}{d}(1) '-'],'linewidth',1.5); hold on;
        errorbar(hdgs, squeeze(parsedData.pRight(3,c,d,:)), squeeze(parsedData.pRightSE(3,c,d,:)), clr{c}{d},'linewidth',1.5);
        L{d} = sprintf('\\Delta=%d',deltas(d));
        ylim([0 1]);
        %if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
    end
    legend(h,L,'location','northwest');
    xlabel('heading angle (deg)'); %ylabel('proportion rightward choices');
    try, changeAxesFontSize(gca,15,15); end    
end

end