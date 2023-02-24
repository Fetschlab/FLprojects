% plot curves split by high/low confidence
% bit of a kluge, just splitting by median confidence rating, but eventual
% plotting will depend on exact analyses, so just use this as a sanity
% check that confidence ratings are genuinely saying something about
% perceived p(correct)

% first, for all trials irrespective of delta
D = length(deltas)+1; % (the extra column we made for pooling across deltas)
% OR select just delta=0:
% D = find(deltas==0);

modlabels = {'Ves','Vis','Comb'};

%ves %vis %comb
clr{1} = {'ko','mo','co'};
clr{2} = {'ko','ro','bo'};
clr{3} = {'ko','yo','go'};

figure(111+D);
set(gcf,'Color',[1 1 1],'Position',[300 500 950+300*(length(cohs)-2) 800],'PaperPositionMode','auto'); clf;

for c=1:length(cohs)
    for m=1:length(mods)
        subplot(length(mods),length(cohs),c+(m-1)*length(cohs))
        if m==1 && c~=1, delete(gca); continue, end
        if plotLogistic(m,c,D)
            h1(m) = plot(xVals,squeeze(yVals(m,c,D,:,1)),[clr{c}{m}(1) '-'],'linewidth',1.5); hold on;
            h2(m) = plot(xVals,squeeze(yVals(m,c,D,:,2)),[clr{c}{m}(1) '--'],'linewidth',1.5); hold on;
            
%             errorbar(hdgs, squeeze(pRight(m,c,D,:,1)), squeeze(pRightSE(m,c,D,:,1)), clr{c}{m},'linewidth',1.5,'linestyle','none');
%             errorbar(hdgs, squeeze(pRight(m,c,D,:,2)), squeeze(pRightSE(m,c,D,:,2)), clr{c}{m}(1),'linewidth',1.5,'linestyle','none','marker','d');
            
            plot(hdgs, squeeze(pRight(m,c,D,:,1)), clr{c}{m},'linewidth',1.5,'linestyle','none');
            plot(hdgs, squeeze(pRight(m,c,D,:,2)), clr{c}{m}(1),'linewidth',1.5,'linestyle','none','marker','x');
        
        else
            %h(m) = errorbar(hdgs, squeeze(pRight(m,c,D,:)), squeeze(pRightSE(m,c,D,:)), [clr{c}{m} '-']); hold on;
        end
        ylim([0 1]);
        if m~=1
            title([modlabels{m} ', coh = ' num2str(cohs(c))]); 
        else
            title(modlabels{m})
        end
    legend('High Bet','Low Bet','Location','northwest');
    if m==length(mods), xlabel('heading angle (deg)'); end
    if m==2, ylabel('proportion rightward choices'); end
    changeAxesFontSize(gca,15,15);
    end
end

% first, for all trials irrespective of delta
D = length(deltas)+1; % (the extra column we made for pooling across deltas)
% OR select just delta=0:
% D = find(deltas==0);

modlabels = {'Ves','Vis','Comb'};

%ves %vis %comb
clr{1} = {'ko','mo','co'};
clr{2} = {'ko','ro','bo'};
clr{3} = {'ko','yo','go'};


%%
figure(112+D);
set(gcf,'Color',[1 1 1],'Position',[300 500 700+300*(length(cohs)-2) 500],'PaperPositionMode','auto'); clf;

for c=1:length(cohs)
    for m=1:length(mods)
        subplot(length(mods),length(cohs),c+(m-1)*length(cohs))
        if m==1 && c~=1, delete(gca); continue, end
        if plotLogistic(m,c,D)
            hold on
%             h1(m) = plot(xVals,squeeze(RTmean(m,c,D,:,1)),[clr{c}{m}(1) '-'],'linewidth',1.5); hold on;
%             h2(m) = plot(xVals,squeeze(RTmean(m,c,D,:,2)),[clr{c}{m}(1) '--'],'linewidth',1.5); hold on;
            
%             errorbar(hdgs, squeeze(pRight(m,c,D,:,1)), squeeze(pRightSE(m,c,D,:,1)), clr{c}{m},'linewidth',1.5,'linestyle','none');
%             errorbar(hdgs, squeeze(pRight(m,c,D,:,2)), squeeze(pRightSE(m,c,D,:,2)), clr{c}{m}(1),'linewidth',1.5,'linestyle','none','marker','d');
            
            plot(hdgs, squeeze(RTmean(m,c,D,:,1)), clr{c}{m}(1),'linestyle','-','linewidth',1.5,'marker',clr{c}{m}(2));
            plot(hdgs, squeeze(RTmean(m,c,D,:,2)), clr{c}{m}(1),'linestyle','--','linewidth',1.5,'marker','x');
        
%             plot(hdgs, squeeze(confMean(m,c,D,:,1)), clr{c}{m}(1),'linestyle','-','linewidth',1.5,'marker',clr{c}{m}(2));
%             plot(hdgs, squeeze(confMean(m,c,D,:,2)), clr{c}{m}(1),'linestyle','--','linewidth',1.5,'marker','x');
%         
        else
            %h(m) = errorbar(hdgs, squeeze(pRight(m,c,D,:)), squeeze(pRightSE(m,c,D,:)), [clr{c}{m} '-']); hold on;
        end
        
        if m~=1
            title([modlabels{m} ', coh = ' num2str(cohs(c))]); 
        else
            title(modlabels{m})
        end
    legend('High Bet','Low Bet','Location','northwest');
    if m==length(mods), xlabel('heading angle (deg)'); end
    if m==2, ylabel('Mean RT (s)'); ylim([0.5 1.0]);
    else, ylim([0.5 0.8]); end
    changeAxesFontSize(gca,15,15);
    end
end

%{
if conftask==1
    
    for c = 1:length(cohs)
        subplot(2+(~isnan(RTmean(1,1,2))),length(cohs),c); box off; hold on;
        for m = 1:length(mods)     % m c d h
            if plotLogistic(m,c,D)
                h1(m) = plot(xVals,squeeze(yVals(m,c,D,:,1)),[clr{c}{m}(1) '-'],'linewidth',1.5); hold on;
                h2(m) = plot(xVals,squeeze(yVals(m,c,D,:,2)),[clr{c}{m}(1) '--'],'linewidth',1.5); hold on;
                
                %errorbar(hdgs, squeeze(pRight(m,c,D,:)), squeeze(pRightSE(m,c,D,:)), clr{c}{m});
            else
                %h(m) = errorbar(hdgs, squeeze(pRight(m,c,D,:)), squeeze(pRightSE(m,c,D,:)), [clr{c}{m} '-']); hold on;
            end
            ylim([0 1]);
            if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
        end
        legend(h1,'vestib','visual','comb','Location','northwest');
        xlabel('heading angle (deg)'); ylabel('proportion rightward choices');
        changeAxesFontSize(gca,15,15);
        
        
        subplot(2+(~isnan(RTmean(1,1,2))),length(cohs),c+length(cohs));  box off; hold on;
        for m = 1:length(mods)
            h1(m) = errorbar(hdgs, squeeze(confMean(m,c,D,:,1)), squeeze(confSE(m,c,D,:,1)), [clr{c}{m} '-'],'linewidth',1.5);
            h2(m) = errorbar(hdgs, squeeze(confMean(m,c,D,:,2)), squeeze(confSE(m,c,D,:,2)), [clr{c}{m} '--'],'linewidth',1.5);
            
            ylim([0 1.6]); hold on;
        end
        xlabel('heading angle (deg)'); ylabel('saccadic endpoint (''confidence'', %)');
        changeAxesFontSize(gca,15,15);
        
        if ~isnan(RTmean(1,1,2))
            subplot(3,length(cohs),c+length(cohs)*2); box off; hold on;
            for m = 1:length(mods)
                h1(m) = errorbar(hdgs, squeeze(RTmean(m,c,D,:,1)), squeeze(RTse(m,c,D,:,1)), [clr{c}{m} '-'],'linewidth',1.5); hold on;
                h2(m) = errorbar(hdgs, squeeze(RTmean(m,c,D,:,2)), squeeze(RTse(m,c,D,:,2)), [clr{c}{m} '--'],'linewidth',1.5); hold on;
                
            end
            xlabel('heading angle (deg)'); ylabel('RT (s)');
        end
        changeAxesFontSize(gca,15,15);
    end
    
elseif conftask==2
    
    for c = 1:length(cohs)
        
        subplot(1+(~isnan(RTmean(1,1,2))),length(cohs),c); box off; hold on;
        for m = 1:length(mods)     % m c d h
            if plotLogistic(m,c,D)
                h1(m) = plot(xVals,squeeze(yVals(m,c,D,:,1)),[clr{c}{m}(1) '-'],'linewidth',1.5); hold on;
                h2(m) = plot(xVals,squeeze(yVals(m,c,D,:,2)),[clr{c}{m}(1) '--'],'linewidth',1.5); hold on;
                
                errorbar(hdgs, squeeze(pRight(m,c,D,:)), squeeze(pRightSE(m,c,D,:)), clr{c}{m},'linewidth',1.5);
            else
                %h(m) = errorbar(hdgs, squeeze(pRight(m,c,D,:)), squeeze(pRightSE(m,c,D,:)), [clr{c}{m} '-']); hold on;
            end
            ylim([0 1]);
            if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
        end
        legend(h1,'vestib','visual','comb','Location','northwest');
        xlabel('heading angle (deg)'); ylabel('proportion rightward choices');
        changeAxesFontSize(gca,15,15);
        
        if ~isnan(RTmean(1,1,2))
            subplot(2,length(cohs),c+length(cohs)); box off; hold on;
            for m = 1:length(mods)
                h1(m) = errorbar(hdgs, squeeze(RTmean(m,c,D,:,1)), squeeze(RTse(m,c,D,:,1)), [clr{c}{m} '-'],'linewidth',1.5); hold on;
                h2(m) = errorbar(hdgs, squeeze(RTmean(m,c,D,:,2)), squeeze(RTse(m,c,D,:,2)), [clr{c}{m} '--'],'linewidth',1.5); hold on;
                
            end
            xlabel('heading angle (deg)'); ylabel('RT (s)');
        end
        changeAxesFontSize(gca,15,15);
    end
end
%}
