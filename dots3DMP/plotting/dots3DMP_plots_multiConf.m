function dots3DMP_plots_multiConf(parsedData,mods,cohs,deltas,hdgs,conftask,RTtask,splitPDW,isPDW)

% isPDW is 1 if the confGrouping is data.PDW or equivalent
% splitPDW is a logical that says whether to sub-split behavioural outcomes
% by PDW. this only makes sense if the initial split is not already PDW,
% but something else (e.g. available reward for high bet)

if isPDW
    if splitPDW
        disp('cannot sub-split by PDW when already split by PDW');
        splitPDW = 0; 
    end
end
xLab = sprintf('heading angle (%s)',char(176));

if conftask==1
    xt = -10:5:10;
elseif conftask==2
    xt = -12:6:12;
elseif conftask==0
    error('cannot plot confidence-based curves for non-confidence task, doh!');
end

%% first, for all trials irrespective of delta
D = length(deltas)+1; % (the extra column we made for pooling across deltas)
% OR select just delta=0:
% D = find(deltas==0);

modlabels = {'Ves','Vis','Comb'};
clrseqs = {'Greys','Reds','Blues'};

nConfGroups = length(parsedData.confGroups);

for m = 1:length(mods)
    if isPDW
        clr{mods(m)} = [0 0 1; 1 0 0; 0 0 0];
    else
        if splitPDW
            nn = ceil(nConfGroups/2);
        else
            nn = nConfGroups;
        end
        clr{mods(m)} = cbrewer('seq',clrseqs{m},nn+2);
        clr{mods(m)} = clr{mods(m)}(end-nn+1:end,:);
    end
end

if splitPDW || isPDW
    
% CHOICES i.e. pRight
figure(201+D);
set(gcf,'Color',[1 1 1],'Position',[300 500 180+300*(length(cohs)-1) 600],'PaperPositionMode','auto'); clf;
for c=1:length(cohs)
    for m=1:length(mods)
        subplot(length(mods),length(cohs),c+(m-1)*length(cohs))
        %         subplot(1,length(mods),m)
        if m==1 && c~=1, delete(gca); continue, end
       
        for nc=1:nConfGroups
            if splitPDW
                if nc>parsedData.confGroupSplit*2 % 1-targ
                    if splitPDW==4, lnstl = ':'; 
                    elseif splitPDW==3, lnstl = '-';
                    else, continue, 
                    end
                elseif nc>parsedData.confGroupSplit % high bet
                    if splitPDW==4, lnstl = ':';
                    elseif splitPDW==1, lnstl = '-';
                    else, continue, 
                    end
                else % low bet
                    if splitPDW==4, lnstl = ':';
                    elseif splitPDW==2, lnstl = '-';
                    else, continue, 
                    end
                end
            else
                lnstl = '-';
            end
            colind = mod(nc,parsedData.confGroupSplit); colind(colind==0) = parsedData.confGroupSplit;
            if parsedData.plotLogistic(m,c,D,nc)
                h(m) = plot(parsedData.xVals,squeeze(parsedData.yVals(m,c,D,:,nc)),'color',clr{mods(m)}(colind,:),'linestyle',lnstl,'linewidth',1.5); hold on;
                errorbar(hdgs, squeeze(parsedData.pRight(m,c,D,:,nc)), squeeze(parsedData.pRightSE(m,c,D,:,nc)), 'color',clr{mods(m)}(colind,:),'linewidth',1.5,'marker','o','linestyle','none');
            else
                h(m) = errorbar(hdgs, squeeze(parsedData.pRight(m,c,D,:,nc)), squeeze(parsedData.pRightSE(m,c,D,:,nc)), 'color',clr{mods(m)}(colind,:),'linewidth',1.5); hold on;
            end

            text(hdgs(1)+2,0.9-(nc-1)*0.1,sprintf('%.1f',squeeze(parsedData.meanRew(m,c,1,nc))),'color',clr{mods(m)}(colind,:),'fontweight','bold','fontsize',14)
        end
        set(gca,'xtick',xt);
        set(gca,'ytick',0:0.25:1,'yticklabel',{'0','.25','.5','.75','1'});
        ylim([0 1]);
        %             if m>1
        %                 %             if c==2&&m==3,keyboard,end
        %                 ht=title([modlabels{m} ', coh = ' num2str(cohs(c))]); set(ht,'color',clr{1}{m}(1));
        %             else
        %                 ht=title(modlabels{m}); %set(ht,'color',clr{1}{m}(1));
        %                 %hL=legend(h,'High Bet','Low Bet','Location','northwest','box','off');
        %             end
            
        
    if m==length(mods), xlabel(xLab); end
    if c==1, ylabel('P(right)'); end
    try changeAxesFontSize(gca,15,15); tidyaxes; catch; end
    end
end

% RT

if RTtask
    
figure(201+D+2);
set(gcf,'Color',[1 1 1],'Position',[300 500 180+300*(length(cohs)-1) 600],'PaperPositionMode','auto'); clf;
for c=1:length(cohs)
    for m=1:length(mods)
        subplot(length(mods),length(cohs),c+(m-1)*length(cohs))
        %         subplot(1,length(mods),m)
        if m==1 && c~=1, delete(gca); continue, end
        for nc=1:nConfGroups
            if splitPDW
                if nc>parsedData.confGroupSplit*2 % 1-targ
                    if splitPDW==4, lnstl = ':'; 
                    elseif splitPDW==3, lnstl = '-';
                    else, continue, 
                    end
                elseif nc>parsedData.confGroupSplit % high bet
                    if splitPDW==4, lnstl = ':';
                    elseif splitPDW==1, lnstl = '-';
                    else, continue, 
                    end
                else % low bet
                    if splitPDW==4, lnstl = ':';
                    elseif splitPDW==2, lnstl = '-';
                    else, continue, 
                    end
                end
            else
                lnstl = '-';
            end
            colind = mod(nc,parsedData.confGroupSplit); colind(colind==0) = parsedData.confGroupSplit;
            h(m) = errorbar(hdgs, squeeze(parsedData.RTmean(m,c,D,:,nc)), squeeze(parsedData.RTse(m,c,D,:,nc)), 'color',clr{mods(m)}(colind,:),'linestyle',lnstl,'linewidth',1.5,'marker','o'); hold on;
            
            if m==2, yy=1; else yy = 0.8; end
            text(hdgs(1)+2,yy-0.03-(nc-1)*0.03,sprintf('%.1f',squeeze(parsedData.meanRew(m,c,1,nc))),'color',clr{mods(m)}(colind,:),'fontweight','bold','fontsize',14)
            
        end
        
        set(gca,'xtick',xt);
        if conftask==1
            ylim([0.8 2]);
        else
            ylim([-0.3 0]+yy);
        end
        
        if m==length(mods), xlabel(xLab); end
        if c==1, ylabel('RT (s)'); end
        try changeAxesFontSize(gca,15,15); tidyaxes; catch; end
    end
end

end

end

% PDW, only if splitPDW==0

if conftask && ~splitPDW && ~isPDW
    
    figure(201+D+3);
    set(gcf,'Color',[1 1 1],'Position',[300 500 180+300*(length(cohs)-1) 600],'PaperPositionMode','auto'); clf;
    for c=1:length(cohs)
        for m=1:length(mods)
            subplot(length(mods),length(cohs),c+(m-1)*length(cohs))
            if m==1 && c~=1, delete(gca); continue, end
            for nc=1:nConfGroups
                colind = mod(nc,parsedData.confGroupSplit); colind(colind==0) = parsedData.confGroupSplit;
                h(nc) = errorbar(hdgs, squeeze(parsedData.confMean(m,c,D,:,nc)), squeeze(parsedData.confSE(m,c,D,:,nc)), 'color',clr{mods(m)}(colind,:),'linewidth',1.5); hold on;
                
                text(hdgs(1)+2,0.4-(nc-1)*0.1,sprintf('%.1f',squeeze(parsedData.meanRew(m,c,1,nc))),'color',clr{mods(m)}(colind,:),'fontweight','bold','fontsize',14)
            end
            set(gca,'xtick',xt);
            ylim([0 1])
                
            if m==length(mods), xlabel(xLab); end
            if c==1
                if conftask==2, ylabel('P(High Bet');
                else,           ylabel('SaccEP')
                end
            end
            try changeAxesFontSize(gca,15,15); tidyaxes; catch; end
        end
    end
    
end
