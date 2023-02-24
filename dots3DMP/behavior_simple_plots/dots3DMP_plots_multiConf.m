function fh = dots3DMP_plots_multiConf(parsedData,mods,cohs,deltas,hdgs,conftask,RTtask,splitPDW,splitcols)


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
%     if splitPDW
%         clr{mods(m)} = cbrewer('seq',clrseqs{m},ceil(nConfGroups/2));
%     else
%         clr{mods(m)} = cbrewer('seq',clrseqs{m},nConfGroups);
        clr{mods(m)} = cbrewer('qual','Dark2',nConfGroups);

%     end
end

% CHOICES i.e. pRight
fh(1) = figure(201+D);
set(gcf,'Color',[1 1 1],'Position',[300 500 180+300*(length(cohs)-1) 600],'PaperPositionMode','auto'); clf;
for c=1:length(cohs)
    for m=1:length(mods)
        subplot(length(mods),length(cohs),c+(m-1)*length(cohs))
        %         subplot(1,length(mods),m)
        if m==1 && c~=1, delete(gca); continue, end
        % plot all trials here?
        %             beta = [gfit.choice.mu(end,1,D,cc) gfit.choice.sigma(end,1,D,cc)];
        %             h(cc) = plot(parsedData.xVals, gfit.choice.func(beta,parsedData.xVals), 'color',[1 1 1]*0.5,'linestyle',clnstl(cc),'linewidth',1.5); hold on;
        %             errorbar(hdgs, squeeze(parsedData.pRight(end,1,D,:,cc)), squeeze(parsedData.pRightSE(end,1,D,:,cc)), 'color',[1 1 1]*0.5,'linewidth',1.5,'markerfacecolor',[1 1 1]*0.5);
        
        for nc=1:nConfGroups
%             if splitPDW
%                 if nc>parsedData.confGroupSplit*2
%                     continue
%                 elseif nc>parsedData.confGroupSplit
%                     lnstl = '-';
%                 else
%                     lnstl = '--';
%                 end
%             else
                lnstl = '-';
%             end
            colind = mod(nc,parsedData.confGroupSplit); colind(colind==0) = parsedData.confGroupSplit;
            if parsedData.plotLogistic(m,c,D,nc)
                try
                h(m) = plot(parsedData.xVals,squeeze(parsedData.yVals(m,c,D,:,nc)),'color',splitcols{m}(colind,:),'linestyle',lnstl,'linewidth',2); hold on;
                errorbar(hdgs, squeeze(parsedData.pRight(m,c,D,:,nc)), squeeze(parsedData.pRightSE(m,c,D,:,nc)), 'color',splitcols{m}(colind,:),'linestyle','none','linewidth',1.5);
                catch
                    keyboard
                end
            else
                h(m) = errorbar(hdgs, squeeze(parsedData.pRight(m,c,D,:,nc)), squeeze(parsedData.pRightSE(m,c,D,:,nc)), 'color',splitcols{m}(colind,:),'linewidth',2); hold on;
            end

            if m==1 && c==1
                text(xt(1)+1,1-(nc*0.1),parsedData.confGroupLabels{nc},'color',splitcols{m}(colind,:),'fontsize',14,'fontweight','bold');
            end

        end
        set(gca,'xtick',hdgs);
        set(gca,'ytick',0:0.25:1,'yticklabel',{'0','.25','.5','.75','1'});
        axis([xt(1) xt(end) 0 1]);
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
    
fh(2) = figure(201+D+2);
set(gcf,'Color',[1 1 1],'Position',[300 500 180+300*(length(cohs)-1) 600],'PaperPositionMode','auto'); clf;
for c=1:length(cohs)
    for m=1:length(mods)
        subplot(length(mods),length(cohs),c+(m-1)*length(cohs))
        %         subplot(1,length(mods),m)
        if m==1 && c~=1, delete(gca); continue, end
        for nc=1:nConfGroups
%             if splitPDW
%                 if nc>parsedData.confGroupSplit*2
%                     continue
%                 elseif nc>parsedData.confGroupSplit
%                     lnstl = '-';
%                 else
%                     lnstl = '--';
%                 end
%             else
                lnstl = '-';
%             end
            colind = mod(nc,parsedData.confGroupSplit); colind(colind==0) = parsedData.confGroupSplit;
            h(m) = errorbar(hdgs, squeeze(parsedData.RTmean(m,c,D,:,nc)), squeeze(parsedData.RTse(m,c,D,:,nc)), 'color',clr{mods(m)}(colind,:),'linewidth',2); hold on;
            set(gca,'xtick',hdgs);
            if conftask==1
                ylim([0.8 2]);
            else
                ylim([0.5 1]);
            end

            if m==1
                for ll=1:length(parsedData.confGroupLabels)
                    text(xt(2),0.9-(nc*0.1),parsedData.confGroupLabels{ll},'color',splitcols{m}(colind,:),'fontsize',12);
                end
            end
        end
        if m==length(mods), xlabel(xLab); end
        if c==1, ylabel('RT (s)'); end
        try changeAxesFontSize(gca,15,15); tidyaxes; catch; end
    end
end

end



% PDW, if splitPDW==0

if conftask && splitPDW==0
    
    cols = 'rbk';
    fh(3) = figure(201+D+3);
    set(gcf,'Color',[1 1 1],'Position',[300 500 180+300*(length(cohs)-1) 600],'PaperPositionMode','auto'); clf;
    for c=1:length(cohs)
        for m=1:length(mods)
            subplot(length(mods),length(cohs),c+(m-1)*length(cohs))
            if m==1 && c~=1, delete(gca); continue, end
            for nc=1:nConfGroups
                colind = mod(nc,parsedData.confGroupSplit); colind(colind==0) = parsedData.confGroupSplit;
                h(m) = errorbar(hdgs, squeeze(parsedData.confMean(m,c,D,:,nc)), squeeze(parsedData.confSE(m,c,D,:,nc)), 'color',clr{mods(m)}(colind,:),'linewidth',2); hold on;
                set(gca,'xtick',hdgs);
                if conftask==1
                    ylim([0.8 2]);
                else
                    ylim([0.25 1]);
                end
                %             if m>1
                %                 ht=title([modlabels{m} ', coh = ' num2str(cohs(c))]); set(ht,'color',clr{1}{m}(1));
                %             else
                %                 ht=title(modlabels{m}); %set(ht,'color',clr{1}{m}(1));
                %                 legend(h,'High Bet','Low Bet','Location','northwest','box','off');
                %             end
            end
            if m==length(mods), xlabel(xLab); end
            if c==1, ylabel('P(High Bet)'); end
            try changeAxesFontSize(gca,15,15); tidyaxes; catch; end
        end
    end
    
end
