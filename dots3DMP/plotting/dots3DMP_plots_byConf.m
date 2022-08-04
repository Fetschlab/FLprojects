function dots3DMP_plots_byConf(parsedData,mods,cohs,deltas,hdgs,conftask,RTtask)

% few options for how to plot this...
% same way as for regular fits i.e. one plot with all three modalities for
% given coh
% or, since we are mostly interested in differences within condition, let's
% plot ves,vis,comb x cohs on separate subplots
% and we can do a separate figure for deltas again

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

         %ves %vis %comb
clr{1} = {'ko','ro','bo'};
clr{2} = {'kd','rd','bd'};
fclr{1} = 'krb';
fclr{2} = 'www';
clnstl = '-:'; % high low conf line style



% CHOICES i.e. pRight
figure(201+D);
set(gcf,'Color',[1 1 1],'Position',[300 500 180+300*(length(cohs)-1) 600],'PaperPositionMode','auto'); clf;
for c=1:length(cohs)
    for m=1:length(mods)
        subplot(length(mods),length(cohs),c+(m-1)*length(cohs))
%         subplot(1,length(mods),m)
        if m==1 && c~=1
            delete(gca);
            % plot all trials here?
%             beta = [gfit.choice.mu(end,1,D,cc) gfit.choice.sigma(end,1,D,cc)];
%             h(cc) = plot(parsedData.xVals, gfit.choice.func(beta,parsedData.xVals), 'color',[1 1 1]*0.5,'linestyle',clnstl(cc),'linewidth',1.5); hold on;
%             errorbar(hdgs, squeeze(parsedData.pRight(end,1,D,:,cc)), squeeze(parsedData.pRightSE(end,1,D,:,cc)), 'color',[1 1 1]*0.5,'linewidth',1.5,'markerfacecolor',[1 1 1]*0.5);
        else
            
        for cc=1:2
            if parsedData.plotLogistic(m,c,D,cc)
                h(m) = plot(parsedData.xVals,squeeze(parsedData.yVals(m,c,D,:,cc)),clr{cc}{m}(1),'linestyle',clnstl(cc),'linewidth',1.5); hold on;
                errorbar(hdgs, squeeze(parsedData.pRight(m,c,D,:,cc)), squeeze(parsedData.pRightSE(m,c,D,:,cc)), clr{cc}{m},'linewidth',1.5);
            else
                h(m) = errorbar(hdgs, squeeze(parsedData.pRight(m,c,D,:,cc)), squeeze(parsedData.pRightSE(m,c,D,:,cc)), 'color',clr{cc}{m}(1),'linestyle',clnstl(cc),'linewidth',1.5,'markerfacecolor',fclr{cc}(m));
            end
            
            
        end
        set(gca,'xtick',xt);
        set(gca,'ytick',0:0.25:1,'yticklabel',{'0','.25','.5','.75','1'});
        ylim([0 1]);
        if m>1
            ht=title([modlabels{m} ', coh = ' num2str(cohs(c))]); set(ht,'color',clr{1}{m}(1)); 
        else
            ht=title(modlabels{m}); set(ht,'color',clr{1}{m}(1)); 
            hL=legend(h,'High Bet','Low Bet','Location','northwest','box','off');
        end
        
        end
        
    if m==length(mods), xlabel(xLab); end
    if c==1, ylabel('P(right)'); end
    try changeAxesFontSize(gca,15,15); tidyaxes; catch; end
    end
end


% RT

if RTtask
    
figure(201+D+1);
set(gcf,'Color',[1 1 1],'Position',[300 500 180+300*(length(cohs)-1) 600],'PaperPositionMode','auto'); clf;
for c=1:length(cohs)
    for m=1:length(mods)
        subplot(length(mods),length(cohs),c+(m-1)*length(cohs))
%         subplot(1,length(mods),m)
        if m==1 && c~=1, delete(gca); continue, end
        
        for cc=1:2
            h(cc) = errorbar(hdgs, squeeze(parsedData.RTmean(m,c,D,:,cc)), squeeze(parsedData.RTse(m,c,D,:,cc)), 'color',clr{cc}{m}(1),'linewidth',1.5,'linestyle',clnstl(cc)); hold on;
        end
        set(gca,'xtick',xt);
        if conftask==1
            ylim([0.8 2]);
        else    
            ylim([0.5 1]);
        end
        if m>1
            ht=title([modlabels{m} ', coh = ' num2str(cohs(c))]); set(ht,'color',clr{1}{m}(1)); 
        else
            ht=title(modlabels{m}); set(ht,'color',clr{1}{m}(1)); 
            legend(h,'High Bet','Low Bet','Location','northwest','box','off');
        end
    if m==length(mods), xlabel(xLab); end
    if c==1, ylabel('RT (s)'); end
    try changeAxesFontSize(gca,15,15); tidyaxes; catch; end
    end
end

end

% CONFIDENCE
figure(201+D+2);
set(gcf,'Color',[1 1 1],'Position',[300 500 950+300*(length(cohs)-2) 800],'PaperPositionMode','auto'); clf;
for c=1:length(cohs)
    for m=1:length(mods)
        subplot(length(mods),length(cohs),c+(m-1)*length(cohs)); hold on
        if m==1 && c~=1, delete(gca); continue, end
        
        for cc=1:2
            errorbar(hdgs, squeeze(parsedData.confMean(m,c,D,:,cc)), squeeze(parsedData.confSE(m,c,D,:,cc)), clr{cc}{m}(1),'linewidth',1.5,'linestyle',clnstl(cc));
        end
           
        ylim([0 1]);
        if m~=1
            title([modlabels{m} ', coh = ' num2str(cohs(c))]); 
        else
            title(modlabels{m})
        end
    if m==length(mods), xlabel('heading angle (deg)'); end
    if m==2
        if conftask==1, ylabel('P(High Bet'); 
        else,           ylabel('SEP')
        end
    end
    try changeAxesFontSize(gca,15,15); catch; end
    end
end

%{

%% now separate by delta

if length(deltas)>1
    
clr{1} = {'bs','cs','gs'};
clr{2} = {'b^','c^','g^'};
clr{3} = {'bo','co','go'};

% CHOICES i.e. pRight
clear L;
figure(208);
set(gcf,'Color',[1 1 1],'Position',[50 20 950+300*(length(cohs)-2) 800],'PaperPositionMode','auto'); clf;
for c = 1:length(cohs)
    for d = 1:length(deltas)
        subplot(length(cohs),length(deltas),d+length(deltas)*(c-1)); box off; hold on;
        for cc=1:2
            beta = [gfit.choice.mu(3,c,d,cc) gfit.choice.sigma(3,c,d,cc)];            
            h(cc) = plot(parsedData.xVals, gfit.choice.func(beta,parsedData.xVals), [clr{cc}{d}(1)],'linewidth',1.5); hold on;
            errorbar(hdgs, squeeze(parsedData.pRight(3,c,d,:,cc)), squeeze(parsedData.pRightSE(3,c,d,:,cc)), clr{cc}{d},'linewidth',1.5);
        end
           
        ylim([0 1]);
        legend(h,'High Bet','Low Bet','Location','northwest');
        xlabel('heading angle (deg)'); ylabel('P(right)');
        try changeAxesFontSize(gca,15,15); catch; end
    end
end

% RT data for individual confidences and deltas is too noisy?
%{
if RTtask
figure(209);
set(gcf,'Color',[1 1 1],'Position',[50 20 950+300*(length(cohs)-2) 800],'PaperPositionMode','auto'); clf;
for c = 1:length(cohs)
    for d = 1:length(deltas)
        subplot(length(cohs),length(deltas),d+length(deltas)*(c-1)); box off; hold on;
        for cc=1:2
            beta = [gfit.RT.ampl(3,c,d,cc) gfit.RT.mu(3,c,d,cc) gfit.RT.sigma(3,c,d,cc) gfit.RT.bsln(3,c,d,cc)];
            h(cc) = plot(parsedData.xVals, gfit.RT.func(beta,parsedData.xVals), clr{cc}{d}(1),'linewidth',1.5); hold on;
            errorbar(hdgs, squeeze(parsedData.RTmean(3,c,d,:,cc)), squeeze(parsedData.RTmean(3,c,d,:,cc)), clr{cc}{d},'linewidth',1.5);
        end
           
%         ylim([0 1]);
        legend(h,'High Bet','Low Bet','Location','northwest');
        xlabel('heading angle (deg)'); ylabel('P(right)');
        try changeAxesFontSize(gca,15,15); catch; end
    end
end

end
%}

end
%}
