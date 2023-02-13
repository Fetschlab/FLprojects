function fh = dots3DMP_plots_cgauss_byConf(gfit,parsedData,mods,cohs,deltas,hdgs,conftask,RTtask,D)

% few options for how to plot this...
% same way as for regular fits i.e. one plot with all three modalities for
% given coh
% or, since we are mostly interested in differences within condition, let's
% plot ves,vis,comb x cohs on separate subplots
% and we can do a separate figure for deltas again

xLab = sprintf('heading angle (%s)',char(176));

if conftask==2
    xtlab  = {'-12','-6','-3','','0','','3','6','12'}; 
    yt  = 0:0.25:1;
    ytl = {'0','.25','.5','.75','1'};
elseif conftask==1
    if sum(hdgs==0), xtlab = {'-10','-3.5','','0','','3.5','10'};
    else,            xtlab = {'-10','-3.5','-1.25','1.25','3.5','10'};
    end
    yt = 0:0.5:2;
    ytl = yt;
elseif conftask==0
    error('cannot plot confidence-based curves for non-confidence task, doh!');
end

fsz = 20;

%% first, for all trials irrespective of delta

% SJ add D as option, 11/21/2022
if nargin < 9 || isempty(D)
    D = length(deltas)+1; % (the extra column we made for pooling across deltas)
    % OR select just delta=0:
%     D = find(deltas==0);
end


modlabels = {'Ves','Vis','Comb'};

         %ves %vis %comb
clr{1} = {'ko','ro','bo'};
clr{2} = {'kd','rd','bd'};
fclr{1} = 'krb';
fclr{2} = 'www';
clnstl = {'-','--'}; % high low conf line style


% CHOICES i.e. pRight
fh(1) = figure(201+D);
set(gcf,'Color',[1 1 1],'Position',[300 500 250+300*(length(mods)-1) 300],'PaperPositionMode','auto'); clf;
for c=1:length(cohs)
    for m=1:length(mods)
        subplot(1,length(mods),m)

        for cc=1:2
            beta = [gfit.choice.mu(m,c,D,cc) gfit.choice.sigma(m,c,D,cc)];
            h(cc) = plot(parsedData.xVals, gfit.choice.func(beta,parsedData.xVals), clr{cc}{m}(1),'linestyle',clnstl{cc},'linewidth',2); hold on;
            errorbar(hdgs, squeeze(parsedData.pRight(m,c,D,:,cc)), squeeze(parsedData.pRightSE(m,c,D,:,cc)), clr{cc}{m},'linewidth',2,'markerfacecolor',fclr{cc}(m));
        end
        set(gca,'xtick',hdgs,'xticklabel',xtlab);
        set(gca,'ytick',yt,'yticklabel',ytl);
        ylim([0 1]);
        if m>1
            %             if c==2&&m==3,keyboard,end
            ht=title([modlabels{m} ', coh = ' num2str(cohs(c))]); set(ht,'color',clr{1}{m}(1));
        else
            ht=title(modlabels{m}); set(ht,'color',clr{1}{m}(1));
            %             hL=legend(h,'High Bet','Low Bet','Location','northwest','box','off','fontsize',fsz);
        end

        if m==ceil(length(mods)/2), xlabel(xLab); end
        if m==1, ylabel('P(right)'); end
        try changeAxesFontSize(gca,fsz,fsz); tidyaxes; catch; end
    end
end


% RT

if RTtask
    
fh(2) = figure(201+D+1);
set(gcf,'Color',[1 1 1],'Position',[300 500 250+300*(length(mods)-1) 300],'PaperPositionMode','auto'); clf;
for c=1:length(cohs)
    for m=1:length(mods)
%         subplot(length(mods),length(cohs),c+(m-1)*length(cohs))
        subplot(1,length(mods),m)
        %if m==1 && c~=1, delete(gca); continue, end
        
        for cc=1:2
            beta = [gfit.RT.ampl(m,c,D,cc) gfit.RT.mu(m,c,D,cc) gfit.RT.sigma(m,c,D,cc) gfit.RT.bsln(m,c,D,cc)];
            h(cc) = plot(parsedData.xVals, gfit.RT.func(beta,parsedData.xVals), clr{cc}{m}(1),'linestyle',clnstl{cc},'linewidth',2); hold on;
            errorbar(hdgs, squeeze(parsedData.RTmean(m,c,D,:,cc)), squeeze(parsedData.RTse(m,c,D,:,cc)), clr{cc}{m},'linewidth',2,'markerfacecolor',fclr{cc}(m));
        end
        set(gca,'xtick',hdgs,'xticklabel',xtlab);
        set(gca,'ytick',0:0.1:2,'yticklabel',{'0','','.2','','.4','','.6','','.8','','1','','1.2','','1.4','','1.6','','1.8','','2'});

        if conftask==1
            ylim([0.8 2]);
        else   
            if mods(m)==2, ylim([0.7 1.1])
            else, ylim([0.5 0.85]); end
        end
        if m>1
            ht=title([modlabels{m} ', coh = ' num2str(cohs(c))]); set(ht,'color',clr{1}{m}(1)); 
        else
            ht=title(modlabels{m}); set(ht,'color',clr{1}{m}(1)); 
            legend(h,'High','Low','Location','northwest','box','off','fontsize',fsz);
        end
    
    if m==ceil(length(mods)/2), xlabel(xLab); end
    if m==1, ylabel('RT'); end
    try changeAxesFontSize(gca,fsz,fsz); tidyaxes; catch; end

    end
end

end

%{
% CONFIDENCE (SANITY CHECK)
figure(201+D+2);
set(gcf,'Color',[1 1 1],'Position',[300 500 950+300*(length(cohs)-2) 800],'PaperPositionMode','auto'); clf;
for c=1:length(cohs)
    for m=1:length(mods)
        subplot(length(mods),length(cohs),c+(m-1)*length(cohs)); hold on
        if m==1 && c~=1, delete(gca); continue, end
        
        for cc=1:2
            errorbar(hdgs, squeeze(parsedData.confMean(m,c,D,:,cc)), squeeze(parsedData.confSE(m,c,D,:,cc)), clr{c}{m},'linewidth',1.5);
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
%}

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
