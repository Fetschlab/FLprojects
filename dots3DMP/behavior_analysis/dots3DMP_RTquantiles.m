function [fh]=dots3DMP_RTquantiles(data,conftask,plotOption)

if conftask==0, error('Need a confidence task!'); end
if nargin<3, plotOption = 2; end

% 2 figures, each with 3x2 config (all mods, ves, visLo, visHigh, combLo,
% combHigh)

% plotOption == 0 - plot errors/low bet only
% plotOption == 1 - plot correct/high bet only
% plotOption == 2 - plot correct/error or high/low bet separately
% plotOption == -1 - plot all trials

% 1 - P(High Bet) as function of RT quantiles
  
% loop over 'conditions' indexed as [modality coherence], with the
% familiar manual kluge necessary to avoid invalid combinations:
ucoh = unique(data.coherence);
ucond = [1 ucoh(1); 2 ucoh(1); 2 ucoh(2); 3 ucoh(1); 3 ucoh(2)];
titles = {'Ves';'Vis (Low Coh)';'Vis (High Coh)';'Comb (Low Coh)';'Comb (High Coh)';'All'};

uhdg  = unique(abs(data.heading));

if conftask==1
    confdata = data.conf >= median(data.conf);
    errfun   = @(x,n) std(x) / sqrt(n);
    yLab = 'Sacc EP';
    yL = [0 1];
    nbins = 4; % number of RT quantiles
    xRange = [0.4 2.2];  % assume human for RT purposes
elseif conftask==2 
    confdata = data.PDW;
    errfun   = @(x,n) sqrt( (x.*(1-x)) ./ n);
    yLab = 'P(High Bet)';
    yL = [0 1];
    nbins = 5;
    xRange = [0.35 1.5];
end

for c = 1:size(ucond,1)+1 % the extra one is for all conditions pooled
    
 
    for h = 1:length(uhdg)
        if c==size(ucond,1)+1
            I = abs(data.heading)==uhdg(h);
        else
            I = abs(data.heading)==uhdg(h) & data.modality==ucond(c,1) & data.coherence==ucond(c,2);
        end
        I = I & ~data.oneTargConf;

        if plotOption==-1
            theseRT = data.RT(I);
            theseConf = confdata(I);
            theseCorr = data.correct(I);
            
            rtQ = [0 quantile(theseRT,nbins-1) inf]; 
            for q = 1:length(rtQ)-1
                J = theseRT>=rtQ(q) & theseRT<rtQ(q+1);
                X(c,h,q) = mean(theseRT(J));
                Y(c,h,q) = mean(theseConf(J));
                Yc(c,h,q) = mean(theseCorr(J));
                
                Ye(c,h,q) = errfun(Y(c,h,q),sum(J));
                Yce(c,h,q) = errfun(Yc(c,h,q),sum(J));
            end

        else
            
            % RTs for correct and error trials under each condition
            corrRT   = data.RT(I & (data.correct | data.heading==0));
            errRT    = data.RT(I & (~data.correct | data.heading==0));

            corrRT   = data.RT(I & data.correct);
            errRT    = data.RT(I & ~data.correct);

            % RTs for high and low bets under each condition
            highRT   = data.RT(I & confdata==1);
            lowRT    = data.RT(I & confdata==0);
            
            % PDW for correct and error under each condition
            corrConf = confdata(I & (data.correct | data.heading==0));
            errConf  = confdata(I & (~data.correct | data.heading==0));

            corrConf = confdata(I & data.correct);
            errConf  = confdata(I & ~data.correct);
            
            % accuracy for high and low bets under each condition
            highCorr  = data.correct(I & confdata==1);
            lowCorr   = data.correct(I & confdata==0);
            
            rtQ_corr  = [0 quantile(corrRT,nbins-1) inf];
            rtQ_err   = [0 quantile(errRT,nbins-1) inf];
            rtQ_high  = [0 quantile(highRT,nbins-1) inf];
            rtQ_low   = [0 quantile(lowRT,nbins-1) inf];
            
            for q = 1:length(rtQ_corr)-1
                J = corrRT>=rtQ_corr(q) & corrRT<rtQ_corr(q+1);
                X(c,h,q,1) = mean(corrRT(J));
                Y(c,h,q,1) = mean(corrConf(J));
                Ye(c,h,q,1) = errfun(Y(c,h,q,1),sum(J));
                
                J = errRT>=rtQ_err(q) & errRT<rtQ_err(q+1);
                X(c,h,q,2) = mean(errRT(J));
                Y(c,h,q,2) = mean(errConf(J));
                Ye(c,h,q,2) = errfun(Y(c,h,q,2),sum(J));

                J = highRT>=rtQ_high(q) & highRT<rtQ_high(q+1);
                Xc(c,h,q,1) = mean(highRT(J));
                Yc(c,h,q,1) = mean(highCorr(J));
                Yce(c,h,q,1) = errfun(Yc(c,h,q,1),sum(J));

                J = lowRT>=rtQ_low(q) & lowRT<rtQ_low(q+1);
                Xc(c,h,q,2) = mean(lowRT(J));
                Yc(c,h,q,2) = mean(lowCorr(J));
                Yce(c,h,q,2) = errfun(Yc(c,h,q,2),sum(J));

            end
            
            
        end
        
    end
end

if plotOption==2
    Ye = zeros(size(Ye));
    Yce = zeros(size(Yce));
end
                    
    
%% plot conf vs RT, for correct / error trials


subplotInd = [2 3 4 5 6 1];
mcols = {'Greys','Reds','Reds','Blues','Blues','Purples'};
fsz = 16;

fh(1)=figure(16);
set(gcf,'Color',[1 1 1],'Position',[200 200 950 950],'PaperPositionMode','auto');

for c = 1:size(ucond,1)+1 % the extra one is for all conditions pooled
    
    cmap = cbrewer('seq',mcols{c},length(uhdg)*2);
    cmap = cmap(length(uhdg)+1:end,:);
    subplot(3,2,subplotInd(c));
    
    clear g L
    for h = 1:length(uhdg)  
        
%         len = 0.2;
%         hdgVec = len .* [sind(uhdg(h)) cosd(uhdg(h))];
%         startPoint = [1.3 0.7];
%         xVec = startPoint(1)+[0 hdgVec(1)];
%         yVec = startPoint(2)+[0 hdgVec(2)];
            
        if plotOption==-1 % plot all trials
            g(h) = errorbar(squeeze(X(c,h,:)),squeeze(Y(c,h,:)),squeeze(Ye(c,h,:)),'color',cmap(h,:),'LineWidth', 2,...
                'LineStyle','-','Marker','o','MarkerSize',6,'MarkerFaceColor',cmap(h,:)); hold on;
        elseif plotOption==0 % plot error trials only
            if h<=3
                g(h) = errorbar(squeeze(X(c,h,:,2)),squeeze(Y(c,h,:,2)),squeeze(Ye(c,h,:,2)),'color',cmap(h,:),'LineWidth', 2,...
                    'LineStyle','-','Marker','o','MarkerSize',6,'MarkerFaceColor',cmap(h,:)); hold on;
            end
        elseif plotOption==1 % plot correct trials only
            g(h) = errorbar(squeeze(X(c,h,:,1)),squeeze(Y(c,h,:,1)),squeeze(Ye(c,h,:,1)),'color',cmap(h,:),'LineWidth', 2,...
                'LineStyle','-','Marker','o', 'MarkerSize',6,'MarkerFaceColor',cmap(h,:)); hold on;
        else % plot correct and errors, separately
            if h<=3 % 
                g(h) = errorbar(squeeze(X(c,h,:,2)),squeeze(Y(c,h,:,2)),squeeze(Ye(c,h,:,2)),'color',cmap(h,:),'LineWidth', 2,...
                    'LineStyle',':','Marker','o','MarkerSize',6,'MarkerFaceColor','w'); hold on;
            end
            k(h) = errorbar(squeeze(X(c,h,:,1)),squeeze(Y(c,h,:,1)),squeeze(Ye(c,h,:,1)),'color',cmap(h,:),'LineWidth', 2,...
                'LineStyle','-','Marker','o','MarkerSize',6,'MarkerFaceColor',cmap(h,:)); hold on;
        end
        
%         plot(xVec,yVec,'color',cmap(h,:))
%         text(xVec(2),yVec(2),num2str(uhdg(h)),'color',cmap(h,:),'fontsize',12,'horizo','center');

%         if c==size(ucond,1)+1
            plot(xRange(2)*0.8+[0.08 0.15],0.9-0.15*h*ones(1,2),'color',cmap(h,:),'linewidth',3);
            text(xRange(2)*0.8+0.25,0.9-0.15*h,sprintf('%.2g%s',uhdg(h),char(176)),'color','k','fontsize',fsz,'horizo','center');
%         end
       
    end
    

    xlim(xRange);
    ylim([0 1])
    
    if c<size(ucond,1)+1 
        if ucond(c,1)==3,xlabel('RT (s)');
        %else, set(gca,'xticklabel',[]);
        end
        if ucond(c,1)==2 && ucond(c,2)==ucoh(1)
            ylabel(yLab)
        end
    else
%         set(gca,'xticklabel',[]);
    end
    text(xRange(2)*0.8+0.25,0.9,'|hdg|','color',cmap(end,:),'fontsize',14,'horizo','center','fontweight','bold');

    if mod(c,2)==1
        set(gca,'yticklabel',[]);
    end
    changeAxesFontSize(gca,fsz,fsz); tidyaxes(gca,fsz); set(gca,'box','off');
    set(gca,'ytick',0:0.25:1,'yticklabel',{'0','','.5','','1'});
    set(gca,'xtick',0:0.25:2.25);
    if conftask==1,set(gca,'xticklabel',{'0','','.5','','1','','1.5','','2',''});
    end
    title(titles{c});
end
% sh=suptitle('Confidence-RT'); set(sh,'fontsize',fsz,'fontweight','bold');

%% repeat for accuracy vs RT, for high / low bets

fh(2)=figure(17);
set(gcf,'Color',[1 1 1],'Position',[200 200 950 950],'PaperPositionMode','auto');

for c = 1:size(ucond,1)+1 % the extra one is for all conditions pooled
    
    cmap = cbrewer('seq',mcols{c},length(uhdg)*2);
    cmap = cmap(length(uhdg)+1:end,:);
    subplot(3,2,subplotInd(c));
    
    clear g L
    for h = 1:length(uhdg)      
        
        if plotOption==-1 % plot all trials
            g(h) = errorbar(squeeze(X(c,h,:)),squeeze(Yc(c,h,:)),squeeze(Yce(c,h,:)),'color',cmap(h,:),'LineWidth', 2,...
                'LineStyle','-','Marker','o','MarkerSize',6,'MarkerFaceColor',cmap(h,:)); hold on;
        elseif plotOption==0 % plot low only
            if h<=3
                g(h) = errorbar(squeeze(Xc(c,h,:,2)),squeeze(Yc(c,h,:,2)),squeeze(Yce(c,h,:,2)),'color',cmap(h,:),'LineWidth', 2,...
                    'LineStyle','-','Marker','o','MarkerSize',6,'MarkerFaceColor',cmap(h,:)); hold on;
            end
        elseif plotOption==1 % plot high only
            g(h) = errorbar(squeeze(Xc(c,h,:,1)),squeeze(Yc(c,h,:,1)),squeeze(Yce(c,h,:,1)),'color',cmap(h,:),'LineWidth', 2,...
                'LineStyle','-','Marker','o', 'MarkerSize',6,'MarkerFaceColor',cmap(h,:)); hold on;
        else % plot high and low separately
            if h<=3
                g(h) = errorbar(squeeze(Xc(c,h,:,2)),squeeze(Yc(c,h,:,2)),squeeze(Yce(c,h,:,2)),'color',cmap(h,:),'LineWidth', 2,...
                    'LineStyle',':','Marker','o','MarkerSize',6,'MarkerFaceColor','w'); hold on;
            end
            k(h) = errorbar(squeeze(Xc(c,h,:,1)),squeeze(Yc(c,h,:,1)),squeeze(Yce(c,h,:,1)),'color',cmap(h,:),'LineWidth', 2,...
                'LineStyle','-','Marker','o','MarkerSize',6,'MarkerFaceColor',cmap(h,:)); hold on;
        end
%         if c==size(ucond,1)+1
            plot(xRange(2)*0.8+[0.08 0.15],0.9-0.15*h*ones(1,2),'color',cmap(h,:),'linewidth',3);
            text(xRange(2)*0.8+0.3,0.9-0.15*h,sprintf('%.2g%s',uhdg(h),char(176)),'color','k','fontsize',fsz,'horizo','center');
%         end
    end
    

    xlim(xRange);
    ylim([0 1])
    
    if c<size(ucond,1)+1 
        if ucond(c,1)==3,xlabel('RT (s)');
        %else, set(gca,'xticklabel',[]);
        end
        if ucond(c,1)==2 && ucond(c,2)==ucoh(1)
            ylabel('Accuracy')
        end
    else
%         set(gca,'xticklabel',[]);
    end
    text(xRange(2)*0.8+0.25,0.9,'|hdg|','color',cmap(end,:),'fontsize',14,'horizo','center','fontweight','bold');

    if mod(c,2)==1
        set(gca,'yticklabel',[]);
    end
    changeAxesFontSize(gca,fsz,fsz); tidyaxes(gca,fsz); set(gca,'box','off');
    set(gca,'ytick',0:0.25:1,'yticklabel',{'0','','.5','','1'});
    set(gca,'xtick',0:0.25:2.25);
    if conftask==1,set(gca,'xticklabel',{'0','','.5','','1','','1.5','','2',''});
    end
    title(titles{c});
end
% sh=suptitle('Accuracy-RT'); set(sh,'fontsize',fsz,'fontweight','bold');
 