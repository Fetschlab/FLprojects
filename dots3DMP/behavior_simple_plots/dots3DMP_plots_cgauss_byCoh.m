function fh=dots3DMP_plots_cgauss_byCoh(gfit,parsedData,mods,cohs,deltas,hdgs,conftask,RTtask)
% SJ 07-2021 converted to function for cleaner workspace

% first, for all trials irrespective of delta
% D = length(deltas)+1; % (the extra column we made for pooling across deltas)
% OR select just delta=0:
% D = find(deltas==0);
D = gfit.D;

if all(mods==1), cohs=1; end
fsz = 14; % fontsize

spRows = 1 + double(conftask>0) + double(RTtask);

modlabels = {'Vestibular','Visual','Combined'};
xLab = sprintf('heading angle (%s)',char(176));

choiceYlab = 'prop. rightward';
if length(mods)>1, cohlabs = {'Low Coherence','High Coherence'}; end


if conftask==1
    yLab = 'Confidence'; confYlims = [0.25 0.8];
    if all(mods==1), RTylims = [0.9 1.5];
    else,            RTylims = [0.9 1.7];
    end
    xt = -10:5:10;
    xtlab = xt;
    %     if sum(hdgs==0)
    %         xtlab = {'-10','-3.5','','0','','3.5','10'};
    %     else
    %         xtlab = {'-10','-3.5','-1.25','1.25','3.5','10'};
    %     end
    yt = 0:0.2:2; ytlab = {'0','','0.4','','0.8','','1.2','','1.6','','2.0'}; % for RT plots

elseif conftask==2
    yLab = 'P(High Bet)'; confYlims = [0.4 1];
    if all(mods==1), RTylims = [0.5 0.72];
    else,            RTylims = [0.5 1.0];
    end
    xt = hdgs;
    xtlab = {'-12','-6','-3','','0','','3','6','12'};
    yt = 0:0.1:1;  ytlab = {'0','','0.2','','0.4','','0.6','','0.8','','1.0'}; % for RT plots
    if length(mods)>1
        %         cohlabs = {sprintf('coh = %.1f',cohs(1)),sprintf('coh = %.1f',cohs(2))};
    end
end


%ves %vis %comb
% clr{1} = {'ko','mo','co'};
clr{1} = {'ko','ro','bo'};
clr{2} = {'ko','ro','bo'};
clr{3} = {'ko','ro','bo'};
clr{4} = clr{1};

fh(1)=figure(101+D);
set(fh(1),'Color',[1 1 1],'Position',[200 80 500 100+250*(double(conftask>0)+double(RTtask))],'PaperPositionMode','auto'); clf;
prettyfig(fh(1));

for c = 1:length(cohs)

    % *** choice ***
    try
        smplot(spRows,length(cohs),c,'axis','on');
    catch
        subplot(spRows,length(cohs),c)
    end
    box off; hold on;
    for m = 1:length(mods)     % m c d h
        beta = [gfit.choice.mu(m,c,D) gfit.choice.sigma(m,c,D)];
        h(m) = plot(parsedData.xVals, gfit.choice.func(beta,parsedData.xVals), [clr{c}{m}(1) '-'],'linewidth',2); hold on;
        errorbar(hdgs, squeeze(parsedData.pRight(m,c,D,:)), squeeze(parsedData.pRightSE(m,c,D,:)), clr{c}{m},'linewidth',2);
        if c==length(cohs), text(hdgs(1)+1,1.0-m*0.16,modlabels{m},'color',clr{c}{m}(1),'fontsize',fsz,'fontweight','bold'); end
        %         text(hdgs(end),0.5-m*0.12,sprintf('%+2.2f, %2.2f',beta(1),beta(2)),'color',clr{c}{m}(1),...
        %             'fontsize',fsz-6,'fontweight','bold','HorizontalAlignment','right');
    end
    if length(mods)>1; title(cohlabs{c}); end
    ylim([0 1]);
    set(gca,'xtick',xt,'xticklabel',xtlab);
    set(gca,'ytick',0:0.25:1,'yticklabel',{'0','','.5','','1'});
    if ~conftask && ~RTtask, xlabel(xLab); end
    if c==1, ylabel(choiceYlab); end
    try changeAxesFontSize(gca,fsz,fsz); tidyaxes(gca,fsz); end


    % *** conf ***
    if conftask
        try
            smplot(spRows,length(cohs),c+length(cohs),'axis','on');
        catch
            subplot(spRows,length(cohs),c+length(cohs))
        end
        box off; hold on;
        for m = 1:length(mods)
            beta = [gfit.conf.ampl(m,c,D) gfit.conf.mu(m,c,D) gfit.conf.sigma(m,c,D) gfit.conf.bsln(m,c,D)];
            h(m) = plot(parsedData.xVals, gfit.conf.func(beta,parsedData.xVals), [clr{c}{m}(1) '-'],'linewidth',2); hold on;
            errorbar(hdgs, squeeze(parsedData.confMean(m,c,D,:)), squeeze(parsedData.confSE(m,c,D,:)), clr{c}{m},'linewidth',2);
        end
        %         if length(mods)>1; title(cohlabs{1}); end
        set(gca,'xtick',xt,'xticklabel',xtlab);
        set(gca,'ytick',0:0.1:1,'yticklabel',{'0','','.2','','.4','','.6','','.8','','1.0'});

        ylim(confYlims);
        xlabel(xLab);
        if c==1, ylabel(yLab); end

        try, changeAxesFontSize(gca,fsz,fsz); tidyaxes(gca,fsz);  end

    end

    % *** RT ***
    if RTtask
        try
            smplot(spRows,length(cohs),c+length(cohs)*(2-double(conftask==0)),'axis','on');
        catch
            subplot(spRows,length(cohs),c+length(cohs)*(2-double(conftask==0)));
        end
        box off; hold on;
        for m = 1:length(mods)
            beta = [gfit.RT.ampl(m,c,D) gfit.RT.mu(m,c,D) gfit.RT.sigma(m,c,D) gfit.RT.bsln(m,c,D)];
            h(m) = plot(parsedData.xVals, gfit.RT.func(beta,parsedData.xVals), [clr{c}{m}(1) '-'],'linewidth',2); hold on;
            errorbar(hdgs, squeeze(parsedData.RTmean(m,c,D,:)), squeeze(parsedData.RTse(m,c,D,:)), clr{c}{m},'linewidth',2);
        end
        %         if length(mods)>1; title(cohlabs{1}); end
        set(gca,'xtick',xt,'xticklabel',xtlab);
        set(gca,'ytick',yt,'yticklabel',ytlab);
        ylim(RTylims)
        xlabel(xLab);
        if c==1, ylabel('RT (s)'); end
        try changeAxesFontSize(gca,fsz,fsz); tidyaxes(gca,fsz);  end

    end
end


%% now repeat, but just separate combined condition by delta

if length(deltas)>1 && any(mods==3)

    % clr{1} = {'bs','cs','gs'};
    % clr{2} = {'b^','c^','g^'};
    % clr{3} = {'bo','co','go'};
    % clr{4} = {'bd','cd','gd'};

    clr = {'bo','co','go'};
    if sum(deltas==0)==0
        clr = {'bo','go'};
    end


    % clr{1} = {'bs','gs'};
    % clr{2} = {'b^','g^'};
    % clr{3} = {'bo','go'};
    % clr{4} = {'bd','gd'};

    clear L;
    fh(2) = figure(208);
    set(fh(2),'Color',[1 1 1],'Position',[800 80 500 100+250*(double(conftask>0)+double(RTtask))],'PaperPositionMode','auto'); clf;
    prettyfig(fh(2));

    for c = 1:length(cohs)
        % *** choice ***
        try
            smplot(spRows,length(cohs),c,'axis','on');
        catch
            subplot(spRows,length(cohs),c)
        end
        box off; hold on;
        for d = 1:length(deltas)     % m c d h
            beta = [gfit.choice.mu(3,c,d) gfit.choice.sigma(3,c,d)];
            h(d) = plot(parsedData.xVals, gfit.choice.func(beta,parsedData.xVals), [clr{d}(1) '-'],'linewidth',2); hold on;
            errorbar(hdgs, squeeze(parsedData.pRight(3,c,d,:)), squeeze(parsedData.pRightSE(3,c,d,:)), clr{d},'linewidth',2);
            L{d} = sprintf('\\Delta=%+d',deltas(d));
            text(hdgs(1)+1,1.0-d*0.16,L{d},'color',clr{d}(1),'fontsize',fsz,'fontweight','bold');
        end
        if length(mods)>1; title(cohlabs{c}); end
        set(gca,'xtick',xt,'xticklabel',xtlab);
        set(gca,'ytick',0:0.25:1,'yticklabel',{'0','','.5','','1'});
        ylim([0 1]);
        %     lh=legend(h,L,'location','southeast'); set(lh,'box','off');
        %     xlabel(xLab);
        if c==1, ylabel(choiceYlab); end
        try changeAxesFontSize(gca,fsz,fsz); tidyaxes(gca,fsz);  end



        % *** conf ***
        if conftask
            try
                smplot(spRows,length(cohs),c+length(cohs),'axis','on');
            catch
                subplot(spRows,length(cohs),c+length(cohs))
            end
            box off; hold on;
            for d = 1:length(deltas)
                beta = [gfit.conf.ampl(3,c,d) gfit.conf.mu(3,c,d) gfit.conf.sigma(3,c,d) gfit.conf.bsln(3,c,d)];
                h(d) = plot(parsedData.xVals, gfit.conf.func(beta,parsedData.xVals), [clr{d}(1) '-'],'linewidth',2); hold on;
                errorbar(hdgs, squeeze(parsedData.confMean(3,c,d,:)), squeeze(parsedData.confSE(3,c,d,:)), clr{d},'linewidth',2);
                %             L{d} = sprintf('\\Delta=%+d',deltas(d));
            end
            %         if length(mods)>1; title(cohlabs{c}); end
            set(gca,'xtick',xt,'xticklabel',xtlab);
            set(gca,'ytick',0:0.1:1,'yticklabel',{'0','','.2','','.4','','.6','','.8','','1.0'});
            ylim(confYlims);
            xlabel(xLab);

            if c==1, ylabel(yLab); end
            try changeAxesFontSize(gca,fsz,fsz); tidyaxes(gca,fsz);  end

        end

        % *** RT ***
        if RTtask
            try
                smplot(spRows,length(cohs),c+length(cohs)*(2-double(conftask==0)),'axis','on');
            catch
                subplot(spRows,length(cohs),c+length(cohs)*(2-double(conftask==0)));
            end
            box off; hold on;
            for d = 1:length(deltas)
                beta = [gfit.RT.ampl(3,c,d) gfit.RT.mu(3,c,d) gfit.RT.sigma(3,c,d) gfit.RT.bsln(3,c,d)];
                h(d) = plot(parsedData.xVals, gfit.RT.func(beta,parsedData.xVals), [clr{d}(1) '-'],'linewidth',2); hold on;
                errorbar(hdgs, squeeze(parsedData.RTmean(3,c,d,:)), squeeze(parsedData.RTse(3,c,d,:)), clr{d},'linewidth',2);
                %             L{d} = sprintf('\\Delta=%+d',deltas(d));
            end
            %         if length(mods)>1; title(cohlabs{c}); end
            set(gca,'xtick',xt,'xticklabel',xtlab);
            set(gca,'ytick',yt,'yticklabel',ytlab);
            xlabel(xLab);
            if c==1, ylabel('RT (s)'); end
            ylim(RTylims);
            try changeAxesFontSize(gca,fsz,fsz); tidyaxes(gca,fsz);  end
        end

    end

    % suptitle('Combined condition: Cue Conflict')

end