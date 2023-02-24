%% the X-shape
% confidence for correct and error trials as a function of heading angle

data.corr = (data.heading>0 & data.choice==2) | (data.heading<0 & data.choice==1);
data.corr(data.heading==0) = 999;

for h = 1:length(hdgs)
    I = data.heading==hdgs(h) & data.corr;
    confCorr(h) = mean(data.conf(I));
    coffCorrSE(h) = std(data.conf(I))/sqrt(sum(I));
    J = data.heading==hdgs(h) & ~data.corr;
    confErr(h) = mean(data.conf(J));
    coffErrSE(h) = std(data.conf(J))/sqrt(sum(J));
end
figure('position',[300 300 800 300]);subplot(131); plot(hdgs,confCorr,'b-o',hdgs,confErr,'r-o','linew',1.5);
xlabel('heading (deg)');
changeAxesFontSize(gca,14,14);
ylabel('confidence');
legend('corrects','errors','Location','Southeast');

clear confCorr confErr
ushdgs = hdgs(hdgs>0);
for h = 1:length(ushdgs)
    I = abs(data.heading)==ushdgs(h) & data.corr;
    confCorr(h) = mean(data.conf(I));
    coffCorrSE(h) = std(data.conf(I))/sqrt(sum(I));
    J = abs(data.heading)==ushdgs(h) & ~data.corr;
    confErr(h) = mean(data.conf(J));
    coffErrSE(h) = std(data.conf(J))/sqrt(sum(J));
    
    RTCorr(h) = mean(data.RT(I));
    RTCorrSE(h) = std(data.RT(I))/sqrt(sum(I));
    
    RTErr(h) = mean(data.RT(J));
    RTErrSE(h) = std(data.RT(J))/sqrt(sum(J));
    
end
subplot(132); plot(ushdgs,confCorr,'b-o',ushdgs,confErr,'r-o','linew',1.5);
xlabel('|heading| (deg)');
changeAxesFontSize(gca,14,14);

subplot(133); plot(RTCorr,confCorr,'b-o',RTErr,confErr,'r-o','linew',1.5);
xlabel('RT');
legend('correct','error')
changeAxesFontSize(gca,14,14);
ylabel('confidence)');
legend('corrects','errors','Location','Northeast');

%% compare choice shifts vs conf shifts

exportfigs = 0;

% first just raw biases
n = 1;
for c=1:length(cohs)
    for d = 1:length(deltas)
        choiceBias(n) = muPMF(3,c,d);
        choiceBiasSE(n) = muPMFse(3,c,d);
        confBias(n) = muConf(3,c,d);
        confBiasSE(n) = muConfse(3,c,d);
        n = n+1;
    end
end
figure(20); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
[hsym,hxe,hye] = errorbar2(choiceBias, confBias, choiceBiasSE, confBiasSE, 'o', 2);
set(hsym,'MarkerSize',10,'MarkerFaceColor','w','Color','k');
hold on; plot([-2 2],[-2 2],'k--','LineWidth',2); axis square;
xlim([-2 2]); ylim([-2 2]);
set(gca,'Xtick',-2:1:2,'Ytick',-2:1:2);
xlabel('Choice shift (deg)'); ylabel('Confidence shift (deg)');
changeAxesFontSize(gca,20,20); set(gca,'box','off')
if exportfigs; export_fig('biasVsConf','-eps'); end



%% calculate predicted/empirical weights (requires dots3DMP_fit_cgauss to be run first)

for c = 1:length(cohs)      % m c d
    wvesPred(c) = (1/sigmaPMF(1,1,D)^2) / ((1/sigmaPMF(1,1,D)^2) + (1/sigmaPMF(2,c,D)^2));
                     % m c d
    actual(1) = (muPMF(3,c,1)-muPMF(3,c,2)+(deltas(1)/2)) / deltas(1);
    actual(2) = (muPMF(3,c,3)-muPMF(3,c,2)+(deltas(3)/2)) / deltas(3);    
    wvesEmp(c) = mean(actual);

    actual(1) = (muConf(3,c,1)-muConf(3,c,2)+(deltas(1)/2)) / deltas(1);
    actual(2) = (muConf(3,c,3)-muConf(3,c,2)+(deltas(3)/2)) / deltas(3);        
    wvesConfBased(c) = mean(actual);
end

nboots = 100;
tic
[muPMFboot,sigmaPMFboot,wvesEmpboot,wvesPredboot] = dots3DMP_cgauss_bootstrap(data,mods,cohs,deltas,nboots);
toc

% looky there!
figure(810); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
er = errorbar(cohs,wvesPred,std(wvesPredboot),std(wvesPredboot),'k-o','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','k'); hold on;
set(er,'Color','k');
ylim([0 1.1]);
set(gca,'XTick',cohs,'Xlim',[cohs(1)-0.04 cohs(end)+0.04],'XTickLabel',cohs*100,'Ytick',0:0.2:1);
xlabel('Visual coherence (%)'); ylabel('Vestibular weight');
changeAxesFontSize(gca,20,20); set(gca,'box','off');
l = legend('Predicted (from single-cues)'); legend('boxoff');
set(l,'Position',[0.2484    0.9250    0.7444    0.0812]);
if exportfigs; export_fig('weights1','-eps'); end

er = errorbar(cohs,wvesEmp,std(wvesEmpboot),std(wvesEmpboot),'m-o','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','m');    
set(er,'Color','m')
l = legend('Predicted', 'Empirical');
set(l,'Position',[0.2317    0.8492    0.7444    0.1516]);
if exportfigs; export_fig('weights2','-eps'); end

er = errorbar(cohs,wvesConfBased,std(wvesConfBasedboot),std(wvesConfBasedboot),'c-o','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','c');    
set(er,'Color','c')
l = legend('Predicted', 'Empirical', 'From confidence curves');
set(l,'Position',[0.3095    0.7890    0.7444    0.2219]);
if exportfigs; export_fig('weights3','-eps'); end


%% does overall conf (or RT) depend on delta?

data.conf = data.PDW;
for c = 1:length(cohs)
    I = data.modality==3 & data.coherence==cohs(c) & data.delta==0;
    I = I & abs(data.heading<2);
    confDzero(c) = mean(data.conf(I));
    confDzeroSE(c) = std(data.conf(I))/sqrt(sum(I));
    rtDzero(c) = mean(data.RT(I));
    rtDzeroSE(c) = std(data.RT(I))/sqrt(sum(I));

    J = data.modality==3 & data.coherence==cohs(c) & data.delta~=0;
    J = J & abs(data.heading<2);
    confDnonzero(c) = mean(data.conf(J));
    confDnonzeroSE(c) = std(data.conf(J))/sqrt(sum(J));
    rtDnonzero(c) = mean(data.RT(J));
    rtDnonzeroSE(c) = std(data.RT(J))/sqrt(sum(J));

    [~,pvalConf(c)] = ttest2(data.conf(I),data.conf(J));
    [~,pvalRT(c)] = ttest2(data.RT(I),data.RT(J));
end

barx = 1:length(cohs);
bary = [confDzero ; confDnonzero]';
errlow = [confDzeroSE ; confDnonzeroSE]';
errhigh = errlow;

figure(809); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
bar(barx,bary); hold on;
er = errorbar([barx-0.14;barx+0.14]',bary,errlow,errhigh,'k.');    
ylim([0 1]);
set(gca,'XTickLabel',[10 50 90],'Ytick',0:0.2:1);
xlabel('Visual coherence (%)'); ylabel('Avg confidence rating');
legend('no conflict', 'conflict'); legend('boxoff');
changeAxesFontSize(gca,20,20); set(gca,'box','off');
if exportfigs; export_fig('confDeltaZeroNonzero','-eps'); end

bary = [rtDzero ; rtDnonzero]';
errlow = [rtDzeroSE ; rtDnonzeroSE]';
errhigh = errlow;

figure(810); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
bar(barx,bary); hold on;
er = errorbar([barx-0.14;barx+0.14]',bary,errlow,errhigh,'k.');    
set(gca,'XTickLabel',[10 50 90]);
xlabel('Visual coherence (%)'); ylabel('Mean RT (s)');
% legend('no conflict', 'conflict'); legend('boxoff');
changeAxesFontSize(gca,20,20); set(gca,'box','off');
if exportfigs; export_fig('RTDeltaZeroNonzero','-eps'); end




%% relationship bw conf and weights

% for each coh, bin trials by conf (try median split first), then calc
% weights separately for low and high conf trials

% this logic doesn't work, because confidence in the combined condition per
% se would not map onto weights even under the hypothesis

% BUT we can do it piecewise:

% high coh, left side of the graph, delta +3,
I = data.modality==3 & ismember(data.coherence,cohs(2:end)) & data.heading<=hdgs(3) & data.delta==deltas(end);
% when conf is lower, Pright should be higher (anticorrelated)
confLow = data.conf(I) < median(data.conf(I));
confHigh = data.conf(I) >= median(data.conf(I));
PrightLow = sum(confLow & data.choice(I)==2) / sum(confLow);
PrightHigh = sum(confHigh & data.choice(I)==2) / sum(confHigh);

clr = {[255 180 0]./255, [175 238 238]./255, [255 140 255]./255};
figure(21); set(gcf,'Color',[1 1 1],'Position',[50 20 250 200],'PaperPositionMode','auto'); clf;
h = plot([1 2], [PrightLow PrightHigh], 'ok-', 'MarkerSize',10, 'MarkerFaceColor','k');
set(h,'Color', clr{3}, 'MarkerFaceColor',clr{3});
xlim([0.5 2.5]); ylim([0 1]);
set(gca,'Xtick',[1 2],'XTickLabel',{'Low Conf','High Conf'},'Ytick',0:0.2:1);
ylabel('Proportion Rightward');
changeAxesFontSize(gca,16,16); set(gca,'box','off');
hold on;
if exportfigs; export_fig('medianSplit1','-eps'); end 

%% low coh, left side of the graph, delta -3, 
I = data.modality==3 & data.coherence==cohs(1) & data.heading<=hdgs(3) & data.delta==deltas(1);
% when conf lower, Pright should be higher (anticorrelated)
confLow = data.conf(I) < median(data.conf(I));
confHigh = data.conf(I) >= median(data.conf(I));
PrightLow = sum(confLow & data.choice(I)==2) / sum(confLow);
PrightHigh = sum(confHigh & data.choice(I)==2) / sum(confHigh);

h = plot([1 2], [PrightLow PrightHigh], 'o-', 'MarkerSize', 10);
set(h,'Color', clr{3}, 'MarkerFaceColor', clr{3});
if exportfigs; export_fig('medianSplit2','-eps'); end


%% high coh, right side of the graph, delta -3, 
I = data.modality==3 & ismember(data.coherence,cohs(2:end)) & data.heading>=hdgs(end-2) & data.delta==deltas(1);
% when conf is lower, Pright should be lower (correlated)
confLow = data.conf(I) < median(data.conf(I));
confHigh = data.conf(I) >= median(data.conf(I));
PrightLow = sum(confLow & data.choice(I)==2) / sum(confLow);
PrightHigh = sum(confHigh & data.choice(I)==2) / sum(confHigh);

h = plot([1 2], [PrightLow PrightHigh], 'ob-', 'MarkerSize',10, 'MarkerFaceColor','b');
set(h,'Color', clr{1}, 'MarkerFaceColor',clr{1});
if exportfigs; export_fig('medianSplit3','-eps'); end


%% low coh, right side of the graph, delta +3, 
I = data.modality==3 & data.coherence==cohs(1) & data.heading>=hdgs(end-2) & data.delta==deltas(end);
% when conf lower, Pright should be lower (correlated)
confLow = data.conf(I) < median(data.conf(I));
confHigh = data.conf(I) >= median(data.conf(I));
PrightLow = sum(confLow & data.choice(I)==2) / sum(confLow);
PrightHigh = sum(confHigh & data.choice(I)==2) / sum(confHigh);

h = plot([1 2], [PrightLow PrightHigh], 'ob-', 'MarkerSize',10, 'MarkerFaceColor','b');
set(h,'Color', clr{1}, 'MarkerFaceColor',clr{1});
if exportfigs; export_fig('medianSplit4','-eps'); end





%% try sliding window across trials (or session by session)






%% for RT, try kiani 2014 analysis:
% plot conf as a function of RT quantile, separately for abs(hdg)

sim=0;

% FOR SIM ONLY: exclude capped values (maxDur+Tnd, often 2300)
if sim
    dataBackup = data;
    removethese = data.RT==mode(data.RT);
    fnames = fieldnames(data);
    for F = 1:length(fnames)
        eval(['data.' fnames{F} '(removethese) = [];']);
    end
end


clear X Y
uhdg = unique(abs(data.heading));

% also loop over 'conditions' indexed as [modality coherence], with the
% familiar manual kluge necessary to avoid invalid combinations:
ucoh = unique(data.coherence);
ucond = [1 ucoh(1); 2 ucoh(1); 2 ucoh(2); 3 ucoh(1); 3 ucoh(2)];
subplotInd = [1 3 4 5 6];
titles = {'Ves-Low';'Vis-lo';'Vis-hi';'Comb-lo';'Comb-hi'};

for c = 1:size(ucond,1)+1 % the extra one is for all conditions pooled
    for h = 1:length(uhdg)
        if c==size(ucond,1)+1
            I = abs(data.heading)==uhdg(h) ;%& data.corr==0;
        else
            I = abs(data.heading)==uhdg(h) & data.modality==ucond(c,1) & data.coherence==ucond(c,2);%& data.corr==0;
        end
        theseRT = data.RT(I);
        theseConf = data.PDW(I);
%         rtQ = [0 quantile(theseRT,3) inf]; % try quartiles
        rtQ = [0 quantile(theseRT,4) inf]; % or quintiles
        for q = 1:length(rtQ)-1
            J = theseRT>=rtQ(q) & theseRT<rtQ(q+1);
            X(c,h,q) = mean(theseRT(J));
            Y(c,h,q) = mean(theseConf(J));
            Ye(c,h,q) = std(theseConf(J)) / sqrt(sum(J));
        end
    end

    if sim 
        Ltxt = {num2str(uhdg(1)), num2str(uhdg(2)), num2str(uhdg(3)), num2str(uhdg(4)), num2str(uhdg(5))};
        clr = {'bo-','rs-','c^-','gd-','kv-'};
    else
%         Ltxt = {num2str(uhdg(1)), num2str(uhdg(2)), num2str(uhdg(3))};
        Ltxt = {num2str(1.25), num2str(3.5), num2str(10)};
        clr = {'bo-','rs-','gd-'};
    end

    if c==size(ucond,1)+1
        figure(18+sim);
        set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
    else    
        figure(16+sim);
        set(gcf,'Color',[1 1 1],'Position',[200 200 360*2 320*3],'PaperPositionMode','auto');
        subplot(3,2,subplotInd(c));
    end

    clear g L
    % (loop over h to allow animating figures one heading at a time)
    for h = 1:length(uhdg)       
%         g(h) = plot(squeeze(X(c,h,:)),squeeze(Y(c,h,:)),clr{h},'LineWidth', 2); hold on;
        g(h) = errorbar(squeeze(X(c,h,:)),squeeze(Y(c,h,:)),squeeze(Ye(c,h,:)),clr{h},'LineWidth', 2); hold on;
        
        set(g(h),'MarkerSize',10,'MarkerFaceColor',clr{h}(1));
        xlim([0.5 2.5]);
        ylim([0 1.1]);
        % set(gca,'Xtick',-2:1:2,'Ytick',-2:1:2);
        xlabel('Reaction time (s)'); ylabel('Confidence');
        changeAxesFontSize(gca,20,20); set(gca,'box','off');
        L{h} = Ltxt{h};
        if c==1
            l = legend(g,L,'Location',[0.55 0.8 .4 .05],'Orientation','Horizontal');
            legend('boxoff')
            text(3.75+sim*0.10,0.7,'Heading (deg)','Fontsize',14)
        end
        if c==size(ucond,1)+1
            l = legend(g,L,'Location','South','Orientation','Horizontal');
            legend('boxoff');
%             export_fig(['kianiFig2-' num2str(h)],'-eps');
        else
            title(titles{c});
        end
    end
    
end

if sim
    data = dataBackup;
end

%% calculate subject's performance ('score') based on accuracy and metacog accuracy

% what we don't want is simply to reward confidence scaling with heading
% angle, or coherence
% rather, within a coh+hdg, how well does confidence predict accuracy?

corr = double((data.heading>0 & data.choice==2) | (data.heading<0 & data.choice==1));
corr(data.heading==0) = NaN;
q=1;
for m = 1:length(mods)
    for c = 1:length(cohs)
        for h = 1:length(hdgs)
            J = data.modality==mods(m) & data.coherence==cohs(c) & data.heading==hdgs(h); % all trials irrespective of delta
            confCorr = mean(data.conf(J & corr==1));
            confErr = mean(data.conf(J & corr==0));
%             diff(m,c,h) = confCorr - confErr;
            diff(q) = confCorr - confErr;

            cc = corrcoef(corr(J),data.conf(J));
%             phi(m,c,h) = cc(1,2);
            phi(q) = cc(1,2);
            
            q = q+1;
        end
    end
end



meanphi = nanmean(phi);
bestphi = 0.25; % hypothetical best performance
worstphi = 0; % not sure about this

totalPcorr = nansum(corr)/sum(data.heading~=0);
bestPcorr = 0.9; % hypothetical best performance
worstPcorr = 0.5;

scorePcorr = (min([totalPcorr bestPcorr]) - worstPcorr) / (bestPcorr-worstPcorr);
scorePhi = (min([meanphi bestphi]) - worstphi) / (bestphi-worstphi);

score =  scorePcorr*0.5 + scorePhi*0.5
payout = score*8 + 12



