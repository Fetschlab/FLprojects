function fh=dots3DMP_plots_confchoice(parsedData,cohs,deltas,conftask)

if conftask == 0
    disp('Cannot plot conf vs choice for non-conf task, doh!')
    return
end

fsz = 15; % fontsize

ucond = [1 cohs(1); 2 cohs(1); 2 cohs(2); 3 cohs(1); 3 cohs(2)];
titles = {'Ves';'Vis (Low Coh)';'Vis (High Coh)';'Comb (Low Coh)';'Comb (High Coh)';'All'};

% first, for all trials irrespective of delta
D = length(deltas)+1; % (the extra column we made for pooling across deltas)
% OR select just delta=0:
% D = find(deltas==0);

         %ves %vis %comb
clr{1} = {'ko','mo','co'};
% clr{1} = {'ko','ro','bo'};
clr{2} = {'ko','ro','bo'};
clr{3} = {'ko','yo','go'};
fh=figure(111);
% set(gcf,'Color',[1 1 1],'Position',[300 1000 450+300*(length(cohs)-2) 200+150*(conftask>0)+150*RTtask],'PaperPositionMode','auto'); clf;
set(gcf,'Color',[1 1 1],'Position',[200 80 600 250],'PaperPositionMode','auto'); clf;
subplot(121);
for c = 1:size(ucond,1)
    h(c) = plot(squeeze(parsedData.pRight(ucond(c,1),ucond(c,2),D,:)), squeeze(parsedData.confMean(ucond(c,1),ucond(c,2),D,:)), [clr{ucond(c,2)}{ucond(c,1)} '-'],'linewidth',1.5); hold on;
%     text(0.4,0.95-(c-1)*0.05,titles{c},'color',clr{ucond(c,2)}{ucond(c,1)}(1),'fontsize',fsz)
end
axis([-0.05 1.05 0.3 0.8]);
xlabel('P(right)');
if conftask==1, ylabel('SaccEP');
elseif conftask==2, ylabel('P(High Bet)');
end
try changeAxesFontSize(gca,fsz,fsz); tidyaxes(gca,fsz); catch; disp('plot clean up skipped'); end

% now do deltas

if length(deltas)>1
    
    subplot(122);
    
    clr{1} = {'bs','cs','gs'};
    clr{2} = {'b^','c^','g^'};
    clr{3} = {'bo','co','go'};
    linstl = '-:';
    
    clear L;
    for d = 1:length(deltas)
        for c = 1:length(cohs)
            
            hd(c) = plot(squeeze(parsedData.pRight(3,c,d,:)), squeeze(parsedData.confMean(3,c,d,:)), [clr{c}{d} linstl(c)],'linewidth',1.5); hold on;
            
        end
    end
    
    axis([-0.05 1.05 0.3 0.8]);       
    xlabel('P(right)');
    
    try changeAxesFontSize(gca,fsz,fsz); tidyaxes(gca,fsz); catch; disp('plot clean up skipped'); end

end