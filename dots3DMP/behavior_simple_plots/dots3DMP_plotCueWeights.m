function dots3DMP_plotCueWeights(wves,wvesBoot,cohs,conftask)

if conftask==1
    xlabpos = 1.3; jit = 0.1; cohlabs = {'Low','High'}; ylabpos = 0.9;
else
    xlabpos = 0.4; jit = 0.02; cohlabs = cohs; ylabpos = 0.4;
end

figure(810); set(gcf,'Color',[1 1 1],'Position',[50 20 250 250],'PaperPositionMode','auto'); clf;
if nargin==5,subplot(121); end
hold on;

fs = 14;
h1 = gca;
erPred = errorbar(h1,cohs-jit,wves.choice.pred,std(wvesBoot.choice.pred),'k-o','LineWidth',1.5,'MarkerSize',5,'MarkerFaceColor','w'); 
erEmp = errorbar(h1,cohs+jit,wves.choice.emp,std(wvesBoot.choice.emp),'r-o','LineWidth',1.5,'MarkerSize',5,'MarkerFaceColor','w');    
set(erPred,'Color','k');
set(erEmp,'Color','r');


text(h1,xlabpos,ylabpos,'Choice (predicted)','color','k','fontsize',fs)
text(h1,xlabpos,ylabpos-0.1,'Choice (empirical)','color','r','fontsize',fs)

if conftask
    erConf = errorbar(h1,cohs,wves.conf.emp,std(wvesBoot.conf.emp),'b-o','LineWidth',1.5,'MarkerSize',5,'MarkerFaceColor','w');
    set(erConf,'Color','b')
    text(h1,xlabpos,ylabpos-0.2,'Conf (empirical)','color','b','fontsize',fs)
end

ylim(h1,[0 1]);
set(h1,'XTick',cohs,'Xlim',[cohs(1)-jit*3 cohs(end)+jit*3],'XTickLabel',cohlabs,'Ytick',0:0.125:1,'YTickLabel',{'0','','','','.5','','','','1'});
xlabel(h1,'Visual coherence'); ylabel(h1,'Vestibular weight');
changeAxesFontSize(h1,16,16); set(h1,'box','off'); tidyaxes(h1);


% h2 = axes;
% set(h2,'XLim',get(h1,'XLim'));
% set(h2,'Color','none');
% set(h2,'XTick',[]);
% set(h2,'YDir','reverse');
% set(h2,'YAxisLocation','Right');
% set(h2,'YLabel','Visual weight');

% plot thresholds

% if nargin==5
%     cols = 'krb';
%     subplot(122);
%     hold on;
%     
%     D = 4;
% 
%     % empirical thresholds
%     choice_thres = squeeze(gfit.choice.sigma(:,:,D));
%     conf_thres   = squeeze(gfit.conf.sigma(:,:,D));
%     
%     % normalize to ves
%     choice_thres = choice_thres ./ choice_thres(1,:);
%     conf_thres = conf_thres ./ conf_thres(1,:);
%     
%     temp = choice_thres.^2;
% 
%     choice_thres_pred = squeeze( (temp(1,:).*temp(2,:)) ./ (temp(1,:) + temp(2,:)) );
%     
%     for m=1:size(choice_thres,1)
%         plot(cohs,choice_thres(m,:),cols(m))
%     end
% end
    
    
    
