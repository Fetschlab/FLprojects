function confRT_distrs(data,mods,cohs,conftask,RTtask)

% confidence and RT distributions by correct/error
% for each condition
% overlapping analysis with RTquantiles...


% data input can be grouped subjects, or individual subject data as deemed
% necessary. seems more parsimonious to structure it this way than
% explicitly extract each subject here

if conftask==1
    coBins = 0:0.1:1; % assumes normalized!
    nco = nan(length(mods),length(cohs),length(coBins)-1);
    [nConfCorrect,nConfError] = deal(nco);
end

if RTtask
    rtBins = 0.5:0.1:2.5;
    nrt = nan(length(mods),length(cohs),length(rtBins)-1);
    [nRTCorrect,nRTError] = deal(nrt);
end
;


for m=1:length(mods)
    for c = 1:length(cohs)
        
        % zero delta only?
        J = data.modality==mods(m) & data.coherence==cohs(c) & data.delta==0;
        
        if RTtask
            rtData_Correct = data.RT(J & data.correct);
            rtData_Error   = data.RT(J & ~data.correct);
            
            nRTCorrect(m,c,:) = histcounts(rtData_Correct,rtBins,'Normalization','Probability');
            nRTError(m,c,:)   = histcounts(rtData_Error,rtBins,'Normalization','Probability');
        end
        
        if conftask==1
            confData_Correct = data.conf(J & data.correct);
            confData_Error   = data.conf(J & ~data.correct);
            
            nConfCorrect(m,c,:) = histcounts(confData_Correct,coBins,'Normalization','Probability');
            nConfError(m,c,:)   = histcounts(confData_Error,coBins,'Normalization','Probability');
        end
    end
    
   
end

% all trials
if conftask==1
    nConfCorrectALL = histcounts(data.conf(data.correct==1),coBins,'Normalization','Probability');
    nConfErrorALL   = histcounts(data.conf(data.correct==0),coBins,'Normalization','Probability');
end

if RTtask
    nRTCorrectALL = histcounts(data.RT(data.correct==1),rtBins,'Normalization','Probability');
    nRTErrorALL   = histcounts(data.RT(data.correct==0),rtBins,'Normalization','Probability');
end

% if 0
% figure;
% for m=1:length(mods)
%     for c = 1:length(cohs)
%         ax=subplot(length(mods),length(cohs),c+(m-1)*2);
%         if m==1 && c==2
%             plot(coBins(1:end-1),nConfCorrectALL,coBins(1:end-1),nConfErrorALL);
%             title('ALL trials')
%         else
%             plot(coBins(1:end-1),squeeze(nConfCorrect(m,c,:)),coBins(1:end-1),squeeze(nConfError(m,c,:)));
%             title(sprintf('m=%d,c=%d',m,c))
%         end
%     end
% end

% figure;
% for m=1:length(mods)
%     for c = 1:length(cohs)
%         ax=subplot(length(mods),length(cohs),c+(m-1)*2);
%         if m==1 && c==2
%             plot(rtBins(1:end-1),nRTCorrectALL,rtBins(1:end-1),nRTErrorALL);
%             title('ALL trials')
%         else
%             plot(rtBins(1:end-1),squeeze(nRTCorrect(m,c,:)),rtBins(1:end-1),squeeze(nRTError(m,c,:)));
%             title(sprintf('m=%d,c=%d',m,c))
%         end
%     end
% end
% end

% figure; 
% subplot(211); hold on;
% plot(coBins(1:end-1)+diff(coBins),nConfCorrectALL,'k','linew',1.5);
% plot(coBins(1:end-1)+diff(coBins),-nConfErrorALL,'r','linew',1.5);
% subplot(212); hold on;
% plot(rtBins(1:end-1)+diff(rtBins),nRTCorrectALL,'k','linew',1.5);
% plot(rtBins(1:end-1)+diff(rtBins),-nRTErrorALL,'r','linew',1.5);

%%
ucoh = unique(data.coherence);
ucond = [1 ucoh(1); 2 ucoh(1); 2 ucoh(2); 3 ucoh(1); 3 ucoh(2)];
titles = {'Ves';'Vis (Low Coh)';'Vis (High Coh)';'Comb (Low Coh)';'Comb (High Coh)';'All'};
mcols = {'Greys','Reds','Reds','Blues','Blues','Purples'};

fsz = 16;

subplotInd = [2 3 4 5 6 1];

uhdg  = unique(abs(data.heading));

figure;
for c = 1:size(ucond,1)+1
    
    cmap = flipud(cbrewer('seq',mcols{c},length(uhdg)*2));
    cmap = cmap(end-1:-2:1,:);
    
    subplot(3,2,subplotInd(c)); hold on
    axis([0.25 2.75 -0.1 1.1]);
    for h = 1:length(uhdg)
        if c==size(ucond,1)+1
            I = abs(data.heading)==uhdg(h);
        else
            I = abs(data.heading)==uhdg(h) & data.modality==ucond(c,1) & data.coherence==ucond(c,2);
        end
        
        scatter(data.RT(I & data.correct),data.conf(I & data.correct),18,cmap(h,:),'o','filled');
        scatter(data.RT(I & ~data.correct),data.conf(I & ~data.correct),18,cmap(h,:),'o');

        
    end
end
%  [b,bint,r,rint,stats] = regress(data.conf,[ones(size(data.RT)) data.RT]);


% TO DO
% make conf and RT on one figure, with one flipped 
% show all trials, with individual trials overlaid?
    
keyboard

