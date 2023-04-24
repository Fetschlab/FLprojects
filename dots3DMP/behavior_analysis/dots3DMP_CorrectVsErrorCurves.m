function fh = dots3DMP_CorrectVsErrorCurves(data,conftask,RTtask,useAbsHdg)
% look at PDW and RT separately for correct and incorrect trials

if useAbsHdg, uhdg = unique(abs(data.heading));
else,         uhdg = unique(data.heading);
end

if conftask==1, yLab = 'Confidence';
elseif conftask==2, yLab = 'P(high bet)';
end

% ucoh = unique(data.coherence);
% ucond = [1 ucoh(1); 2 ucoh(1); 2 ucoh(2); 3 ucoh(1); 3 ucoh(2)];
% mcols = [0 0 0 ; 1 0 1 ; 1 0 0; 0 1 1; 0 0 1; .5 0 .5]; 

ucond = [1;2;3];
mcols = [0 0 0; 1 0 0; 0 0 1];

% nconds = size(ucond,1)+1; % include 'all conds' trace
nconds = size(ucond,1);

for c = 1:nconds+1 % the extra one is for all conditions pooled
    
    for h = 1:length(uhdg)
        if c==size(ucond,1)+1
            I = data.heading==uhdg(h);
        else
            I = data.heading==uhdg(h) & data.modality==ucond(c,1);% & data.coherence==ucond(c,2);
        end
        
        J = I & data.correct;
        nCorr(h,c) = sum(J);
        
        if conftask==1
            confCorr(h,c) = mean(data.conf(J));
            coffCorrSE(h,c) = std(data.conf(J))/sqrt(sum(J));
        elseif conftask==2
            confCorr(h,c) = mean(data.PDW(J));
        end
        
        K = I & ~data.correct;
        nErr(h,c) = sum(K);
        if conftask==1
            confErr(h,c) = mean(data.conf(K));
            coffErrSE(h,c) = std(data.conf(K))/sqrt(sum(K));
        elseif conftask==2
            confErr(h,c) = mean(data.PDW(K));
        end
        
        if RTtask
            RTCorr(h,c) = mean(data.RT(J));
            RTCorrSE(h,c) = std(data.RT(J))/sqrt(sum(J));
            
            RTErr(h,c) = mean(data.RT(K));
            RTErrSE(h,c) = std(data.RT(K))/sqrt(sum(K));
        end
    end
end

if conftask==2
    confCorrSE = sqrt( (confCorr.*(1-confCorr)) ./ nCorr );
    confErrSE = sqrt( (confErr.*(1-confErr)) ./ nErr );
else
    disp('not showing standard errors')
    confCorrSE = zeros(size(confCorr));
    confErrSE = zeros(size(confErr));
end

if RTtask
    xRange = [min([RTCorr(:);RTErr(:)]) max([RTCorr(:);RTErr(:)])] .* [0.95 1.05];
end


% PDW vs heading, correct trials
fh = figure('color','white','position',[300 300 400 400]);
subplot(221); hold on
for c = 1:nconds
    errorbar(uhdg,confCorr(:,c),confCorrSE(:,c),'color',mcols(c,:),'linestyle','-','linew',1.5,'marker','o');
end
ylabel(yLab);
axis([min(uhdg)-1 max(uhdg)+1 0.25 1])
set(gca,'xtick',uhdg);
try changeAxesFontSize(gca,12,12); tidyaxes; catch; end
box off;
title('Correct')

% PDW vs heading, error trials
subplot(223); hold on
for c = 1:nconds
    errorbar(uhdg(1:3),confErr(1:3,c),confErrSE(1:3,c),'color',mcols(c,:),'linestyle','-','linew',1.5,'marker','o');
end
xlabel(sprintf('heading angle (%s)',char(176))); 
ylabel(yLab);
axis([min(uhdg)-1 max(uhdg)+1 0.25 1])
set(gca,'xtick',uhdg);
try changeAxesFontSize(gca,12,12); tidyaxes; catch; end
box off;
title('Error')

% RT vs heading, correct trials
if RTtask
subplot(222); hold on
for c = 1:nconds
    errorbar(RTCorr(:,c),confCorr(:,c),confCorrSE(:,c),'color',mcols(c,:),'linestyle','-','linew',1.5,'marker','o');
    if c==2
        for h=1:length(uhdg)
            text(RTCorr(h,c)+0.02,confCorr(h,c)+0.02,sprintf('%g',uhdg(h)),'color',mcols(c,:),'fontsize',12)
        end
    end
end
axis([xRange 0 1])
try changeAxesFontSize(gca,12,12); tidyaxes; catch; end
title('Correct')

% RT vs heading, error trials
subplot(224); hold on
for c = 1:nconds
    errorbar(RTErr(1:3,c),confErr(1:3,c),confErrSE(1:3,c),'color',mcols(c,:),'linestyle','-','linew',1.5,'marker','o');
end
xlabel('RT (s)')
axis([xRange 0 1])
try changeAxesFontSize(gca,12,12); tidyaxes; catch; end
end
title('Error')