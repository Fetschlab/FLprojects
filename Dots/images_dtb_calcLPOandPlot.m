function logOddsCorrMap = images_dtb_calcLPOandPlot(R,P)

% Calculate log posterior odds of a correct response (Eq. 3 in Kiani et al.
% 2014), and plots it (plus a few things) depending on plotflag

% Formerly a part of images_dtb_2D.m, and requires P and R structs from
% that function.

% CF started it circa 2018
% SJ modified 02-2023 to remove reassigning logOdds map to P inside function
% SJ 06-2023 re-implemented 2-D drift rate (i.e. for time-varying drift)

I = mean(R.drift,2)>=0; % In case drift is signed, calculate only for positives,
                % then the kluge with separate marginals (Pxt's) can be done elsewhere
% unlike 1D code, here we'll marginalize over drift in one step               
odds = (squeeze(sum(P.up.distr_loser(I,:,:),1)) ./ size(R.drift(I),1)) ./ ...
       (squeeze(sum(P.lo.distr_loser(I,:,:),1)) ./ size(R.drift(I),1));
odds(odds<1) = 1; % fix some stray negatives/zeros (what about infs?)
odds(isinf(odds)) = nan;
logOddsCorrMap = log(odds)';

if R.plotflag

    % (1) plot choice and RT
    figure(111); set(gcf,'Color',[1 1 1],'Position',[600 600 450 700],'PaperPositionMode','auto'); clf;
    subplot(3,1,1); plot(mean(R.drift,2),P.up.p,'o-'); title('Prob correct bound crossed before tmax') % meaningless w signed drift
    subplot(3,1,2); plot(mean(R.drift,2),P.up.p./(P.up.p+P.lo.p),'o-'); title('Relative prob of corr vs. incorr bound crossed (Pcorr)'); % actually pRight with signed drift
    subplot(3,1,3); plot(mean(R.drift,2),P.up.mean_t,'o-'); title('mean RT'); xlabel('drift rate');

    % remake Fig. 5c-d of Kiani et al 2014 (steps analogous to makeLogOddsCorrMap_*)
    n = 100; % set n to 100+ for smooth plots, lower for faster plotting
    
    % (2) first an example PDF
    c = round(size(R.drift(I),1)/2) - 1 + sum(I==0); % pick an intermediate drift rate, or make a loop to see all of them
    q = 30; % exponent for log cutoff (redefine zero as 10^-q, for better plots)
    Pmap = squeeze(P.up.distr_loser(c,:,:))';
    Pmap(Pmap<10^-q) = 10^-q;
    logPmap = (log10(Pmap)+q)/q;

    figure(c*100); set(gcf, 'Color', [1 1 1], 'Position', [200 300 500 850/R.plotflag], 'PaperPositionMode', 'auto');
    if R.plotflag==1; subplot(2,1,1); end
    [~,h] = contourf(R.t,P.y,logPmap,n); colormap(jet);
    caxis([0 1]); % because we log transformed such that 10^-q is 0 and 1 is 1
    colorbar('YTick',0:.2:1,'YTickLabel',{['10^-^{' num2str(q) '}']; ['10^-^{' num2str(q*.8) '}']; ['10^-^{' num2str(q*.6) '}']; ['10^-^{' num2str(q*.4) '}']; ['10^-^{' num2str(q*.2) '}']; '1'}); 
    set(gca,'XLim',[0 R.t(end)],'XTick',0:0.5:floor(R.t(end)*2)/2,'TickDir','out');
    set(gca,'YLim',[-4*R.Bup 0],'YTick',floor(-4*R.Bup):round(R.Bup*2)/2:0);
    set(h,'LineColor','none');
    xlabel('Time (s)'); ylabel('Accumulated evidence of losing accumulator');
    title('Probability density of losing accumulator');
    changeAxesFontSize(gca,15,15);
    if R.plotflag==2
        ylim([-2.5 0]);
        xlabel([]);ylabel([]);title([]);
        changeAxesFontSize(gca,24,24);
        export_fig('2dacc_PDF','-eps');
    end
        
    % (3) then the log odds corr map
    if R.plotflag==1
        subplot(2,1,2);
    else
        figure(c*1000); set(gcf, 'Color', [1 1 1], 'Position', [700 400 450 700/R.plotflag], 'PaperPositionMode', 'auto');
    end
    [~,h] = contourf(R.t,P.y,logOddsCorrMap,n); % colormap(parula);
    caxis([0 3]);
    colorbar('YTick',0:0.5:3); 
    set(gca,'XLim',[0 R.t(end)],'XTick',0:0.5:floor(R.t(end)*2)/2,'TickDir','out');
    set(gca,'YLim',[-4*R.Bup 0],'YTick',floor(-4*R.Bup):round(R.Bup*2)/2:0); 
    set(h,'LineColor','none');
    xlabel('Time (s)'); ylabel('Accumulated evidence of losing accumulator');
    title('Log odds correct vs. state of losing accumulator');
    changeAxesFontSize(gca,15,15);
    
    if R.plotflag==2
        ylim([-2.5 0])
        xlabel([]);ylabel([]);title([]);
        changeAxesFontSize(gca,24,24);
        export_fig('2dacc_logoddscorr','-eps')
    end

end
