function [logOddsMapR, logOddsMapL, logOddsCorrMap, tAxis, vAxis] = dots3DMP_makeLogOddsCorrMap(hdgs,k,B,sigma,maxdur,plotflag)

% CF wrote it ca 2013
% essentially a wrapper for runFPonSeries aka FP4, by RK 2007
% also recreates the figures in Kiani & Shadlen 2009

% 9-2019 modified for 3DMP, heading is the difficulty variable, momentary
% evidence is k*sin(hdg) not k*C, and duration is longer

uvect = k*sind(hdgs); % vector of mean drift rates
dT = 1; % time step (ms)
buf = 4; % buffer for bound crossings
tAxis = (0:dT:maxdur)'; % time axis

y0 = 0;
Bup = ones(1,length(tAxis))*(B+buf); % upper bound
Blo = ones(1,length(tAxis))*-(B+buf); % lower bound

% use Fokker-Planck equation to propagate the probability density of DV;
% xtDist is the probability as a function of x (DV) and t (time)
[~,~,~,~,~,~,~,~,xtDist] = runFPonSeries(uvect,tAxis,Bup,Blo,y0,sigma);


% remake Fig. 4a-b of Kiani and Shadlen 2009
vAxis = linspace(-(B+buf), B+buf, size(xtDist{1},1));
mapLeft = 0; mapRight = 0;
% for each coherence, the "map" is the probability density conditioned on
% the stimulus being either "right" or "left" (Fig. 4A)

% this marginalizes over motion strength for a given direction (see text below Equation in paper)
for c = 1:length(hdgs)
    if hdgs(c)>0
        mapRight = mapRight + xtDist{c};
    else  % careful if using coh=0
        mapLeft = mapLeft + xtDist{c};
    end
end
mapRight = mapRight / sum(hdgs>0); 
mapLeft = mapLeft / sum(hdgs<0);

% clip PDFs at some minimum, for plotting purposes
mapRightClipped = mapRight; mapRightClipped(mapRightClipped<1e-6) = 1e-6;
mapLeftClipped = mapLeft; mapLeftClipped(mapLeftClipped<1e-6) = 1e-6;

n = 100; % set n to 100+ for Roozbeh's smooth plots, lower for faster plotting
if plotflag==2
    % log-transform P to enable log scale for contourf
    logPmap = (log10(mapRightClipped)+6)/6;
    figure; set(gcf, 'Color', [1 1 1], 'Position', [100 100 450 350], 'PaperPositionMode', 'auto'); hold on;
    [~,h] = contourf(tAxis,vAxis,logPmap,n); caxis([0 1]); % log transformed such that 1e-6 is 0 and 1 is 1.
    colorbar('YTick',[0 1/6 2/6 3/6 4/6 5/6 1],'YTickLabel',{'1e-6';'1e-5';'1e-4';'1e-3';'1e-2';'1e-1';'1'}); 
    set(h,'LineColor','none'); xlabel('time (ms)'); ylabel('DV'); title('Prob of DV for rightward motion');
    colormap(jet);

%     logPmap = (log10(mapLeftClipped)+6)/6;
%     figure; set(gcf, 'Color', [1 1 1], 'Position', [100 100 450 350], 'PaperPositionMode', 'auto'); hold on;
%     [~,h] = contourf(tAxis,dvAxis,logPmap,n); caxis([0 1]); % log transformed such that 1e-6 is 0 and 1 is 1.
%     colorbar('YTick',[0 1/6 2/6 3/6 4/6 5/6 1],'YTickLabel',{'1e-6';'1e-5';'1e-4';'1e-3';'1e-2';'1e-1';'1'}); 
%     set(h,'LineColor','none'); xlabel('time (ms)'); ylabel('DV');  title('PDF of DV for leftward motion');
end

% log odds is log(P/(1-P)), which in this case is Right/Left or vice versa
logOddsMapR = log(mapRight./mapLeft); 
logOddsMapL = log(mapLeft./mapRight);
logOddsCorrMap = [logOddsMapL(1:(size(vAxis,2)+1)/2,:) ; logOddsMapR((size(vAxis,2)+1)/2+1:end,:)];
if plotflag
    figure; set(gcf, 'Color', [1 1 1], 'Position', [100 100 450 350], 'PaperPositionMode', 'auto'); hold on;
    [~,h] = contourf(tAxis,vAxis,logOddsCorrMap,n);
    colorbar; set(h,'LineColor','none');
%     figure; [c,h] = contourf(tAxis,vAxis,logOddsMapL,n);
%     figure; [c,h] = contourf(tAxis,vAxis,logOddsMapR,n);
    caxis([0 5]);
% % %     [~,hh] = contour(tAxis,vAxis,logOddsCorrMap,[theta -theta]);
% % %     set(hh,'LineColor','k','LineWidth',2);
    colorbar('YTick',0:1:5); 
    xlabel('Time (ms)'); ylabel('DV');
    set(gca,'XTick',0:round(maxdur/5):maxdur,'YTick',-B:round(B/2):B,'XTickLabel',0:round(maxdur/5):maxdur,'TickDir','out');
%     changeAxesFontSize(gca,24,24);
    title('Log odds correct');
    colormap(jet);
end

