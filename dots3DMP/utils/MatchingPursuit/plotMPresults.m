clear all
close all
load MPresults

% set frequency range of interest
freqMin = 0; freqMax = 65;
tempmin = abs(f-freqMin); tempmax = abs(f-freqMax);
freqRange = find(tempmin==min(tempmin)) : find(tempmax==min(tempmax));

% log-transform, or not
% rEnergy1(rEnergy1<1e-4) = 1e-4;
% rEnergy2(rEnergy2<1e-4) = 1e-4;
% rEnergy1 = log(rEnergy1);
% rEnergy2 = log(rEnergy2);


figure; set(gcf,'Position',[100,100,800,800]);
subplot(421);
plot(t,inputSignal(:,trialNum,1)); 
title(['Condition 1, trial: ' num2str(trialNum)]); xlabel('Time (s)');
axis tight

subplot(423);
plot(t,rSignal1,'k');
title('Reconstruction'); xlabel('Time (s)');
axis tight

subplot(223);
pcolor(t,f(freqRange),rEnergy1(freqRange,:)); shading interp;
xlabel('Time (s)'); ylabel('Frequency (Hz)');
axis tight; colorbar;

subplot(422);
plot(t,inputSignal(:,trialNum,2));
title(['Condition 2, trial: ' num2str(trialNum)]); xlabel('Time (s)');
axis tight

subplot(424);
plot(t,rSignal2,'k');
title('Reconstruction'); xlabel('Time (s)');
axis tight

subplot(224);
pcolor(t,f(freqRange),rEnergy2(freqRange,:)); shading interp;
xlabel('Time (s)'); ylabel('Frequency (Hz)');
axis tight; colorbar;


figure; set(gcf,'Position',[200,100,800,500]);
subplot(1,2,1); title('signal 1');
plot(f(freqRange),mean(rEnergy1(freqRange,:),2),'LineWidth',1.1);
xlabel('Frequency');
ylabel('Mean energy');subplot(1,2,2); title('signal 2');
plot(f(freqRange),mean(rEnergy2(freqRange,:),2),'LineWidth',1.1);
xlabel('Frequency');
ylabel('Mean energy');