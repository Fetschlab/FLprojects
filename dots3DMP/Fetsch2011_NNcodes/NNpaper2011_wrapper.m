% Using the dataset from Fetsch et al. 2011 Nature Neuroscience, to
% recreate figures/analyses

% SJ updated 06-2021

% some figures look slightly different from figures in paper...why?
% likelihood calculations (Fig. 4) seem to be off

% CF 11-2024: fixed a few things, prob related to a version mismatch with 
% dots3DMP_likelihood_decoding.m.
% Fig 4 (single-trial likelihoods) still seems off, could be a bug
% somewhere, and congruency section throws an error because tuning_curves
% is missing a 5th dim...
% but the basic components are here if anyone wants to update it some day
 

%%
clear; clc
% load the data (change current folder as needed)
load Fetsch_et_al_NatNeuro_2011.mat

% struct "data" with following fields, each with 79029 rows (one for each
% trial, concatenated across sessions):

% filename: string format is m=monkey number (m18 is monkey Y, I think,
%           so m24 is monkey W), c=cell number, r = 'run' (block) number
% heading: heading angle in degrees
% coherence: visual stimulus strength, 16% ('low') or 60% ('high');
%            arbitrarily assigned 16 for vestib trials
% modality: aka stimulus condition, 1=vestib, 2=visual, 3=combined
% delta: confict angle (deg), positive means visual to the right of vestib
% choice: the monkey's saccadic decision, 1=rightward, 0=leftward
% spikes: raster array where each column is a 1-ms bin containing 1 for a spike and 0 for no spike;
%         length is 2200 because it includes 100 ms before and 100 ms after the 2-s stimulus epoch
% spikeRates: spike rate for each trial, typically counted in a 1-s bin centered on the peak velocity
%             of the stimulus, which turns out to be a bit later than t=1100 in the spikes matrix,
%             something like spikes(:,636:1636)(?). I recall that for a few cells this window was 
%             manually adjusted based on idiosyncratic firing rate dynamics, the details of
%             which are lost to time I'm afraid.


%%
mods   = unique(data.modality);
cohs   = unique(data.coherence);
deltas = unique(data.delta);
hdgs   = unique(data.heading);

% select monkey (optional), [] = select both, 'W' = m18, 'Y' = m24
monkey = []; % [], 'W','Y'

if ~isempty(monkey)
    switch monkey
        case 'W', monkID = 18;
        case 'Y', monkID = 24;
        otherwise
            error('no monkey with that ID');
    end
    removethese = ~startsWith(data.filename,['m' num2str(monkID)]);
    fnames = fieldnames(data);
    for F = 1:length(fnames)
        if strcmp(fnames(F), 'spikes')
            data.spikes(removethese,:) = [];
        else
            eval(['data.' fnames{F} '(removethese) = [];']);
        end
    end
else
    
end
%% Psychometric curves (Fig 1)

% logistic fits
parsedData = dots3DMP_parseData_func(data,mods,cohs,deltas,hdgs);
dots3DMP_plots_NN(parsedData,mods,cohs,deltas,hdgs);

% cum Gaussian fits... use gfit results for the weights and thresholds
gfit = dots3DMP_fit_cgauss_NN(data,mods,cohs,deltas);
dots3DMP_plots_cgauss_NN(gfit,parsedData,mods,cohs,deltas,hdgs)

%% weights and psycophysical thresholds (~Fig 2)
% NOTE: Fig. 2 in paper is plotted separately for each monkey

[wvesEmp,wvesPred] = dots3DMP_wgts_thres_NN(gfit.muPMF,gfit.sigmaPMF,cohs,deltas);

% bootstrapping for errorbars (resample data with replacement nboots times)
nboots = 100;
[muPMFboot,sigmaPMFboot,wvesEmpboot,wvesPredboot] = dots3DMP_cgauss_bootstrap_NN(data,mods,cohs,deltas,nboots);
dots3DMP_plot_wgts_bootstrap(wvesPred,wvesEmp,wvesEmpboot,wvesPredboot,sigmaPMFboot,gfit,cohs);

%% example MSTd neuron (Fig. 3)

[ufile,~,data.unitnum] = unique(data.filename,'stable');

clear monkUnit
monkUnit(startsWith(ufile,'m18'),1) = 'W'; % 48 units
monkUnit(startsWith(ufile,'m24'),1) = 'Y'; % 60 units

% calculate mean firing rate for each heading and condition
[meanFRs,semFRs] = dots3DMP_neuron_tuning(data,mods,cohs,deltas,hdgs); % all units

unit = 8; % 8 seems to be the one in the paper! 3 also looks nice
dots3DMP_plot_neuron_tuning(meanFRs(:,:,:,:,unit),semFRs(:,:,:,:,unit),cohs,hdgs);

%% Likelihood based decoding (Fig. 4)

step = 0.1; % heading interpolation (deg steps)
numtrs = 100; % num simulated trials per hdg/condition
[tuning_curves, posterior, pop_lh, xq, simdata] = ...
    dots3DMP_likelihood_decoding(data,meanFRs,hdgs,mods,cohs,deltas,numtrs,step);

%% Fig 4 in paper...

figure('position',[300 300 600 400],'color','w'); hold on

% select a few simulated trials to plot
simTrs2plot = randperm(numtrs,5);

% vestibular only
ax=subplot(231); hold on;
m=1; c=1; d=2; h=5;
% boolean index for condition in simulated trials data struct                        
K = simdata.modality==m & simdata.coherence==cohs(c) & simdata.heading==hdgs(h) & simdata.delta==deltas(d);
temp = posterior(K,:);
plot(xq,temp(simTrs2plot,:),'linew',1,'color','k')
ax.XLim = [-6 6]; ax.TickDir = 'out'; ax.XTick = -5:2.5:5;
ylab = ylabel('Likelihood (p(r|\theta))');
ylab.Position(2) = -0.1;
ylim([0 0.4]);

% visual only
ax=subplot(232); hold on;         
m=2; c=1; d=2; h=5;
K = simdata.modality==m & simdata.coherence==cohs(c) & simdata.heading==hdgs(h) & simdata.delta==deltas(d);
temp = posterior(K,:);
plot(xq,temp(simTrs2plot,:),'linew',1,'color','r')
m=2; c=2; d=2; h=5;
K = simdata.modality==m & simdata.coherence==cohs(c) & simdata.heading==hdgs(h) & simdata.delta==deltas(d);
temp = posterior(K,:);
plot(xq,temp(simTrs2plot,:),'linew',1,'color','m','linestyle',':')
ax.XLim = [-6 6]; ax.TickDir = 'out'; ax.XTick = -5:2.5:5;
ylim([0 0.4]);

% combined, delta = +4, low coh
% visual is to the right of vestibular, at low coh, likelihood decoding is
% biased towards vestibular (left)
ax=subplot(234); hold on;                   
m=3; c=1; d=3; h=5;
K = simdata.modality==m & simdata.coherence==cohs(c) & simdata.heading==hdgs(h) & simdata.delta==deltas(d);
temp = posterior(K,:);
plot(xq,temp(simTrs2plot,:),'linew',1,'color','c','linestyle','-')
ax.XLim = [-6 6]; ax.TickDir = 'out'; ax.XTick = -5:2.5:5;
ylim([0 0.4]);

% combined, delta = +4, high coh
% visual is to the right of vestibular, at high coh, likelihood decoding is
% biased towards visual (right)
ax=subplot(235); hold on;         
m=3; c=2; d=3; h=5;
K = simdata.modality==m & simdata.coherence==cohs(c) & simdata.heading==hdgs(h) & simdata.delta==deltas(d);
temp = posterior(K,:);
plot(xq,temp(simTrs2plot,:),'linew',1,'color','b','linestyle','-')
ax.XLim = [-6 6]; ax.TickDir = 'out'; ax.XTick = -5:2.5:5;
ylim([0 0.4]);

% #TODO - cumulative Gaussian fits to simChoices for psychometric curves
% (4e + f)
ax = subplot(233);
ax = subplot(236);

%% Fig. 5 congruency index and tuning curves
% correlation between firing rate and heading for ves and vis separately,
% then product of these two

% use non-interpolated mean firing rates
vesFRs = squeeze(meanFRs(mods==1,1,2,:,:));
visFRs = permute(squeeze(meanFRs(mods==2,:,2,:,:)),[2 3 1]);
corrhdgs = hdgs;

% use interpolated tuning curves instead? no, probably artificially inflates
% correlations
% vesFRs = squeeze(tuning_curves(mods==1,1,2,:,:));
% visFRs = permute(squeeze(tuning_curves(mods==2,:,2,:,:)),[2 3 1]);
% corrhdgs = xq';

clear vesCorr visCorr* vesP visP
for u=1:size(vesFRs,2)
    [vesCorr(u,1),vesP(u,1)] = corr(vesFRs(:,u),corrhdgs);
    for c=1:length(cohs)
        [visCorr(u,c),visP(u,c)] = corr(visFRs(:,u,c),corrhdgs);
    end
end

[congruencyIndex,isSignificant] = dots3DMP_plot_congruency_NN(...
    monkUnit,tuning_curves,vesCorr,visCorr,vesP,visP,xq);

%% 




