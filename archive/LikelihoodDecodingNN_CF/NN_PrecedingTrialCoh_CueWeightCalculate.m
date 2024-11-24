%% Label indexing variables
clear;
% addpath(genpath('C:\Users\yhaile2\Documents\CODE\GitHubCodes\Fetschlab\FLprojects\dots3DMP\Fetsch2011_NNcodes\'))
% cd('C:\Users\yhaile2\Documents\CODE\MATLAB\3DMP\SensoryExpectation\NN2011\Fetsch2011NatNeuro_Data\')
load Fetsch_et_al_2011_behavOnly  % ~17500 trials of each nonzero delta
% cd('C:\Users\yhaile2\Documents\CODE\MATLAB\3DMP\3DMP_Lucio\BehaviorDataFiles_L\CleanDataFiles_L\')
% load lucio_20220512-20230605_clean

data.delta(data.delta==4) = 3;
data.delta(data.delta==-4) = -3;
mod = data.modality;
coh = data.coherence;
delta = data.delta;
mods = unique(data.modality);
hdngs = unique(data.heading);
cohs = unique(data.coherence);
deltas = unique(data.delta);
% conf = data.PDW; Uncomment for Lucio Data

% builds index to trials preceded by specified coh (< or > .5) PTC x1
% for tr = 2: length(coh)
%     if coh(tr-1) == cohs(1)
%         precLoCoh(tr) = true;
%     else
%         precLoCoh(tr) = false;
%     end
%     if coh(tr-1) == cohs(2)
%         precHiCoh(tr) = true;
%     else
%         precHiCoh(tr) = false;
%     end
% end

 %% PTCx1 excluding PTCx2 trials (Result: No consistent PTC effect observed)
% for tr = 3: length(coh)
%     if coh(tr-1) == cohs(1) && coh(tr-2) ~= cohs(1) && mod~=1  % if trial-1 and trial-2 coh are both low and not Vestib only trial
%         precLoCoh(tr) = true;
%     else
%         precLoCoh(tr) = false;
%     end
%     if coh(tr-1) == cohs(2) && coh(tr-2) ~= cohs(2) % PTC Hi x2, visual or comb modality prec trial, any delta
%         precHiCoh(tr) = true;
%     else
%         precHiCoh(tr) = false;
%     end
% end
% 
% % x1 preceding trial coh influences cue weights on current trial. Prec coh as Low vs. High vs. Either (control/baseline)
% trPerDelta = 10000; % Heading and coherence of trial are sampled randomly sampled in new boostrapped dataset, relying on their even distribution in the data being resampled
% tic
% for b = 1:10 % Number of boostrap resamplings
%     % Low coh preceding
%     d1_locohPrec = randsample(find(precLoCoh' & mod==3 & data.delta==deltas(1)),trPerDelta,1);
%     d2_locohPrec = randsample(find(precLoCoh' & mod==3 & data.delta==deltas(2)),trPerDelta,1);
%     d3_locohPrec = randsample(find(precLoCoh' & mod==3 & data.delta==deltas(3)),trPerDelta,1);
%     trt_locohPrec = [d1_locohPrec; d2_locohPrec; d3_locohPrec]; % list of the trials selected for this sample, to be used for indexing behaviora data
%     % High coh preceding
%     d1_hicohPrec = randsample(find(precHiCoh' & mod==3 & data.delta==deltas(1)),trPerDelta,1);
%     d2_hicohPrec = randsample(find(precHiCoh' & mod==3 & data.delta==deltas(2)),trPerDelta,1);
%     d3_hicohPrec = randsample(find(precHiCoh' & mod==3 & data.delta==deltas(3)),trPerDelta,1);
%     trt_hicohPrec = [d1_hicohPrec; d2_hicohPrec; d3_hicohPrec]; % list of the trials selected for this sample, to be used for indexing behaviora data
%     % Any preceding coh (control comparison)
%     precCoh = precLoCoh | precHiCoh;
%     d1_anyPrec = randsample(find(precCoh' & mod==3 & data.delta==deltas(1)),trPerDelta,1);
%     d2_anyPrec = randsample(find(precCoh' & mod==3 & data.delta==deltas(2)),trPerDelta,1);
%     d3_anyPrec = randsample(find(precCoh' & mod==3 & data.delta==deltas(3)),trPerDelta,1);
%     trt_AnyPrec = [d1_anyPrec; d2_anyPrec; d3_anyPrec]; % list of the trials selected for this sample, to be used for indexing behaviora data
% 
%     fn = fieldnames(data);
%     for f = 1:length(fn)
%         locohPrecData.(fn{f}) = data.(fn{f})(trt_locohPrec);
%         hicohPrecData.(fn{f}) = data.(fn{f})(trt_hicohPrec);
%         anyPrecData.(fn{f}) = data.(fn{f})(trt_AnyPrec);
%     end
% 
%     % weights following Low coherence trials
%     gfit_locohPrec = dots3DMP_fit_cgauss_NN(locohPrecData,mods,cohs,deltas);
%     wVes_loCohP(b,:) = dots3DMP_wgts_thres_NN(gfit_locohPrec.muPMF,gfit_locohPrec.sigmaPMF,cohs,deltas);
%     % cue weight following High coherence trials
%     gfit_hicohPrec = dots3DMP_fit_cgauss_NN(hicohPrecData,mods,cohs,deltas);
%     wVes_hiCohP(b,:) = dots3DMP_wgts_thres_NN(gfit_hicohPrec.muPMF,gfit_hicohPrec.sigmaPMF,cohs,deltas);
%     % weights for trials in general
%     gfit_anyPrec = dots3DMP_fit_cgauss_NN(anyPrecData,mods,cohs,deltas);
%     wVes_anyCohP(b,:) = dots3DMP_wgts_thres_NN(gfit_anyPrec.muPMF,gfit_anyPrec.sigmaPMF,cohs,deltas);
% end
% wVes_loCohPrec_boostrapAvg = mean(wVes_loCohP,1);
% wVes_hiCohPrec_boostrapAvg = mean(wVes_hiCohP,1);
% wVes_anyCohPrec_boostrapAvg = mean(wVes_anyCohP,1);
% 
% % Does vestibular cue weight increase following Low coh trial x1?
% wVes_loPvHiP = wVes_loCohPrec_boostrapAvg > wVes_hiCohPrec_boostrapAvg
% wVes_loPvAnyP = wVes_loCohPrec_boostrapAvg > wVes_anyCohPrec_boostrapAvg
% wVes_HiPvAnyP = wVes_hiCohPrec_boostrapAvg < wVes_anyCohPrec_boostrapAvg
% % Statistical significance
% loPvsHiP_test = ttest2(wVes_loCohP,wVes_hiCohP)
% loPvsAnyP_ttest = ttest2(wVes_loCohP,wVes_anyCohP)
% hiPvsAnyP_ttest = ttest2(wVes_hiCohP,wVes_anyCohP)
% 
% % x1 vs x2 Preceding comparison
% 
% % ttest2(wVes_loCohP,wVes_loCohPX2)
% % mean(wVes_loCohP,1) < mean(wVes_loCohPX2,1)
% % ttest2(wVes_hiCohP,wVes_hiCohPX2)
% % mean(wVes_hiCohP,1) > mean(wVes_hiCohPX2,1)
% % ttest2(wVes_anyCohP,wVes_anyCohPX2)
% 
% % Plotting vesWeights
% wVes_loCohPrec_boostrapAvg = mean(wVes_loCohP,1); % vesWeight means (recalc each time since shared var with x2)
% wVes_hiCohPrec_boostrapAvg = mean(wVes_hiCohP,1);
% wVes_anyCohPrec_boostrapAvg = mean(wVes_anyCohP,1);
% % x-axes
% xLoPrecHi = ones(size(wVes_loCohP(:,1)))*1;
% xHiPrecHi = ones(size(wVes_loCohP(:,1)))*1.2;
% xAnyPrecHi = ones(size(wVes_loCohP(:,1)))*1.1;
% xLoPrecLo = ones(size(wVes_loCohP(:,1)))*2;
% xHiPrecLo = ones(size(wVes_loCohP(:,1)))*2.2;
% xAnyPrecLo = ones(size(wVes_loCohP(:,1)))*2.1;
% % Consolidate plotting variables
% xHiCohTr = [xLoPrecHi xAnyPrecHi xHiPrecHi];
% xLoCohTr = [xLoPrecLo xAnyPrecLo xHiPrecLo];
% HighCohTrialsW = [wVes_loCohP(:,2) wVes_anyCohP(:,2) wVes_hiCohP(:,2)];
% LowCohTrialsW = [wVes_loCohP(:,1) wVes_anyCohP(:,1) wVes_hiCohP(:,1)];
% wVes_boostrpAvgs_LoCohs = [wVes_loCohPrec_boostrapAvg(:,1) wVes_anyCohPrec_boostrapAvg(:,1) wVes_hiCohPrec_boostrapAvg(:,1)];
% wVes_boostrpAvgs_HiCohs = [wVes_loCohPrec_boostrapAvg(:,2) wVes_anyCohPrec_boostrapAvg(:,2) wVes_hiCohPrec_boostrapAvg(:,2)];
% 
% clr = {'b','k','r'};
% hold on;
% for pl = 1:length(clr)
%     plot(xHiCohTr(:,pl),HighCohTrialsW(:,pl),'.','Color',clr{pl},'MarkerSize',8)
%     plot(xLoCohTr(:,pl),LowCohTrialsW(:,pl),'.','Color',clr{pl})
%     plot(xHiCohTr(1,pl),wVes_boostrpAvgs_HiCohs(pl),'o','Color',clr{pl})
%     plot(xLoCohTr(1,pl),wVes_boostrpAvgs_LoCohs(pl),'o','Color',clr{pl})
% end
% hold off;
% 
% ax = gca;
% ax.XTick = [1.1 2.1];
% ax.XTickLabels = {'High coh trials','Low coh trials'};
% title('Previous trial coherence affects cue weighting'); ylabel('Vestibular cue weight');
% ax.XLim = [.75,2.45];
% 
% toc % 100 loops = 11mins run time
% 
%% x2 Preceding trials of same Coh, maintain or magnify effect of 1 prec trial on cue weights?
% Does preceding trial coh influence cue weights on current trial? Prec coh as Low vs. High vs. Either (control/baseline)

% builds index to trials preceded by specified trial type (modality, coherence etc.)

% PTC Low x2 Mod 2|3 vs PTC Hi x2, Mod 2|3
% for tr = 3:length(coh) 
%     if coh(tr-1) == cohs(1) && coh(tr-2) == cohs(1) && mod(tr-1)~=1 && mod(tr-2)~=1 % Neither of 2 consec prec trials can be vestib
%         precLoCoh2x(tr) = true;
%     else
%         precLoCoh2x(tr) = false;
%     end
%     if coh(tr-1) == cohs(2) && coh(tr-2) == cohs(2)
%         precHiCoh2x(tr) = true;
%     else
%         precHiCoh2x(tr) = false;
%     end
% end

% PT Vestib Only x2 
% for tr = 3:length(coh) 
%     if mod(tr-1)==1 && mod(tr-2)==1 % 2 consec vestib only trials
%         precLoCoh2x(tr) = true; % Blue on plot = PT vestib only x2
%     else
%         precLoCoh2x(tr) = false;
%     end
%     if coh(tr-1) == cohs(1) && coh(tr-2) == cohs(1) && mod(tr-1)~=1 && mod(tr-2)~=1 % PTC Low x2
%         precHiCoh2x(tr) = true; %  Red on plot
%     else
%         precHiCoh2x(tr) = false;
%     end
% end

% PTC Low Mod2 vs PTC Low Mod3  || PTC Hi Mod2 vs PTC Hi Mod3
% for tr = 3:length(coh) 
%     if coh(tr-1) == cohs(1) && coh(tr-2) == cohs(1) && mod(tr-1)~=1 && mod(tr-2)~=1 % PTC Lo X2 mod 2 or 3
%         precLoCoh2x(tr) = true; % Blue
%     else
%         precLoCoh2x(tr) = false;
%     end
%     if coh(tr-1) == cohs(1) && coh(tr-2) == cohs(1) && mod(tr-1)==2 && mod(tr-2)==2 % PTC Low, mod 2 only
%         precHiCoh2x(tr) = true; % Red
%     else
%         precHiCoh2x(tr) = false;
%     end
% end

%% Mod2 vs Mod3 PTC effects seperated
for tr = 3:length(coh) 
    if coh(tr-1) == cohs(1) && coh(tr-2) == cohs(1) && mod(tr-1)==3 && mod(tr-2)==3 % PTC Lo Coh X2 mod 3 only
        PTCLoX2_mod3(tr) = true; % Light Green
    else
        PTCLoX2_mod3(tr) = false;
    end
    if coh(tr-1) == cohs(1) && coh(tr-2) == cohs(1) && mod(tr-1)==2 && mod(tr-2)==2 % PTC Low, mod 2 only
        PTCLoX2_mod2(tr) = true; % Light Red
    else
        PTCLoX2_mod2(tr) = false;
    end
    if coh(tr-1) == cohs(2) && coh(tr-2) == cohs(2) && mod(tr-1)==3 && mod(tr-2)==3 % PTC Hi X2 mod 3 only
        PTCHiX2_mod3(tr) = true; % Dark Green
    else
        PTCHiX2_mod3(tr) = false;
    end
    if coh(tr-1) == cohs(2) && coh(tr-2) == cohs(2) && mod(tr-1)==2 && mod(tr-2)==2 % PTC Hi x2, mod 2 only
        PTCHiX2_mod2(tr) = true; % Dark Red
    else
        PTCHiX2_mod2(tr) = false;
    end
    if mod(tr-1)~=mod(tr-2) && coh(tr-1)~=coh(tr-2)
        precX2_cntrl(tr) = true;
    else
        precX2_cntrl(tr) = false;
    end
end

trPerDelta = 1500; % Heading and coherence of trial are sampled randomly sampled in each boostrapped dataset, relying on their even distribution in the data being sampled, while delta is specified due to its uneven sampling in the OG dataset
tic
for b = 1:50 % number of trial resamplings ('bootstraps') 
    % Lo Coh X2 mod 2 only
    d1_locohPrec = randsample(find(PTCLoX2_mod2' & mod==3 & data.delta==deltas(1)),trPerDelta,1);
    d2_locohPrec = randsample(find(PTCLoX2_mod2' & mod==3 & data.delta==deltas(2)),trPerDelta,1);
    d3_locohPrec = randsample(find(PTCLoX2_mod2' & mod==3 & data.delta==deltas(3)),trPerDelta,1);
    trt_PTC_LoCohX2_mod2 = [d1_locohPrec; d2_locohPrec; d3_locohPrec]; % list of the trials selected for this sample, to be used for indexing behaviora data
    % Lo Coh X2 mod 3 only
    d1_hicohPrec = randsample(find(PTCLoX2_mod3' & mod==3 & data.delta==deltas(1)),trPerDelta,1);
    d2_hicohPrec = randsample(find(PTCLoX2_mod3' & mod==3 & data.delta==deltas(2)),trPerDelta,1);
    d3_hicohPrec = randsample(find(PTCLoX2_mod3' & mod==3 & data.delta==deltas(3)),trPerDelta,1);
    trt_PTC_LoCohX2_mod3 = [d1_hicohPrec; d2_hicohPrec; d3_hicohPrec]; % list of the trials selected for this sample, to be used for indexing behaviora data
    % Hi Coh X2 mod 2 only
    d1_locohPrec = randsample(find(PTCHiX2_mod2' & mod==3 & data.delta==deltas(1)),trPerDelta,1);
    d2_locohPrec = randsample(find(PTCHiX2_mod2' & mod==3 & data.delta==deltas(2)),trPerDelta,1);
    d3_locohPrec = randsample(find(PTCHiX2_mod2' & mod==3 & data.delta==deltas(3)),trPerDelta,1);
    trt_PTC_HiCohX2_mod2 = [d1_locohPrec; d2_locohPrec; d3_locohPrec]; % list of the trials selected for this sample, to be used for indexing behaviora data
    % Hi Coh X2 mod 3 only
    d1_hicohPrec = randsample(find(PTCHiX2_mod3' & mod==3 & data.delta==deltas(1)),trPerDelta,1);
    d2_hicohPrec = randsample(find(PTCHiX2_mod3' & mod==3 & data.delta==deltas(2)),trPerDelta,1);
    d3_hicohPrec = randsample(find(PTCHiX2_mod3' & mod==3 & data.delta==deltas(3)),trPerDelta,1);
    trt_PTC_HiCohX2_mod3 = [d1_hicohPrec; d2_hicohPrec; d3_hicohPrec]; % list of the trials selected for this sample, to be used for indexing behaviora data
    
    % (control comparison)
%     precX2_cntrl = PTCLoX2_mod2 | PTCLoX2_mod3 | PTCHiX2_mod2 | PTCHiX2_mod3;
    d1_anyPrec = randsample(find(precX2_cntrl' & mod==3 & data.delta==deltas(1)),trPerDelta,1);
    d2_anyPrec = randsample(find(precX2_cntrl' & mod==3 & data.delta==deltas(2)),trPerDelta,1);
    d3_anyPrec = randsample(find(precX2_cntrl' & mod==3 & data.delta==deltas(3)),trPerDelta,1);
    trt_AnyPrecX2 = [d1_anyPrec; d2_anyPrec; d3_anyPrec]; % list of the trials selected for this sample, to be used for indexing behaviora data

    fn = fieldnames(data);
    for f = 1:length(fn)
        PTC_LoCohDataX2_mod2.(fn{f}) = data.(fn{f})(trt_PTC_LoCohX2_mod2);
        PTC_LoCohDataX2_mod3.(fn{f}) = data.(fn{f})(trt_PTC_LoCohX2_mod3);
        PTC_HiCohDataX2_mod2.(fn{f}) = data.(fn{f})(trt_PTC_HiCohX2_mod2);
        PTC_HiCohDataX2_mod3.(fn{f}) = data.(fn{f})(trt_PTC_HiCohX2_mod3);
        anyPrecDataX2.(fn{f}) = data.(fn{f})(trt_AnyPrecX2);
    end

    % Calculate and store weight for each boostrap iteration of trials for each PTC condition
    % weights following Low coherence trials
    gfit_locohPrec_mod2 = dots3DMP_fit_cgauss_NN(PTC_LoCohDataX2_mod2,mods,cohs,deltas);
    wVes_loCohPX2_mod2(b,:) = dots3DMP_wgts_thres_NN(gfit_locohPrec_mod2.muPMF,gfit_locohPrec_mod2.sigmaPMF,cohs,deltas);
    gfit_locohPrec_mod3 = dots3DMP_fit_cgauss_NN(PTC_LoCohDataX2_mod3,mods,cohs,deltas);
    wVes_loCohPX2_mod3(b,:) = dots3DMP_wgts_thres_NN(gfit_locohPrec_mod3.muPMF,gfit_locohPrec_mod3.sigmaPMF,cohs,deltas);
    % cue weight following High coherence trials
    gfit_HicohPrec_mod2 = dots3DMP_fit_cgauss_NN(PTC_HiCohDataX2_mod2,mods,cohs,deltas);
    wVes_HiCohPX2_mod2(b,:) = dots3DMP_wgts_thres_NN(gfit_HicohPrec_mod2.muPMF,gfit_HicohPrec_mod2.sigmaPMF,cohs,deltas);
    gfit_HicohPrec_mod3 = dots3DMP_fit_cgauss_NN(PTC_HiCohDataX2_mod3,mods,cohs,deltas);
    wVes_HiCohPX2_mod3(b,:) = dots3DMP_wgts_thres_NN(gfit_HicohPrec_mod3.muPMF,gfit_HicohPrec_mod3.sigmaPMF,cohs,deltas);
    % weights for control trials ---- in general following 2 consecutive same coh trials anyX2
    gfit_anyPrec = dots3DMP_fit_cgauss_NN(anyPrecDataX2,mods,cohs,deltas);
    wVes_anyCohPX2(b,:) = dots3DMP_wgts_thres_NN(gfit_anyPrec.muPMF,gfit_anyPrec.sigmaPMF,cohs,deltas);
end

% Avg cue weight across all boostraps, for trials of each PTC condition
% PTC Low
wVes_loCohPrec_boostrapAvg_mod2 = mean(wVes_loCohPX2_mod2,1);
wVes_loCohPrec_boostrapAvg_mod3 = mean(wVes_loCohPX2_mod3,1);
% PTC High
wVes_HiCohPrec_boostrapAvg_mod2 = mean(wVes_HiCohPX2_mod2,1);
wVes_HiCohPrec_boostrapAvg_mod3 = mean(wVes_HiCohPX2_mod3,1);
% Control PTC
wVes_anyCohPrec_boostrapAvg = mean(wVes_anyCohPX2,1);

%     % PTC Low stats
% %vestibular weight comparisons 
% wVes_loPmod2vloPmod3 = wVes_loCohPrec_boostrapAvg_mod2 > wVes_loCohPrec_boostrapAvg_mod3 % mod 2 vis uncertainty = more elevated vesW than mod3 vis uncertainty?
% wVes_loPmod2vAnyPX2 = wVes_loCohPrec_boostrapAvg_mod2 > wVes_anyCohPrec_boostrapAvg % Does vestibular cue weight increase following Low coh trial?
% wVes_loPmod3vAnyPX2 = wVes_loCohPrec_boostrapAvg_mod3 > wVes_anyCohPrec_boostrapAvg % Does vestibular cue weight increase following Low coh trial?
% % Statistical significance
% loPmod2vsloPmod3_testX2 = ttest2(wVes_loCohPX2_mod2,wVes_loCohPX2_mod3)
% loPmod2vsAnyP_ttestX2 = ttest2(wVes_loCohPX2_mod2,wVes_anyCohPX2)
% loPmod3vsAnyP_ttestX2 = ttest2(wVes_loCohPX2_mod3,wVes_anyCohPX2)
%     % PTC High stats
% %vestibular weight comparisons 
% wVes_HiPmod2vHiPmod3 = wVes_HiCohPrec_boostrapAvg_mod2 > wVes_HiCohPrec_boostrapAvg_mod3 % mod 2 vis uncertainty = more elevated vesW than mod3 vis uncertainty?
% wVes_HiPmod2vAnyPX2 = wVes_HiCohPrec_boostrapAvg_mod2 > wVes_anyCohPrec_boostrapAvg % Does vestibular cue weight increase folHiwing Hiw coh trial?
% wVes_HiPmod3vAnyPX2 = wVes_HiCohPrec_boostrapAvg_mod3 > wVes_anyCohPrec_boostrapAvg % Does vestibular cue weight increase folHiwing Hiw coh trial?
% % Statistical significance
% HiPmod2vsHiPmod3_testX2 = ttest2(wVes_HiCohPX2_mod2,wVes_HiCohPX2_mod3)
% HiPmod2vsAnyP_ttestX2 = ttest2(wVes_HiCohPX2_mod2,wVes_anyCohPX2)
% HiPmod3vsAnyP_ttestX2 = ttest2(wVes_HiCohPX2_mod3,wVes_anyCohPX2)

ttest2(wVes_HiCohPX2_mod3(:,1),wVes_loCohPX2_mod3(:,1)) % Hi vs Low PTC result in significantly different vesW on Low coh trials

% Plotting vesWeights
% x-axis
xLoPrecHi_mod2 = ones(size(wVes_loCohPX2_mod2(:,1)))*1;
xLoPrecHi_mod3 = ones(size(wVes_loCohPX2_mod2(:,1)))*1.2;
xHiPrecHi_mod2 = ones(size(wVes_loCohPX2_mod2(:,1)))*1;
xHiPrecHi_mod3 = ones(size(wVes_loCohPX2_mod2(:,1)))*1.2;
xAnyPrecHi = ones(size(wVes_loCohPX2_mod2(:,1)))*1.1;
xLoPrecLo_mod2 = ones(size(wVes_loCohPX2_mod2(:,1)))*2;
xLoPrecLo_mod3 = ones(size(wVes_loCohPX2_mod2(:,1)))*2.2;
xHiPrecLo_mod2 = ones(size(wVes_loCohPX2_mod2(:,1)))*2;
xHiPrecLo_mod3 = ones(size(wVes_loCohPX2_mod2(:,1)))*2.2;
xAnyPrecLo = ones(size(wVes_loCohPX2_mod2(:,1)))*2.1;
% Consolidate plotting variables
xHiCohTr = [xLoPrecHi_mod2 xLoPrecHi_mod3 xHiPrecHi_mod2 xHiPrecHi_mod3 xAnyPrecHi];
xLoCohTr = [xLoPrecLo_mod2 xLoPrecLo_mod3 xHiPrecLo_mod2 xHiPrecLo_mod3 xAnyPrecLo];

HighCohTrialsW = [wVes_loCohPX2_mod2(:,2) wVes_loCohPX2_mod3(:,2) wVes_HiCohPX2_mod2(:,2) wVes_HiCohPX2_mod3(:,2) wVes_anyCohPX2(:,2)];
LowCohTrialsW = [wVes_loCohPX2_mod2(:,1) wVes_loCohPX2_mod3(:,1) wVes_HiCohPX2_mod2(:,1) wVes_HiCohPX2_mod3(:,1) wVes_anyCohPX2(:,1)];

wVes_boostrpAvgs_LoCohs = [ wVes_loCohPrec_boostrapAvg_mod2(:,1) wVes_loCohPrec_boostrapAvg_mod3(:,1) wVes_HiCohPrec_boostrapAvg_mod2(:,1) wVes_HiCohPrec_boostrapAvg_mod3(:,1) wVes_anyCohPrec_boostrapAvg(:,1) ];
wVes_boostrpAvgs_HiCohs = [ wVes_loCohPrec_boostrapAvg_mod2(:,2) wVes_loCohPrec_boostrapAvg_mod3(:,2) wVes_HiCohPrec_boostrapAvg_mod2(:,2) wVes_HiCohPrec_boostrapAvg_mod3(:,2) wVes_anyCohPrec_boostrapAvg(:,2) ];

clr = {[1.00,0.40,0.40], [0.46,0.92,0.39], 'r',[0.15,0.61,0.01],'k'};
% clr = {'b','k','r'}; % LRed, LGreen, DRed, DGreen, Blk
hold on;
for pl = 1:length(clr)
    plot(xHiCohTr(:,pl),HighCohTrialsW(:,pl),'.','Color',clr{pl},'MarkerSize',8)
    plot(xLoCohTr(:,pl),LowCohTrialsW(:,pl),'.','Color',clr{pl})
    plot(xHiCohTr(1,pl),wVes_boostrpAvgs_HiCohs(pl),'o','Color',clr{pl})
    plot(xLoCohTr(1,pl),wVes_boostrpAvgs_LoCohs(pl),'o','Color',clr{pl})
end
hold off;

ax = gca;
ax.XTick = [1.1 2.1];
ax.XTickLabels = {'High coh trials','Low coh trials'};
title('Previous trial coherence and modality affect cue weighting x2'); ylabel('Vestibular cue weight');
ax.XLim = [.75,2.45];

toc % 100 loops = 11mins run time 


%% Vesibular only vs. Visual only PT x2 effect on Cue weights
% Clear Vars !!!
for tr = 3:length(coh) 
   
    if coh(tr-1) == cohs(1) && coh(tr-2) == cohs(1) && mod(tr-1)==2 && mod(tr-2)==2 % PTC Low, mod 2 only
        PTCLoX2_mod2(tr) = true; % Light Red
    else
        PTCLoX2_mod2(tr) = false;
    end
    if coh(tr-1) == cohs(2) && coh(tr-2) == cohs(2) && mod(tr-1)==2 && mod(tr-2)==2 % PTC Hi x2, mod 2 only
        PTCHiX2_mod2(tr) = true; % Dark Red
    else
        PTCHiX2_mod2(tr) = false;
    end
    if mod(tr-1) == mods(1) && mod(tr-2) == mods(1)
        PTvesOnly_mod1X2(tr) = true; % Blue
    else
        PTvesOnly_mod1X2(tr) = false;
    end
    if mod(tr-1)~=mod(tr-2) && coh(tr-1)~=coh(tr-2)
        precX2_cntrl(tr) = true;
    else
        precX2_cntrl(tr) = false;
    end

end

trPerDelta = 1500; % Heading and coherence of trial are sampled randomly sampled in each boostrapped dataset, relying on their even distribution in the data being sampled, while delta is specified due to its uneven sampling in the OG dataset
tic
for b = 1:50 % number of trial resamplings ('bootstraps') 
    % Lo Coh X2 mod 2 only
    d1_locohPrec = randsample(find(PTCLoX2_mod2' & mod==3 & data.delta==deltas(1)),trPerDelta,1);
    d2_locohPrec = randsample(find(PTCLoX2_mod2' & mod==3 & data.delta==deltas(2)),trPerDelta,1);
    d3_locohPrec = randsample(find(PTCLoX2_mod2' & mod==3 & data.delta==deltas(3)),trPerDelta,1);
    trt_PTC_LoCohX2_mod2 = [d1_locohPrec; d2_locohPrec; d3_locohPrec]; % list of the trials selected for this sample, to be used for indexing behaviora data
    % Hi Coh X2 mod 2 only
    d1_HicohPrec = randsample(find(PTCHiX2_mod2' & mod==3 & data.delta==deltas(1)),trPerDelta,1);
    d2_HicohPrec = randsample(find(PTCHiX2_mod2' & mod==3 & data.delta==deltas(2)),trPerDelta,1);
    d3_HicohPrec = randsample(find(PTCHiX2_mod2' & mod==3 & data.delta==deltas(3)),trPerDelta,1);
    trt_PTC_HiCohX2_mod2 = [d1_HicohPrec; d2_HicohPrec; d3_HicohPrec];
    % Ves only x2 mod 1 only
    d1_VesOnlyPrec = randsample(find(PTvesOnly_mod1X2' & mod==3 & data.delta==deltas(1)),trPerDelta,1);
    d2_VesOnlyPrec = randsample(find(PTvesOnly_mod1X2' & mod==3 & data.delta==deltas(2)),trPerDelta,1);
    d3_VesOnlyPrec = randsample(find(PTvesOnly_mod1X2' & mod==3 & data.delta==deltas(3)),trPerDelta,1);
    trt_PT_VesOnlyX2_mod1 = [d1_VesOnlyPrec; d2_VesOnlyPrec; d3_VesOnlyPrec];
    % (control comparison)
%     precX2_cntrl = PTCLoX2_mod2 | PTCHiX2_mod2 | PTvesOnly_mod1;
    d1_anyPrec = randsample(find(precX2_cntrl' & mod==3 & data.delta==deltas(1)),trPerDelta,1);
    d2_anyPrec = randsample(find(precX2_cntrl' & mod==3 & data.delta==deltas(2)),trPerDelta,1);
    d3_anyPrec = randsample(find(precX2_cntrl' & mod==3 & data.delta==deltas(3)),trPerDelta,1);
    trt_AnyPrecX2 = [d1_anyPrec; d2_anyPrec; d3_anyPrec]; % list of the trials selected for this sample, to be used for indexing behaviora data

    fn = fieldnames(data);
    for f = 1:length(fn)
        PTC_LoCohDataX2_mod2.(fn{f}) = data.(fn{f})(trt_PTC_LoCohX2_mod2);
        PTC_HiCohDataX2_mod2.(fn{f}) = data.(fn{f})(trt_PTC_HiCohX2_mod2);
        PTC_VesOnlyDataX2_mod1.(fn{f}) = data.(fn{f})(trt_PT_VesOnlyX2_mod1);
        anyPrecDataX2.(fn{f}) = data.(fn{f})(trt_AnyPrecX2);
    end

% Calculate and store weight for each boostrap iteration of trials for each PTC condition
    % weights following Mod 2 trials
    gfit_locohPrec_mod2 = dots3DMP_fit_cgauss_NN(PTC_LoCohDataX2_mod2,mods,cohs,deltas);
    wVes_loCohPX2_mod2(b,:) = dots3DMP_wgts_thres_NN(gfit_locohPrec_mod2.muPMF,gfit_locohPrec_mod2.sigmaPMF,cohs,deltas);
    gfit_HicohPrec_mod2 = dots3DMP_fit_cgauss_NN(PTC_HiCohDataX2_mod2,mods,cohs,deltas);
    wVes_HiCohPX2_mod2(b,:) = dots3DMP_wgts_thres_NN(gfit_HicohPrec_mod2.muPMF,gfit_HicohPrec_mod2.sigmaPMF,cohs,deltas);
    % weights following Mod 1 trials
    gfit_VesOnlyPrec_mod1X2 = dots3DMP_fit_cgauss_NN(PTC_VesOnlyDataX2_mod1,mods,cohs,deltas);
    wVes_VesOnlyPX2_mod1(b,:) = dots3DMP_wgts_thres_NN(gfit_VesOnlyPrec_mod1X2.muPMF,gfit_VesOnlyPrec_mod1X2.sigmaPMF,cohs,deltas);
    % weights following assortment of trials (control)
    gfit_anyPrec = dots3DMP_fit_cgauss_NN(anyPrecDataX2,mods,cohs,deltas);
    wVes_anyCohPX2(b,:) = dots3DMP_wgts_thres_NN(gfit_anyPrec.muPMF,gfit_anyPrec.sigmaPMF,cohs,deltas);
end

% Avg cue weight across all boostraps, for trials of each PTC condition
% Mod 2 x2
wVes_loCohPrec_boostrapAvg_mod2 = mean(wVes_loCohPX2_mod2,1);
wVes_HiCohPrec_boostrapAvg_mod2= mean(wVes_HiCohPX2_mod2,1);
% Mod 1 x2
wVes_VesOnlyP_boostrapAvg_mod1= mean(wVes_VesOnlyPX2_mod1,1);
% Control PT
wVes_anyCohPrec_boostrapAvg = mean(wVes_anyCohPX2,1);

% Plotting vesWeights
% x-axis
% Prec Hi
xLoPrecHi_mod2 = ones(size(wVes_loCohPX2_mod2(:,1)))*1;
xHiPrecHi_mod2 = ones(size(wVes_loCohPX2_mod2(:,1)))*1;
xVesOnlyPrecHi_mod1X2 = ones(size(wVes_loCohPX2_mod2(:,1)))*1.2;
xAnyPrecHi = ones(size(wVes_loCohPX2_mod2(:,1)))*1.1;
% Prec Lo
xLoPrecLo_mod2 = ones(size(wVes_loCohPX2_mod2(:,1)))*2;
xHiPrecLo_mod2 = ones(size(wVes_loCohPX2_mod2(:,1)))*2;
xVesOnlyPrecLo_mod1X2 = ones(size(wVes_loCohPX2_mod2(:,1)))*2.2;
xAnyPrecLo = ones(size(wVes_loCohPX2_mod2(:,1)))*2.1;
% Consolidate plotting variables
xHiCohTr = [xLoPrecHi_mod2  xHiPrecHi_mod2 xVesOnlyPrecHi_mod1X2 xAnyPrecHi];
xLoCohTr = [xLoPrecLo_mod2  xHiPrecLo_mod2 xVesOnlyPrecLo_mod1X2 xAnyPrecLo];

HighCohTrialsW = [wVes_loCohPX2_mod2(:,2) wVes_HiCohPX2_mod2(:,2) wVes_VesOnlyPX2_mod1(:,2) wVes_anyCohPX2(:,2)];
LowCohTrialsW = [wVes_loCohPX2_mod2(:,1) wVes_HiCohPX2_mod2(:,1) wVes_VesOnlyPX2_mod1(:,1) wVes_anyCohPX2(:,1)];

wVes_boostrpAvgs_LoCohs = [ wVes_loCohPrec_boostrapAvg_mod2(:,1) wVes_HiCohPrec_boostrapAvg_mod2(:,1) wVes_VesOnlyP_boostrapAvg_mod1(:,1) wVes_anyCohPrec_boostrapAvg(:,1) ];
wVes_boostrpAvgs_HiCohs = [ wVes_loCohPrec_boostrapAvg_mod2(:,2) wVes_HiCohPrec_boostrapAvg_mod2(:,2) wVes_VesOnlyP_boostrapAvg_mod1(:,2) wVes_anyCohPrec_boostrapAvg(:,2) ];

clr = {[1.00,0.40,0.40], 'r', 'b', 'k'}; % LRed, DRed, Blue, Blk
figure(); hold on;
for pl = 1:length(clr)
    plot(xHiCohTr(:,pl),HighCohTrialsW(:,pl),'.','Color',clr{pl},'MarkerSize',8)
    plot(xLoCohTr(:,pl),LowCohTrialsW(:,pl),'.','Color',clr{pl})
    plot(xHiCohTr(1,pl),wVes_boostrpAvgs_HiCohs(pl),'o','Color',clr{pl})
    plot(xLoCohTr(1,pl),wVes_boostrpAvgs_LoCohs(pl),'o','Color',clr{pl})
end
hold off;

ax = gca;
ax.XTick = [1.1 2.1];
ax.XTickLabels = {'High coh trials','Low coh trials'};
title('Previous trial coherence and modality affect cue weighting x2'); ylabel('Vestibular cue weight');
ax.XLim = [.75,2.45];

toc % 100 loops = 11mins run time 

%% Modality 1 - x2 vs x1 consec, effect magnify?
% Control is all trials not preceded by mod 1
for tr = 3:length(coh) 
    
    if mod(tr-1) == mods(1) && mod(tr-2) == mods(1)
        PTvesOnly_mod1X2(tr) = true; % Dark Blue
    else
        PTvesOnly_mod1X2(tr) = false;
    end
    if mod(tr-1) == mods(1) % x1 Mod1 - Light blue
        PTvesOnly_mod1X1(tr) = true;
    else
        PTvesOnly_mod1X1(tr) = false;
    end
    if mod(tr-1)~=mod(tr-2) && coh(tr-1)~=coh(tr-2)
        precX2_cntrlMod1(tr) = true;
    else
        precX2_cntrlMod1(tr) = false;
    end

end

trPerDelta = 1500; % Heading and coherence of trial are sampled randomly sampled in each boostrapped dataset, relying on their even distribution in the data being sampled, while delta is specified due to its uneven sampling in the OG dataset
tic
for b = 1:5%:50 % number of trial resamplings ('bootstraps') 

    % Ves only x2 mod 1 only
    d1_VesOnlyPrec = randsample(find(PTvesOnly_mod1X2' & mod==3 & data.delta==deltas(1)),trPerDelta,1);
    d2_VesOnlyPrec = randsample(find(PTvesOnly_mod1X2' & mod==3 & data.delta==deltas(2)),trPerDelta,1);
    d3_VesOnlyPrec = randsample(find(PTvesOnly_mod1X2' & mod==3 & data.delta==deltas(3)),trPerDelta,1);
    trt_PT_VesOnlyX2_mod1 = [d1_VesOnlyPrec; d2_VesOnlyPrec; d3_VesOnlyPrec];
    % Ves only x1 mod 1 only
    d1_VesOnlyPrec = randsample(find(PTvesOnly_mod1X1' & mod==3 & data.delta==deltas(1)),trPerDelta,1);
    d2_VesOnlyPrec = randsample(find(PTvesOnly_mod1X1' & mod==3 & data.delta==deltas(2)),trPerDelta,1);
    d3_VesOnlyPrec = randsample(find(PTvesOnly_mod1X1' & mod==3 & data.delta==deltas(3)),trPerDelta,1);
    trt_PT_VesOnlyX1_mod1 = [d1_VesOnlyPrec; d2_VesOnlyPrec; d3_VesOnlyPrec];
    % (control comparison)
    d1_anyPrec = randsample(find(precX2_cntrlMod1' & mod==3 & data.delta==deltas(1)),trPerDelta,1);
    d2_anyPrec = randsample(find(precX2_cntrlMod1' & mod==3 & data.delta==deltas(2)),trPerDelta,1);
    d3_anyPrec = randsample(find(precX2_cntrlMod1' & mod==3 & data.delta==deltas(3)),trPerDelta,1);
    trt_AnyPrecX2 = [d1_anyPrec; d2_anyPrec; d3_anyPrec]; % list of the trials selected for this sample, to be used for indexing behaviora data

    fn = fieldnames(data);
    for f = 1:length(fn)
        PTC_VesOnlyDataX1_mod1.(fn{f}) = data.(fn{f})(trt_PT_VesOnlyX1_mod1);
        PTC_VesOnlyDataX2_mod1.(fn{f}) = data.(fn{f})(trt_PT_VesOnlyX2_mod1);
        anyPrecDataX2.(fn{f}) = data.(fn{f})(trt_AnyPrecX2);
    end

% Calculate and store weight for each boostrap iteration of trials for each PTC condition
    % weights following Mod 1 x2 trials
    gfit_VesOnlyPrec_mod1X2 = dots3DMP_fit_cgauss_NN(PTC_VesOnlyDataX2_mod1,mods,cohs,deltas);
    wVes_VesOnlyPX2_mod1(b,:) = dots3DMP_wgts_thres_NN(gfit_VesOnlyPrec_mod1X2.muPMF,gfit_VesOnlyPrec_mod1X2.sigmaPMF,cohs,deltas);
    % x1
    gfit_VesOnlyPrec_mod1X1 = dots3DMP_fit_cgauss_NN(PTC_VesOnlyDataX1_mod1,mods,cohs,deltas);
    wVes_VesOnlyPX1_mod1(b,:) = dots3DMP_wgts_thres_NN(gfit_VesOnlyPrec_mod1X1.muPMF,gfit_VesOnlyPrec_mod1X1.sigmaPMF,cohs,deltas);
    % weights following assortment of trials (control)
    gfit_anyPrec = dots3DMP_fit_cgauss_NN(anyPrecDataX2,mods,cohs,deltas);
    wVes_anyCohPX2(b,:) = dots3DMP_wgts_thres_NN(gfit_anyPrec.muPMF,gfit_anyPrec.sigmaPMF,cohs,deltas);
end

% Avg cue weight across all boostraps, for trials of each PTC condition

% Mod 1 x2
wVes_VesOnlyP_boostrapAvg_mod1X2= mean(wVes_VesOnlyPX2_mod1,1);
% Mod 1 x1
wVes_VesOnlyP_boostrapAvg_mod1X1= mean(wVes_VesOnlyPX1_mod1,1);
% Control PT
wVes_anyCohPrec_boostrapAvg = mean(wVes_anyCohPX2,1);

% Plotting vesWeights
% x-axies
% Prec Hi
xVesOnlyPrecHi_mod1X2 = ones(size(wVes_VesOnlyPX2_mod1(:,1)))*1.2;
xVesOnlyPrecHi_mod1X1 = ones(size(wVes_VesOnlyPX2_mod1(:,1)))*1.2;
xAnyPrecHi = ones(size(wVes_VesOnlyPX2_mod1(:,1)))*1.1;
% Prec Lo
xVesOnlyPrecLo_mod1X2 = ones(size(wVes_VesOnlyPX2_mod1(:,1)))*2.2;
xVesOnlyPrecLo_mod1X1 = ones(size(wVes_VesOnlyPX2_mod1(:,1)))*2.2;
xAnyPrecLo = ones(size(wVes_VesOnlyPX2_mod1(:,1)))*2.1;
% Consolidate plotting variables
% X-axis position
xHiCohTr = [xVesOnlyPrecHi_mod1X1 xVesOnlyPrecHi_mod1X2 xAnyPrecHi];
xLoCohTr = [xVesOnlyPrecLo_mod1X1 xVesOnlyPrecLo_mod1X2 xAnyPrecLo];

HighCohTrialsW = [ wVes_VesOnlyPX1_mod1(:,2) wVes_VesOnlyPX2_mod1(:,2) wVes_anyCohPX2(:,2)];
LowCohTrialsW = [wVes_VesOnlyPX1_mod1(:,1) wVes_VesOnlyPX2_mod1(:,1) wVes_anyCohPX2(:,1)];

wVes_boostrpAvgs_LoCohs = [ wVes_VesOnlyP_boostrapAvg_mod1X1(:,1) wVes_VesOnlyP_boostrapAvg_mod1X2(:,1) wVes_anyCohPrec_boostrapAvg(:,1) ];
wVes_boostrpAvgs_HiCohs = [ wVes_VesOnlyP_boostrapAvg_mod1X1(:,2) wVes_VesOnlyP_boostrapAvg_mod1X2(:,2) wVes_anyCohPrec_boostrapAvg(:,2) ];

clr = {'c', 'b', 'k'}; % LBlue, DBlue, Blk
figure(); hold on;
for pl = 1:length(clr)
    plot(xHiCohTr(:,pl),HighCohTrialsW(:,pl),'.','Color',clr{pl},'MarkerSize',8)
    plot(xLoCohTr(:,pl),LowCohTrialsW(:,pl),'.','Color',clr{pl})
    plot(xHiCohTr(1,pl),wVes_boostrpAvgs_HiCohs(pl),'o','Color',clr{pl})
    plot(xLoCohTr(1,pl),wVes_boostrpAvgs_LoCohs(pl),'o','Color',clr{pl})
end
hold off;

ax = gca;
ax.XTick = [1.1 2.1];
ax.XTickLabels = {'High coh trials','Low coh trials'};
title('Previous trial coherence and modality affect cue weighting x2'); ylabel('Vestibular cue weight');
ax.XLim = [.75,2.45];

toc % 100 loops = 11mins run time 





%% Mod 1 x2, Mod2 Lo and High PTC, Mod3 Lo and High PTC

for tr = 3:length(coh) 
    if coh(tr-1) == cohs(1) && coh(tr-2) == cohs(1) && mod(tr-1)==3 && mod(tr-2)==3 % PTC Lo Coh X2 mod 3 only
        PTCLoX2_mod3(tr) = true; % Light Green
    else
        PTCLoX2_mod3(tr) = false;
    end
    if coh(tr-1) == cohs(1) && coh(tr-2) == cohs(1) && mod(tr-1)==2 && mod(tr-2)==2 % PTC Low, mod 2 only
        PTCLoX2_mod2(tr) = true; % Light Red
    else
        PTCLoX2_mod2(tr) = false;
    end
    if mod(tr-1) == mods(1) && mod(tr-2) == mods(1) % VesOnly x2
        PTvesOnly_mod1X2(tr) = true; % Dark Blue
    else
        PTvesOnly_mod1X2(tr) = false;
    end
    if coh(tr-1) == cohs(2) && coh(tr-2) == cohs(2) && mod(tr-1)==3 && mod(tr-2)==3 % PTC Hi X2 mod 3 only
        PTCHiX2_mod3(tr) = true; % Dark Green
    else
        PTCHiX2_mod3(tr) = false;
    end
    if coh(tr-1) == cohs(2) && coh(tr-2) == cohs(2) && mod(tr-1)==2 && mod(tr-2)==2 % PTC Hi x2, mod 2 only
        PTCHiX2_mod2(tr) = true; % Dark Red
    else
        PTCHiX2_mod2(tr) = false;
    end
    if mod(tr-1)~=mod(tr-2) && coh(tr-1)~=coh(tr-2)
        precX2_cntrl(tr) = true;
    else
        precX2_cntrl(tr) = false;
    end
end

trPerDelta = 500; % Heading and coherence of trial are sampled randomly sampled in each boostrapped dataset, relying on their even distribution in the data being sampled, while delta is specified due to its uneven sampling in the OG dataset
tic
for b = 1:50 % number of trial resamplings ('bootstraps') 
    % Lo Coh X2 mod 2 only
    d1_locohPrec = randsample(find(PTCLoX2_mod2' & mod==3 & data.delta==deltas(1)),trPerDelta,1);
    d2_locohPrec = randsample(find(PTCLoX2_mod2' & mod==3 & data.delta==deltas(2)),trPerDelta,1);
    d3_locohPrec = randsample(find(PTCLoX2_mod2' & mod==3 & data.delta==deltas(3)),trPerDelta,1);
    trt_PTC_LoCohX2_mod2 = [d1_locohPrec; d2_locohPrec; d3_locohPrec]; % list of the trials selected for this sample, to be used for indexing behaviora data
    % Lo Coh X2 mod 3 only
    d1_hicohPrec = randsample(find(PTCLoX2_mod3' & mod==3 & data.delta==deltas(1)),trPerDelta,1);
    d2_hicohPrec = randsample(find(PTCLoX2_mod3' & mod==3 & data.delta==deltas(2)),trPerDelta,1);
    d3_hicohPrec = randsample(find(PTCLoX2_mod3' & mod==3 & data.delta==deltas(3)),trPerDelta,1);
    trt_PTC_LoCohX2_mod3 = [d1_hicohPrec; d2_hicohPrec; d3_hicohPrec]; % list of the trials selected for this sample, to be used for indexing behaviora data
    % Ves only x2 mod 1 only
    d1_VesOnlyPrec = randsample(find(PTvesOnly_mod1X2' & mod==3 & data.delta==deltas(1)),trPerDelta,1);
    d2_VesOnlyPrec = randsample(find(PTvesOnly_mod1X2' & mod==3 & data.delta==deltas(2)),trPerDelta,1);
    d3_VesOnlyPrec = randsample(find(PTvesOnly_mod1X2' & mod==3 & data.delta==deltas(3)),trPerDelta,1);
    trt_PT_VesOnlyX2_mod1 = [d1_VesOnlyPrec; d2_VesOnlyPrec; d3_VesOnlyPrec];
    % Hi Coh X2 mod 2 only
    d1_locohPrec = randsample(find(PTCHiX2_mod2' & mod==3 & data.delta==deltas(1)),trPerDelta,1);
    d2_locohPrec = randsample(find(PTCHiX2_mod2' & mod==3 & data.delta==deltas(2)),trPerDelta,1);
    d3_locohPrec = randsample(find(PTCHiX2_mod2' & mod==3 & data.delta==deltas(3)),trPerDelta,1);
    trt_PTC_HiCohX2_mod2 = [d1_locohPrec; d2_locohPrec; d3_locohPrec]; % list of the trials selected for this sample, to be used for indexing behaviora data
    % Hi Coh X2 mod 3 only
    d1_hicohPrec = randsample(find(PTCHiX2_mod3' & mod==3 & data.delta==deltas(1)),trPerDelta,1);
    d2_hicohPrec = randsample(find(PTCHiX2_mod3' & mod==3 & data.delta==deltas(2)),trPerDelta,1);
    d3_hicohPrec = randsample(find(PTCHiX2_mod3' & mod==3 & data.delta==deltas(3)),trPerDelta,1);
    trt_PTC_HiCohX2_mod3 = [d1_hicohPrec; d2_hicohPrec; d3_hicohPrec]; % list of the trials selected for this sample, to be used for indexing behaviora data
    
    % (control comparison)
%     precX2_cntrl = PTCLoX2_mod2 | PTCLoX2_mod3 | PTCHiX2_mod2 | PTCHiX2_mod3;
    d1_anyPrec = randsample(find(precX2_cntrl' & mod==3 & data.delta==deltas(1)),trPerDelta,1);
    d2_anyPrec = randsample(find(precX2_cntrl' & mod==3 & data.delta==deltas(2)),trPerDelta,1);
    d3_anyPrec = randsample(find(precX2_cntrl' & mod==3 & data.delta==deltas(3)),trPerDelta,1);
    trt_AnyPrecX2 = [d1_anyPrec; d2_anyPrec; d3_anyPrec]; % list of the trials selected for this sample, to be used for indexing behaviora data

    fn = fieldnames(data);
    for f = 1:length(fn)
        PTC_LoCohDataX2_mod2.(fn{f}) = data.(fn{f})(trt_PTC_LoCohX2_mod2);
        PTC_LoCohDataX2_mod3.(fn{f}) = data.(fn{f})(trt_PTC_LoCohX2_mod3);
        PTC_VesOnlyDataX2_mod1.(fn{f}) = data.(fn{f})(trt_PT_VesOnlyX2_mod1);
        PTC_HiCohDataX2_mod2.(fn{f}) = data.(fn{f})(trt_PTC_HiCohX2_mod2);
        PTC_HiCohDataX2_mod3.(fn{f}) = data.(fn{f})(trt_PTC_HiCohX2_mod3);
        anyPrecDataX2.(fn{f}) = data.(fn{f})(trt_AnyPrecX2);
    end

    % Calculate and store weight for each boostrap iteration of trials for each PTC condition
    % weights following Low coherence trials
    gfit_locohPrec_mod2 = dots3DMP_fit_cgauss_NN(PTC_LoCohDataX2_mod2,mods,cohs,deltas);
    wVes_loCohPX2_mod2(b,:) = dots3DMP_wgts_thres_NN(gfit_locohPrec_mod2.muPMF,gfit_locohPrec_mod2.sigmaPMF,cohs,deltas);
    gfit_locohPrec_mod3 = dots3DMP_fit_cgauss_NN(PTC_LoCohDataX2_mod3,mods,cohs,deltas);
    wVes_loCohPX2_mod3(b,:) = dots3DMP_wgts_thres_NN(gfit_locohPrec_mod3.muPMF,gfit_locohPrec_mod3.sigmaPMF,cohs,deltas);
    % cue weight following High coherence trials
    gfit_HicohPrec_mod2 = dots3DMP_fit_cgauss_NN(PTC_HiCohDataX2_mod2,mods,cohs,deltas);
    wVes_HiCohPX2_mod2(b,:) = dots3DMP_wgts_thres_NN(gfit_HicohPrec_mod2.muPMF,gfit_HicohPrec_mod2.sigmaPMF,cohs,deltas);
    gfit_HicohPrec_mod3 = dots3DMP_fit_cgauss_NN(PTC_HiCohDataX2_mod3,mods,cohs,deltas);
    wVes_HiCohPX2_mod3(b,:) = dots3DMP_wgts_thres_NN(gfit_HicohPrec_mod3.muPMF,gfit_HicohPrec_mod3.sigmaPMF,cohs,deltas);
    % weights following Mod 1 x2 trials
    gfit_VesOnlyPrec_mod1X2 = dots3DMP_fit_cgauss_NN(PTC_VesOnlyDataX2_mod1,mods,cohs,deltas);
    wVes_VesOnlyPX2_mod1(b,:) = dots3DMP_wgts_thres_NN(gfit_VesOnlyPrec_mod1X2.muPMF,gfit_VesOnlyPrec_mod1X2.sigmaPMF,cohs,deltas);
    % weights for control trials ---- in general following 2 consecutive same coh trials anyX2
    gfit_anyPrec = dots3DMP_fit_cgauss_NN(anyPrecDataX2,mods,cohs,deltas);
    wVes_anyCohPX2(b,:) = dots3DMP_wgts_thres_NN(gfit_anyPrec.muPMF,gfit_anyPrec.sigmaPMF,cohs,deltas);
end

% Avg cue weight across all boostraps, for trials of each PTC condition
% PTC Low
wVes_loCohPrec_boostrapAvg_mod2 = mean(wVes_loCohPX2_mod2,1);
wVes_loCohPrec_boostrapAvg_mod3 = mean(wVes_loCohPX2_mod3,1);
% PTC High
wVes_HiCohPrec_boostrapAvg_mod2 = mean(wVes_HiCohPX2_mod2,1);
wVes_HiCohPrec_boostrapAvg_mod3 = mean(wVes_HiCohPX2_mod3,1);
% Mod 1 x2
wVes_VesOnlyP_boostrapAvg_mod1X2= mean(wVes_VesOnlyPX2_mod1,1);
% Control PTC
wVes_anyCohPrec_boostrapAvg = mean(wVes_anyCohPX2,1);

%     % PTC Low stats
% %vestibular weight comparisons 
% wVes_loPmod2vloPmod3 = wVes_loCohPrec_boostrapAvg_mod2 > wVes_loCohPrec_boostrapAvg_mod3 % mod 2 vis uncertainty = more elevated vesW than mod3 vis uncertainty?
% wVes_loPmod2vAnyPX2 = wVes_loCohPrec_boostrapAvg_mod2 > wVes_anyCohPrec_boostrapAvg % Does vestibular cue weight increase following Low coh trial?
% wVes_loPmod3vAnyPX2 = wVes_loCohPrec_boostrapAvg_mod3 > wVes_anyCohPrec_boostrapAvg % Does vestibular cue weight increase following Low coh trial?
% % Statistical significance
% loPmod2vsloPmod3_testX2 = ttest2(wVes_loCohPX2_mod2,wVes_loCohPX2_mod3)
% loPmod2vsAnyP_ttestX2 = ttest2(wVes_loCohPX2_mod2,wVes_anyCohPX2)
% loPmod3vsAnyP_ttestX2 = ttest2(wVes_loCohPX2_mod3,wVes_anyCohPX2)
%     % PTC High stats
% %vestibular weight comparisons 
% wVes_HiPmod2vHiPmod3 = wVes_HiCohPrec_boostrapAvg_mod2 > wVes_HiCohPrec_boostrapAvg_mod3 % mod 2 vis uncertainty = more elevated vesW than mod3 vis uncertainty?
% wVes_HiPmod2vAnyPX2 = wVes_HiCohPrec_boostrapAvg_mod2 > wVes_anyCohPrec_boostrapAvg % Does vestibular cue weight increase folHiwing Hiw coh trial?
% wVes_HiPmod3vAnyPX2 = wVes_HiCohPrec_boostrapAvg_mod3 > wVes_anyCohPrec_boostrapAvg % Does vestibular cue weight increase folHiwing Hiw coh trial?
% % Statistical significance
% HiPmod2vsHiPmod3_testX2 = ttest2(wVes_HiCohPX2_mod2,wVes_HiCohPX2_mod3)
% HiPmod2vsAnyP_ttestX2 = ttest2(wVes_HiCohPX2_mod2,wVes_anyCohPX2)
% HiPmod3vsAnyP_ttestX2 = ttest2(wVes_HiCohPX2_mod3,wVes_anyCohPX2)

% ttest2(wVes_HiCohPX2_mod3(:,1),wVes_loCohPX2_mod3(:,1)) % Hi vs Low PTC result in significantly different vesW on Low coh trials

% Plotting vesWeights
% x-axis
%Hi Coh trials
xVesOnlyPrecHi_mod1X2 = ones(size(wVes_VesOnlyPX2_mod1(:,1)))*1;
xLoPrecHi_mod2 = ones(size(wVes_loCohPX2_mod2(:,1)))*1.1;
xLoPrecHi_mod3 = ones(size(wVes_loCohPX2_mod2(:,1)))*1.2;
xHiPrecHi_mod2 = ones(size(wVes_loCohPX2_mod2(:,1)))*1.12;
xHiPrecHi_mod3 = ones(size(wVes_loCohPX2_mod2(:,1)))*1.22;
xAnyPrecHi = ones(size(wVes_loCohPX2_mod2(:,1)))*1.3;
%Lo Coh Trials
xVesOnlyPrecLo_mod1X2 = ones(size(wVes_VesOnlyPX2_mod1(:,1)))*2;
xLoPrecLo_mod2 = ones(size(wVes_loCohPX2_mod2(:,1)))*2.1;
xLoPrecLo_mod3 = ones(size(wVes_loCohPX2_mod2(:,1)))*2.2;
xHiPrecLo_mod2 = ones(size(wVes_loCohPX2_mod2(:,1)))*2.12;
xHiPrecLo_mod3 = ones(size(wVes_loCohPX2_mod2(:,1)))*2.22;
xAnyPrecLo = ones(size(wVes_loCohPX2_mod2(:,1)))*2.3;
% Consolidate plotting variables
xHiCohTr = [xLoPrecHi_mod2 xLoPrecHi_mod3 xVesOnlyPrecHi_mod1X2 xHiPrecHi_mod2 xHiPrecHi_mod3 xAnyPrecHi];
xLoCohTr = [xLoPrecLo_mod2 xLoPrecLo_mod3 xVesOnlyPrecLo_mod1X2 xHiPrecLo_mod2 xHiPrecLo_mod3 xAnyPrecLo];

HighCohTrialsW = [wVes_loCohPX2_mod2(:,2) wVes_loCohPX2_mod3(:,2) wVes_VesOnlyPX2_mod1(:,2) wVes_HiCohPX2_mod2(:,2) wVes_HiCohPX2_mod3(:,2) wVes_anyCohPX2(:,2)];
LowCohTrialsW = [wVes_loCohPX2_mod2(:,1) wVes_loCohPX2_mod3(:,1) wVes_VesOnlyPX2_mod1(:,1) wVes_HiCohPX2_mod2(:,1) wVes_HiCohPX2_mod3(:,1) wVes_anyCohPX2(:,1)];

wVes_boostrpAvgs_LoCohs = [ wVes_loCohPrec_boostrapAvg_mod2(:,1) wVes_loCohPrec_boostrapAvg_mod3(:,1) wVes_VesOnlyP_boostrapAvg_mod1X2(:,1) wVes_HiCohPrec_boostrapAvg_mod2(:,1) wVes_HiCohPrec_boostrapAvg_mod3(:,1) wVes_anyCohPrec_boostrapAvg(:,1) ];
wVes_boostrpAvgs_HiCohs = [ wVes_loCohPrec_boostrapAvg_mod2(:,2) wVes_loCohPrec_boostrapAvg_mod3(:,2) wVes_VesOnlyP_boostrapAvg_mod1X2(:,2) wVes_HiCohPrec_boostrapAvg_mod2(:,2) wVes_HiCohPrec_boostrapAvg_mod3(:,2) wVes_anyCohPrec_boostrapAvg(:,2) ];

clr = {[1.00,0.40,0.40], [0.46,0.92,0.39],'b', 'r',[0.15,0.61,0.01],'k'};
% clr = {'b','k','r'}; % LRed, LGreen, Blue, DRed, DGreen, Blk
figure(); hold on;
for pl = 1:length(clr)
    plot(xHiCohTr(:,pl),HighCohTrialsW(:,pl),'.','Color',clr{pl},'MarkerSize',8)
    plot(xLoCohTr(:,pl),LowCohTrialsW(:,pl),'.','Color',clr{pl})
    plot(xHiCohTr(1,pl),wVes_boostrpAvgs_HiCohs(pl),'o','Color',clr{pl})
    plot(xLoCohTr(1,pl),wVes_boostrpAvgs_LoCohs(pl),'o','Color',clr{pl})
end
hold off;

ax = gca;
ax.XTick = [1.15 2.15];
ax.XTickLabels = {'High coh trials','Low coh trials'};
title('Previous trial coherence and modality affect cue weighting'); ylabel('Vestibular cue weight');
ax.XLim = [.75,2.45];

toc % 100 loops = 11mins run time 
%% Past trial VesOnly x1 split by High and Low confidence

for tr = 3:length(conf)
    if mod(tr-1)==mods(1) && conf(tr-1)==1
        vesOnlyX1_HighConf(tr) = true;
    else
        vesOnlyX1_HighConf(tr) = false;
    end
    if mod(tr-1)==mods(1) && conf(tr-1)==0
        vesOnlyX1_LoConf(tr) = true;
    else
        vesOnlyX1_LoConf(tr) = false;
    end
    if mod(tr-1)~=mod(tr-2) && coh(tr-1)~=coh(tr-2)
        precX2_cntrlMod1(tr) = true;
    else
        precX2_cntrlMod1(tr) = false;
    end
end


trPerDelta = 1000; % Heading and coherence of trial are sampled randomly sampled in each boostrapped dataset, relying on their even distribution in the data being sampled, while delta is specified due to its uneven sampling in the OG dataset
tic
for b = 1:50 % number of trial resamplings ('bootstraps') 

    % Ves only x2 mod 1 only
    d1_VesOnlyPrec_HiConf = randsample(find(vesOnlyX1_HighConf' & mod==3 & data.delta==deltas(1)),trPerDelta,1);
    d2_VesOnlyPrec_HiConf = randsample(find(vesOnlyX1_HighConf' & mod==3 & data.delta==deltas(2)),trPerDelta,1);
    d3_VesOnlyPrec_HiConf = randsample(find(vesOnlyX1_HighConf' & mod==3 & data.delta==deltas(3)),trPerDelta,1);
    trt_PT_VesOnly_HiConf = [d1_VesOnlyPrec_HiConf; d2_VesOnlyPrec_HiConf; d3_VesOnlyPrec_HiConf];
    % Ves only x1 mod 1 only
    d1_VesOnlyPrec_LoConf = randsample(find(vesOnlyX1_LoConf' & mod==3 & data.delta==deltas(1)),trPerDelta,1);
    d2_VesOnlyPrec_LoConf = randsample(find(vesOnlyX1_LoConf' & mod==3 & data.delta==deltas(2)),trPerDelta,1);
    d3_VesOnlyPrec_LoConf = randsample(find(vesOnlyX1_LoConf' & mod==3 & data.delta==deltas(3)),trPerDelta,1);
    trt_PT_VesOnly_LoConf = [d1_VesOnlyPrec_LoConf; d2_VesOnlyPrec_LoConf; d3_VesOnlyPrec_LoConf];
    % (control comparison)
    d1_anyPrec = randsample(find(precX2_cntrlMod1' & mod==3 & data.delta==deltas(1)),trPerDelta,1);
    d2_anyPrec = randsample(find(precX2_cntrlMod1' & mod==3 & data.delta==deltas(2)),trPerDelta,1);
    d3_anyPrec = randsample(find(precX2_cntrlMod1' & mod==3 & data.delta==deltas(3)),trPerDelta,1);
    trt_AnyPrecX2 = [d1_anyPrec; d2_anyPrec; d3_anyPrec]; % list of the trials selected for this sample, to be used for indexing behaviora data

    fn = fieldnames(data);
    for f = 1:length(fn)
        PTC_VesOnlyData_HiConf.(fn{f}) = data.(fn{f})(trt_PT_VesOnly_HiConf);
        PTC_VesOnlyData_LoConf.(fn{f}) = data.(fn{f})(trt_PT_VesOnly_LoConf); % Previous trial was VesO and Low Confidence
        anyPrecDataX2.(fn{f}) = data.(fn{f})(trt_AnyPrecX2);
    end

% Calculate and store weight for each boostrap iteration of trials for each PTC condition
    % weights following Mod 2 trials

    % weights following Mod 1 x2 trials
    gfit_VesOnlyPrec_HiConf = dots3DMP_fit_cgauss_NN(PTC_VesOnlyData_HiConf,mods,cohs,deltas);
    wVes_VesOnlyP_HiConf(b,:) = dots3DMP_wgts_thres_NN(gfit_VesOnlyPrec_HiConf.muPMF,gfit_VesOnlyPrec_HiConf.sigmaPMF,cohs,deltas);
    % x1
    gfit_VesOnlyPrec_LoConf = dots3DMP_fit_cgauss_NN(PTC_VesOnlyData_LoConf,mods,cohs,deltas);
    wVes_VesOnlyPX1_LoConf(b,:) = dots3DMP_wgts_thres_NN(gfit_VesOnlyPrec_LoConf.muPMF,gfit_VesOnlyPrec_LoConf.sigmaPMF,cohs,deltas);
    % weights following assortment of trials (control)
    gfit_anyPrec = dots3DMP_fit_cgauss_NN(anyPrecDataX2,mods,cohs,deltas);
    wVes_anyCohPX2(b,:) = dots3DMP_wgts_thres_NN(gfit_anyPrec.muPMF,gfit_anyPrec.sigmaPMF,cohs,deltas);
end

% Avg cue weight across all boostraps, for trials of each PTC condition

% Mod 1 HiConf
wVes_VesOnlyP_boostrapAvg_HiConf= mean(wVes_VesOnlyP_HiConf,1);
% Mod 1 LoConf
wVes_VesOnlyP_boostrapAvg_LoConf= mean(wVes_VesOnlyPX1_LoConf,1);
% Control PT
wVes_anyCohPrec_boostrapAvg = mean(wVes_anyCohPX2,1);

% Plotting vesWeights
% x-axis
% Prec Hi

xVesOnlyPrecHi_mod1X2 = ones(size(wVes_VesOnlyPX1_LoConf(:,1)))*1.1;
xVesOnlyPrecHi_mod1X1 = ones(size(wVes_VesOnlyPX1_LoConf(:,1)))*1;
xAnyPrecHi = ones(size(wVes_VesOnlyPX1_LoConf(:,1)))*1.2;
% Prec Lo
xVesOnlyPrecLo_mod1X2 = ones(size(wVes_VesOnlyPX1_LoConf(:,1)))*2.1;
xVesOnlyPrecLo_mod1X1 = ones(size(wVes_VesOnlyPX1_LoConf(:,1)))*2;
xAnyPrecLo = ones(size(wVes_VesOnlyPX1_LoConf(:,1)))*2.2;
% Consolidate plotting variables
% X-axis position
xHiCohTr = [xVesOnlyPrecHi_mod1X1 xVesOnlyPrecHi_mod1X2 xAnyPrecHi];
xLoCohTr = [xVesOnlyPrecLo_mod1X1 xVesOnlyPrecLo_mod1X2 xAnyPrecLo];

HighCohTrialsW = [ wVes_VesOnlyPX1_LoConf(:,2) wVes_anyCohPX2(:,2) wVes_VesOnlyP_HiConf(:,2)];
LowCohTrialsW = [wVes_VesOnlyPX1_LoConf(:,1) wVes_anyCohPX2(:,1) wVes_VesOnlyP_HiConf(:,1)];

wVes_boostrpAvgs_LoCohs = [ wVes_VesOnlyP_boostrapAvg_LoConf(:,1) wVes_anyCohPrec_boostrapAvg(:,1) wVes_VesOnlyP_boostrapAvg_HiConf(:,1) ];
wVes_boostrpAvgs_HiCohs = [ wVes_VesOnlyP_boostrapAvg_LoConf(:,2) wVes_anyCohPrec_boostrapAvg(:,2) wVes_VesOnlyP_boostrapAvg_HiConf(:,2) ];

clr = {'c', 'k','b'}; % LBlue, DBlue, Blk
figure(); hold on;
for pl = 1:length(clr)
    plot(xHiCohTr(:,pl),HighCohTrialsW(:,pl),'.','Color',clr{pl},'MarkerSize',8)
    plot(xLoCohTr(:,pl),LowCohTrialsW(:,pl),'.','Color',clr{pl})
    plot(xHiCohTr(1,pl),wVes_boostrpAvgs_HiCohs(pl),'o','Color',clr{pl})
    plot(xLoCohTr(1,pl),wVes_boostrpAvgs_LoCohs(pl),'o','Color',clr{pl})
end
hold off;

ax = gca;
ax.XTick = [1.1 2.1];
ax.XTickLabels = {'High coh trials','Low coh trials'};
title('Previous trial coherence and modality affect cue weighting x2'); ylabel('Vestibular cue weight');
ax.XLim = [.75,2.45];

toc % 100 loops = 11mins run time 

%%

addpath('C:\Users\yhaile2\Documents\CODE\MATLAB\3DMP\SensoryExpectation\Functions_CueWeights_PastTrialEffects')
tr = 3:length(conf);
vesOnlyX1_HighConf = mod(tr-1)==mods(1) & conf(tr-1)==1 ;
vesOnlyX1_LoConf = mod(tr-1)==mods(1) & conf(tr-1)==0;
precX2_cntrlMod1 = mod(tr-1)~=mod(tr-2) & coh(tr-1)~=coh(tr-2);

LI = [vesOnlyX1_HighConf precX2_cntrlMod1 vesOnlyX1_LoConf];
pastTrialEffectsOnCueW(tr,data,LI,1000,5)

% PT mod1 x1 vs x2
tr = 3:length(mod);

vesOnlyx2 = mod(tr-1) == mods(1) & mod(tr-2) == mods(1);
vesOnlyx1 = mod(tr-1) == mods(1);
vesOnlyx1_notx2 = mod(tr-1) == mods(1) & mod(tr-2) ~= mods(1);
precX2_cntrlMod1 = mod(tr-1)~=mod(tr-2) & coh(tr-1)~=coh(tr-2);

LI = [precX2_cntrlMod1 vesOnlyx1_notx2 vesOnlyx1 vesOnlyx2];
[vesW vesW_avg] = pastTrialEffectsOnCueW(tr,data,LI,1000,5)




