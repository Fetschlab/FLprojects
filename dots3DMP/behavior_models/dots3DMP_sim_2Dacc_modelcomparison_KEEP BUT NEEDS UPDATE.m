% simulate "CCC" (cue-combination + confidence) experiment
% specifically in the case of the vis-ves heading task

% CF started it 2016
% heavily modified early 2019
% switched to 2D accumulator model 06/2020

% build expt and hand-pick some model params

clear
close all

RTtask = 1;
conftask = 2; % 1 - sacc endpoint, 2 - PDW

confModel     = 3; % 1 - time_only, 2 - evidence_only, 3 - evidence+time
separate_k    = 0; % 0 - kcomb = optimal kves/kvis comb, 1 - kcomb is a free parameter
cue_weighting = 0; % 0 - wVes = wVis = 0.5, -1 - random weights on each trial, 1 - wVes and wVis according to reliabilites
tempweighting = 3; % 0 - no weighting, 1 - weight by vel, 2 - weight by acc, 3 - vestibular by acc, visual by vel

plotExampleTrials = 0;

nreps = 100; % number of repetitions of each unique trial type
            % start small to verify it's working, then increase
            % (ntrials depends on num unique trial types)

cohs = [0.4 0.8]; % visual coherence levels (these are really just labels, since k's are set manually)
% hdgs = [-10 -5 -2.5 -1.25 0 1.25 2.5 5 10]; % heading angles
% hdgs = [-10 -3.5 -1.25 1.25 3.5 10]; % heading angles
hdgs = [-12 -6 -3 -1.5 0 1.5 3 6 12];
deltas = [-3 0 3]; % conflict angle; positive means vis to the right
mods = [1 2 3]; % stimulus modalities: ves, vis, comb
duration = 2000; % stimulus duration (ms)

theta = 0.4; % threshold for high bet in logOdds, ignored if conftask==1

kves  = 25;
kvis  = [15 40];
kcomb = [15 40];
sigmaVes = 0.02;
sigmaVis = [0.02 0.02];
sigmaComb = [0.02 0.02];
BVes     = 0.7; % don't accept negative bound heights
BVis     = 1.2; % fixed across cohs
BComb    = 1.0;
muTndVes = 300;
muTndVis = 300; % fixed across cohs
muTndComb = 300;

sdTnd = 60; % fixed SD
% Tnds = muTnd + randn(ntrials,1).*sdTnd; % fixed for all sims of given trial

% assume the mapping is based on an equal amount of experience with the 
% *three* levels of reliability (ves, vis-low, vis-high) hence k and sigma
% are their averages
k = mean([kves kvis]);

RVes.t = 0.001:0.001:duration/1000;
RVes.Bup = BVes;
RVes.drift = kves * sind(hdgs(hdgs>=0)); % takes only unsigned drift rates
RVes.lose_flag = 1;
RVes.plotflag = 0; % 1 = plot, 2 = plot and export_fig

PVes =  images_dtb_2d(RVes);

RVis.t = 0.001:0.001:duration/1000;
RVis.Bup = BVis;
RVis.drift = mean(kvis) * sind(hdgs(hdgs>=0)); % takes only unsigned drift rates
RVis.lose_flag = 1;
RVis.plotflag = 0; % 1 = plot, 2 = plot and export_fig

PVis =  images_dtb_2d(RVis);

RComb.t = 0.001:0.001:duration/1000;
RComb.Bup = BComb;
RComb.drift = mean(kcomb) * sind(hdgs(hdgs>=0)); % takes only unsigned drift rates
RComb.lose_flag = 1;
RComb.plotflag = 0; % 1 = plot, 2 = plot and export_fig

PComb =  images_dtb_2d(RComb);
% create acceleration and velocity profiles (arbitrary for now)
% SJ 04/2020
% Hou et al. 2019, peak vel = 0.37m/s, SD = 210ms
% 07/2020 lab settings...160cm in 1.3s, sigma=0.14

ampl = 0.16; % movement in metres
pos = normcdf(1:duration,duration/2,0.14*duration)*ampl;
vel = gradient(pos)*1000; % metres/s
acc = gradient(vel);

% vel = normpdf(1:duration,duration/2,210);
% vel = 0.37*vel./max(vel);
% acc = gradient(vel)*1000; % multiply by 1000 to get from m/s/ms to m/s/s

% normalize
vel = vel./max(vel);
acc = abs(acc./max(acc)); % (and abs)

origParams.kves = kves;
origParams.kvis = kvis;
origParams.kcomb = kcomb;
origParams.sigmaVes = sigmaVes;
origParams.sigmaVis = sigmaVis;
origParams.sigmaComb = sigmaComb;
origParams.BVes = BVes;
origParams.BVis = BVis;
origParams.BComb = BComb;
origParams.TndVes = muTndVes;
origParams.TndVis = muTndVis;
origParams.TndComb = muTndComb;
origParams.confModel = confModel; % 1 - time_only, 2 - evidence_only, 3 - evidence+time
origParams.separate_k    = separate_k; % 0 - kcomb = optimal kves/kvis comb, 1 - kcomb is a free parameter
origParams.cue_weighting = cue_weighting; % 0 - wVes = wVis = 0.5, -1 - random weights on each trial, 1 - wVes and wVis according to reliabilites
origParams.tempweighting = tempweighting; % 0 - no weighting, 1 - weight by vel, 2 - weight by acc, 3 - vestibular by acc, visual by vel

if conftask==2
    origParams.theta = theta; % PDW only
end

%% build trial list

[hdg, modality, coh, delta, ntrials] = dots3DMP_create_trial_list(hdgs,mods,cohs,deltas,nreps,0); % don't shuffle

% sample durations from truncated exponential?
    % not practical experimentally.
    % instead, make constant dur, assume RT task, or (equivalently, as far as
    % the simulation is concerned) fixed duration with internal bounded
    % process (Kiani et al. 2008)
dur = ones(ntrials,1) * duration;

% Tnds = muTnd + randn(ntrials,1).*sdTnd;

%% bounded evidence accumulation

% assume momentary evidence is proportional to sin(heading),
% as in drugowitsch et al 2014

% dv_all = cell(ntrials,1); % shouldn't need to store every trial's DV, but if you want to, it's here

choice = nan(ntrials,1);
RT = nan(ntrials,1);
finalV = nan(ntrials,1); % now this is the value of the losing accumulator
hitBound = zeros(1,ntrials);
logOddsCorr = nan(ntrials,1);
expectedPctCorr = nan(ntrials,1);
conf = nan(ntrials,1);

% for plotting
doneWith1 = 0;
doneWith2 = 0;
doneWith3 = 0;

% ME is now a draw from bivariate normal with mean vector Mu and covariance matrix V
% start with correlation matrix:
S = [1 -1/sqrt(2) ; -1/sqrt(2) 1];
% -1/sqrt(2) is the correlation for our version of images_dtb;
% can only be changed with an update to the flux file

tic
for n = 1:ntrials
%     Tnd = Tnds(n) / 1000; % Tnd for nth trial in seconds

    switch modality(n)
        case 1
            Tnd = (muTndVes + randn.*sdTnd)/1000;
            B = BVes; P = PVes; R = RVes;

            mu = acc .* kves * sind(hdg(n)) / 1000; % mean of momentary evidence
                % (I'm guessing drift rate in images_dtb is per second, hence div by 1000)
            s = [sigmaVes sigmaVes]; % standard deviaton vector (see below)
        case 2
            Tnd = (muTndVis + randn.*sdTnd)/1000;
            B = BVis; P = PVis; R = RVis;
            
            mu = vel .* kvis(cohs==coh(n)) * sind(hdg(n)) / 1000;
            s = [sigmaVis(cohs==coh(n)) sigmaVis(cohs==coh(n))];
        case 3
            Tnd = (muTndComb + randn.*sdTnd)/1000;
            B = BComb; P = PComb; R = RComb;
            
            
            if ~separate_k
                % positive delta defined as ves to the left, vis to the right
                
                switch tempweighting
                    case 0
                        tVes = ones(size(acc));
                        tVis = ones(size(acc));
                    case 1
                        tVes = vel;
                        tVis = vel;
                    case 2
                        tVes = acc;
                        tVis = acc;
                    case 3
                        tVes = acc;
                        tVis = vel;
                end
                    
                muVes = tVes .* kves               * sind(hdg(n)-delta(n)/2) / 1000;
                muVis = tVis .* kvis(cohs==coh(n)) * sind(hdg(n)+delta(n)/2) / 1000;
                
                if cue_weighting == 1
                    % optimal weights (Drugo et al.)
                    wVes = sqrt( kves^2 / (kvis(cohs==coh(n))^2 + kves^2) );
                    wVis = sqrt( kvis(cohs==coh(n))^2 / (kvis(cohs==coh(n))^2 + kves^2) );
                else
                    if cue_weighting == 0
                        wVes = 0.5;
                    else
                        wVes = rand;
                    end
                    wVis = 1 - wVes;
                end
                
                mu = wVes.*muVes + wVis.*muVis;
                
            else
                mu = vel .* kcomb(cohs==coh(n)) * sind(hdg(n)) / 1000;
            end
            
            
                
            % what about delta? Drugo task didn't have conflict condition
            % and what about temporal weighting? by acc or vel?
            
            if separate_k == 0
                % the DV is a sample from a dist with mean = weighted sum of
                % means. thus the variance is the weighted sum of variances
                % (error propagation formula):
                sigmaComb = sqrt(wVes.^2 .* sigmaVes^2 + wVis.^2 .* sigmaVis(cohs==coh(n))^2); % assume zero covariance
                s = [sigmaComb sigmaComb];
            else
                s = [sigmaComb(cohs==coh(n)) sigmaComb(cohs==coh(n))];
            end

    end

%     Mu = [mu,-mu]; % mean vector for 2D DV
    Mu = [mu; -mu]';

    % convert correlation to covariance matrix
    V = diag(s)*S*diag(s);
%     dv = [0 0; cumsum(mvnrnd(Mu,V,dur(n)-1))]; % bivariate normrnd
    dv = [0 0; cumsum(mvnrnd(Mu,V))]; % dv is now scaled by physical signal vel/acc

    % because Mu is signed according to heading (positive=right),
    % dv(:,1) corresponds to evidence favoring rightward, not evidence
    % favoring the correct decision (as in Kiani eqn. 3 and images_dtb)

%     dv_all{n} = dv;
    % decision outcome
    cRT1 = find(dv(:,1)>=B, 1);
    cRT2 = find(dv(:,2)>=B, 1);
    % the options are:
    % (1) only right accumulator hits bound,
    if ~isempty(cRT1) && isempty(cRT2)
        RT(n) = cRT1/1000;
        finalV(n) = dv(cRT1,2); % only 1 hit, so 2 is the loser
        hitBound(n) = 1;
        choice(n) = 1;
    % (2) only left accumulator hits bound,
    elseif isempty(cRT1) && ~isempty(cRT2)
        RT(n) = cRT2/1000;
        finalV(n) = dv(cRT2,1); % only 2 hit, so 1 is the loser
        hitBound(n) = 1;
        choice(n) = -1;
    % (3) neither hits bound,
    elseif isempty(cRT1) && isempty(cRT2)
        RT(n) = dur(n)/1000;
            % this is interesting: which DV matters for confidence
            % if neither hits bound? for now take the (abs) maximum,
            % ie whichever was closest to hitting bound. (alternative
            % would be their average?)
            % SJ 07/2020 finalV is fully determined by distance of loser from bound when winner
            % hits, so in this case, should also account for where winner is wrt bound
            % imagine winner 'did' hit bound, then where
            % would loser be relatively speaking (since logOdds map is
            % fixed)

        whichWon = dv(end,:)==max(dv(end,:));
        finalV(n) = dv(end,~whichWon) + B-dv(end,whichWon);
        % ^ effectively shifting the losing dv up by whatever the
        % difference is between the bound and the winning dv
%         finalV(n) = dv(end,~whichWon); % the not-whichWon is the loser
        % % finalV(n) = mean(dvEnds);
        finalV(n) = dv(end,~whichWon) + B-dv(end,whichWon); %

        hitBound(n) = 0;
        a = [1 -1];
        choice(n) = a(whichWon);
    % (4) or both do
    else
        RT(n) = min([cRT1 cRT2])/1000;
        whichWon = [cRT1<=cRT2 cRT1>cRT2];
        finalV(n) = dv(min([cRT1 cRT2]),~whichWon); % the not-whichWon is the loser
        hitBound(n) = 1;
        a = [1 -1];
        choice(n) = a(whichWon);
    end

    diffV = abs((P.y+B)-finalV(n));
    diffT = abs(R.t-RT(n));
            
    switch confModel
        case 1 % time only (RT-dependent)
            if contask==1
                conf(n) = 1 - RT(n) ./ range(P.t);
            elseif conftask==2
                conf(n) = RT(n) < theta;
            end
        case 2 % evidence only (losign accum dependent)
            if conftask==1
                conf(n) = max(diffV) ./ range(P.y);
            elseif conftask==2
                conf(n) = max(diffV) > theta;
            end
        case 3 % evidence and time
            % use map to look up log-odds that the motion is rightward
            
            thisV = find(diffV==min(diffV));
            thisT = find(diffT==min(diffT));
            logOddsCorr(n) = P.logOddsCorrMap(thisV(1), thisT(1));
            
            if conftask==1 % sacc endpoint
                expectedPctCorr(n) = logistic(logOddsCorr(n)); % convert to pct corr
                conf(n) = 2*expectedPctCorr(n) - 1; % convert to 0..1
            elseif conftask==2 % PDW
                conf(n) = logOddsCorr(n) > theta;
            end
       
        
    end
                  
    if isnan(conf(n)), conf(n)=0; end % if dvs are almost overlapping, force conf to zero as it can sometimes come out as NaN
    RT(n) = RT(n) + Tnd; % add NDT

    if plotExampleTrials

        if modality(n)==1 && hdg(n)==3 && choice(n)==1 && doneWith1==0 % make a better plot for talk/poster;
                              % must not shuffle trial list for this to go in correct in order

            figure(1000); set(gcf, 'Color', [1 1 1], 'Position', [100 100 350 375], 'PaperPositionMode', 'auto');
            hold on; box off;
            xlabel('Time (ms)');
            ylabel('Accum. evidence');
            ylim([-1.25 1.1].*B); xlim([0 duration]);
            %             set(gca,'yTick',-0.5:0.5:2);
            set(gca,'xTick',0:500:2000);
            changeAxesFontSize(gca,18,18);

            yl = get(gca,'ylim');
            plot(1:length(vel),vel+yl(1),'k--','linew',0.5)

            plot(dv(1:round((RT(n)-Tnd)*1000),1),'k-','LineWidth',2); hold on; % winning (R, since we preselect for R choice)
            plot(dv(1:round((RT(n)-Tnd)*1000),2),'color',[1 1 1]*0.5,'LineWidth',1); hold on; % losing (L)
            fprintf('Vest hdg = %.1f, RT (-NDT) = %.2f, conf = %.2f',hdg(n),RT(n)-Tnd,conf(n));

            plot(1:length(dv),ones(1,length(dv))*B,'k-','LineWidth',4);
            if conftask==2, text(1, B+0.02, sprintf('Theta = %.g',theta),'verti','bottom'); end

            plot([round((RT(n)-Tnd)*1000) round((RT(n)-Tnd)*1000)],dv(round((RT(n)-Tnd)*1000),:),'k:')
%             plot([round(RT(n)*1000) round(RT(n)*1000)],[dv(round((RT(n)-Tnd)*1000),2) B],'r--')
                
            box off;
            xlabel('Time (ms)');
%             ylabel('Accum. evidence for leftward');
            %ylim([-0.4*B B*1.1]); xlim([0 duration]);
%             set(gca,'yTick',-0.5:0.5:2);
            %set(gca,'xTick',0:500:2000);
            %changeAxesFontSize(gca,18,18);
            %export_fig('ves_acc2','-eps');

            doneWith1=1; n
        end
        if modality(n)==2 && hdg(n)==6 && coh(n)==cohs(2) && choice(n)==1 && doneWith2==0
            figure(1000); hold on; 
            plot(dv(1:round((RT(n)-Tnd)*1000),1),'r-','LineWidth',2); %hold on; % winning (R, since we preselect for R choice)
            plot(dv(1:round((RT(n)-Tnd)*1000),2),'m-','LineWidth',1); % losing (L)            doneWith2 = 1;
            plot([round((RT(n)-Tnd)*1000) round((RT(n)-Tnd)*1000)],dv(round((RT(n)-Tnd)*1000),:),'r:')
            fprintf('Vis hdg = %.1f, RT (-NDT) = %.2f, conf = %.2f',hdg(n),RT(n)-Tnd,conf(n));
            doneWith2 = 1; n
        end
        if modality(n)==3 && hdg(n)==1.5 && coh(n)==cohs(2) && choice(n)==1 && delta(n)==0 && doneWith3==0
            figure(1000); hold on; 
            plot(dv(1:round((RT(n)-Tnd)*1000),1),'b-','LineWidth',2); %hold on; % winning (R, since we preselect for R choice)
            plot(dv(1:round((RT(n)-Tnd)*1000),2),'c-','LineWidth',1); % losing (L)
            plot([round((RT(n)-Tnd)*1000) round((RT(n)-Tnd)*1000)],dv(round((RT(n)-Tnd)*1000),:),'b:')
            fprintf('Comb hdg = %.1f, RT (-NDT) = %.2f, conf = %.2f',hdg(n),RT(n)-Tnd,conf(n));
            doneWith3 = 1; n
%             export_fig('threeTrials','-eps');
        end
    end

end
toc

if plotExampleTrials
figure(1000); %yl = get(gca,'ylim');
plot(1:length(vel),vel+yl(1),'k--','linew',0.5)
end

choice(choice==0) = sign(randn); % not needed under usual circumstances

% sanity check:
pCorrect_total = (sum(choice==1 & hdg>0) + sum(choice==-1 & hdg<0)) / ntrials





%% format data like the real expt

choice(choice==1) = 2; choice(choice==-1) = 1; % 1=left, 2=right

data.modality = modality;
data.heading = hdg;
data.coherence = coh;
data.delta = delta;
data.choice = choice;
data.RT = RT; % already in seconds
data.conf = conf;

if conftask==2
    data.PDW=data.conf;
end
subject = 'simul';

% cd('/Users/stevenjerjian/Desktop/FetschLab/Analysis')
% save(sprintf('2DAccSim_sepk_sepbounds_conftask%d_%dtrs.mat',conftask,ntrials),'P','R','data','cohs','deltas','hdgs','mods','origParams','RTtask','conftask','subject')

%% plots

    
dots3DMP_parseData
dots3DMP_plots
if 0
%%
dots3DMP_parseData_splitConf
dots3DMP_plots_splitConf


%% fit cumulative gaussians

dots3DMP_fit_cgauss
dots3DMP_plots_cgauss

end