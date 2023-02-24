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
useVelAcc = 0 ;

plotExampleTrials = 1;

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

%% 
% set some parameters manually to get reasonable sim outputs

% sensitivity constants for converting heading angle into mean of momentary evidence
ks = 25; % scale factor, for quickly testing different k levels
  % set manually to get reasonable results from images_dtb_2d
kves = ks;
kvis = ks * [2 4.5]/3; % straddles vestibular reliability, by construction

% threshold for high bet in logOdds, ignored if conftask==1
theta = 0.8; 


sigmaVes = 0.01;        % std of momentary evidence
sigmaVis = [0.01 0.01]; % allow for separate sigmas for condition, coherence
    % set manually to get reasonable looking dv trajectories;
    % also affects peakiness/flatness of confidence curve

B = 1.2; % bound height (also hand-tuned)

% draw non-decision times from Gaussian dist
muTnd = [0.3 0.3 0.3]; sdTnd = 0.06; % in seconds
%% use parameters from model fits to data

subject = 'lucio';

kves = 0.23;
kvis = [0.15 0.32];

B = 0.33;
muTnd  = [0.49 0.69 0.56];
confLapse = [0.07 0.38 0.09];
theta = 0.07;

% these are marginalized out of model fitting (?), need them here
sigmaVes = 0.01;
sigmaVis = [0.01 0.01];
sdTnd  = 0; % fixed SD

% kves and kvis were scaled down for minimization, but scaled up 
kves = kves.*100;
kvis = kvis.*100;
%%

% assume the log odds mapping is based on an equal amount of experience with the
% *three* levels of reliability (ves, vis-low, vis-high) hence k and sigma
% are their averages
k = mean([kves kvis]);


% uses method of images to compute PDF of the 2D DV, as in van den Berg et
% al. 2016 (similar to Kiani et al. 2014)
% lose_flag and plot_flag slow things down. 
% lose_flag is needed for logOdds mapping
R.t = 0.001:0.001:duration/1000;
R.Bup = B;
R.drift = k * sind(hdgs(hdgs>=0)); % takes only unsigned drift rates
R.lose_flag = 1;
R.plotflag  = 1; % 1 = plot, 2 = plot and export_fig
P =  images_dtb_2d(R);


if useVelAcc
    % create acceleration and velocity profiles for scaling drift rate - SJ 04/2020
    % cf. Drugowitsch et al. 2014
    % cf. Hou et al. 2019, peak vel = 0.37m/s, SD = 210ms
    % vel = normpdf(1:duration,duration/2,210);
    % vel = 0.37*vel./max(vel);
    % acc = gradient(vel)*1000; % multiply by 1000 to get from m/s/ms to m/s/s
    
    % FL settings...160cm, sigma=0.14 (programmed)
    ampl = 0.16; % movement in metres
    pos = normcdf(1:duration,duration/2,0.14*duration)*ampl;
    vel = gradient(pos)*1000; % metres/s
    acc = gradient(vel);
    
    % normalize
    vel = vel./max(vel);
    acc = abs(acc./max(acc)); % (and abs)
    
else
    vel = ones(1,duration);
    acc = vel;
end


% store parameter settings for bookkeeping
origParams.kves = kves;
origParams.kvis = kvis;
origParams.sigmaVes = sigmaVes;
origParams.sigmaVis = sigmaVis;
origParams.B = B;
origParams.Tnd = muTnd;
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

Tnds = muTnd + randn(ntrials,1).*sdTnd;

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
    
    if length(muTnd)==1
        Tnd = Tnds(n); % Tnd for nth trial, in seconds
    else
        Tnd = Tnds(n,modality(n));
    end
    
    switch modality(n)
        case 1
        
            mu = acc .* kves * sind(hdg(n)) / 1000; % mean of momentary evidence
                % (I'm guessing drift rate in images_dtb is per second, hence div by 1000)
            s = [sigmaVes sigmaVes]; % standard deviaton vector (see below)
        case 2
            
            mu = vel .* kvis(cohs==coh(n)) * sind(hdg(n)) / 1000;
            s = [sigmaVis(cohs==coh(n)) sigmaVis(cohs==coh(n))];
        case 3
            
            % positive delta defined as ves to the left, vis to the right
            muVes = acc .* kves               * sind(hdg(n)-delta(n)/2) / 1000;
            muVis = vel .* kvis(cohs==coh(n)) * sind(hdg(n)+delta(n)/2) / 1000;

            % optimal weights (Drugo et al.)
            wVes = sqrt( kves^2 / (kvis(cohs==coh(n))^2 + kves^2) );
            wVis = sqrt( kvis(cohs==coh(n))^2 / (kvis(cohs==coh(n))^2 + kves^2) );

%             wVes = rand; wVis = 1 - wVes;

            mu = wVes.*muVes + wVis.*muVis;
            
            % clearly brain does not have access to optimal weights at
            % outset of trial, but must come up with some heuristic version
            % of these based on inferred reliability from accumulated
            % evidence? 

            % the DV is a sample from a dist with mean = weighted sum of
            % means. thus the variance is the weighted sum of variances
            % (error propagation formula):
            sigmaComb = sqrt(wVes.^2 .* sigmaVes^2 + wVis.^2 .* sigmaVis(cohs==coh(n))^2); % assume zero covariance
            s = [sigmaComb sigmaComb];
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

    % store each trial's dv, if desired...
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

    % use map to look up log-odds that the motion is rightward
    diffV = abs((P.y+B)-finalV(n));
    diffT = abs(R.t-RT(n));

    thisV = find(diffV==min(diffV));
    thisT = find(diffT==min(diffT));
    logOddsCorr(n) = P.logOddsCorrMap(thisV(1), thisT(1));

    if conftask==1 % sacc endpoint
        expectedPctCorr(n) = logistic(logOddsCorr(n)); % convert to pct corr
        conf(n) = 2*expectedPctCorr(n) - 1; % convert to 0..1
    elseif conftask==2 % PDW
        conf(n) = logOddsCorr(n) > theta;
    end
    if isnan(conf(n)), conf(n)=0; end % if dvs are almost overlapping, force conf to zero as it can sometimes come out as NaN

    RT(n) = RT(n) + Tnd; % add NDT, already in seconds!
    
    if plotExampleTrials

        if modality(n)==1 && hdg(n)==3 && choice(n)==1 && doneWith1==0 % make a better plot for talk/poster;
                              % must not shuffle trial list for this to go in correct in order

            figure(1000); set(gcf, 'Color', [1 1 1], 'Position', [100 100 350 375], 'PaperPositionMode', 'auto');
            hold on; box off;
            xlabel('Time (ms)');
            ylabel('Accum. evidence');
            ylim([-1.1 1.1].*B); xlim([0 duration]);
            %             set(gca,'yTick',-0.5:0.5:2);
            set(gca,'xTick',0:500:2000);
            changeAxesFontSize(gca,18,18);
            
            plot(dv(1:round((RT(n)-Tnd)*1000),1),'k-','LineWidth',2); hold on; % winning (R, since we preselect for R choice)
            plot(dv(1:round((RT(n)-Tnd)*1000),2),'color',[1 1 1]*0.5,'LineWidth',1); hold on; % losing (L)
            fprintf('Vest hdg = %.1f, RT (-NDT) = %.2f, conf = %.2f',hdg(n),RT(n)-Tnd,conf(n));

            plot(1:length(dv),ones(1,length(dv))*B,'k-','LineWidth',4);
            if conftask==2, text(50, B+0.02, sprintf('\\Theta = %.g',theta),'verti','bottom','fontsize',14,'fontweight','bold'); end

%             plot([round((RT(n)-Tnd)*1000) round((RT(n)-Tnd)*1000)],dv(round((RT(n)-Tnd)*1000),:),'k:')
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
        if modality(n)==2 && hdg(n)==3 && coh(n)==cohs(2) && choice(n)==1 && doneWith2==0
            figure(1000); hold on; 
            plot(dv(1:round((RT(n)-Tnd)*1000),1),'r-','LineWidth',2); %hold on; % winning (R, since we preselect for R choice)
            plot(dv(1:round((RT(n)-Tnd)*1000),2),'m-','LineWidth',1); % losing (L)            doneWith2 = 1;
%             plot([round((RT(n)-Tnd)*1000) round((RT(n)-Tnd)*1000)],dv(round((RT(n)-Tnd)*1000),:),'r:')
            fprintf('Vis hdg = %.1f, RT (-NDT) = %.2f, conf = %.2f',hdg(n),RT(n)-Tnd,conf(n));
            doneWith2 = 1; n
        end
        if modality(n)==3 && hdg(n)==3 && coh(n)==cohs(2) && choice(n)==1 && delta(n)==0 && doneWith3==0
            figure(1000); hold on; 
            plot(dv(1:round((RT(n)-Tnd)*1000),1),'b-','LineWidth',2); %hold on; % winning (R, since we preselect for R choice)
            plot(dv(1:round((RT(n)-Tnd)*1000),2),'c-','LineWidth',1); % losing (L)
%             plot([round((RT(n)-Tnd)*1000) round((RT(n)-Tnd)*1000)],dv(round((RT(n)-Tnd)*1000),:),'b:')
            fprintf('Comb hdg = %.1f, RT (-NDT) = %.2f, conf = %.2f',hdg(n),RT(n)-Tnd,conf(n));
            doneWith3 = 1; n
%             export_fig('threeTrials','-eps');
        end
    end

end
toc

% if plotExampleTrials
% figure(1000);
% text(duration*0.75,0+B/2+0.2,'Ves','fontsize',14,'horizo','center');
% text(duration*0.75,0+B/2,'Vis','fontsize',14,'horizo','center');
% text(duration*0.75,0+B/2-0.2,'Comb','fontsize',14,'horizo','center');


% yl = get(gca,'ylim');
% if useVelAcc,plot(1:length(vel),vel+yl(1),'k--','linew',0.5); end
% end

choice(choice==0) = sign(randn); % not needed under usual circumstances

% sanity check:
correct = (choice==1 & hdg>0) | (choice==-1 & hdg<0) | ...
    (rand<0.5 & (hdg==0 | abs(hdg)<abs(delta)));

pCorrect = sum(correct) / ntrials



%% format data like the real expt

choice(choice==1) = 2; choice(choice==-1) = 1; % 1=left, 2=right

data.modality = modality;
data.heading = hdg;
data.coherence = coh;
data.delta = delta;
data.choice = choice;
data.RT = RT; % already in seconds
data.conf = conf;
data.correct = correct;


if conftask==2
    data.PDW=data.conf;
end
% subject = 'simul';

cd('/Users/stevenjerjian/Desktop/FetschLab/Analysis')
save(sprintf('2DAccSim_lucioParams.mat',conftask,ntrials),'data','cohs','deltas','hdgs','mods','origParams','RTtask','conftask','subject')

% cd('/Users/stevenjerjian/Desktop/FetschLab/Analysis')
% save(sprintf('2DAccSim_conftask%d_%dtrs.mat',conftask,ntrials),'P','R','data','cohs','deltas','hdgs','mods','origParams','RTtask','conftask','subject')

%% plots
mods   = unique(data.modality);
cohs   = unique(data.coherence);
deltas = unique(data.delta);
hdgs   = unique(data.heading);

% means per condition, logistic fits
parsedData = dots3DMP_parseData(data,mods,cohs,deltas,hdgs,conftask,RTtask);

% gaussian fits
gfit = dots3DMP_fit_cgauss(data,mods,cohs,deltas,conftask,RTtask);

% plot it
%     dots3DMP_plots(parsedData,mods,cohs,deltas,hdgs,conftask,RTtask)
dots3DMP_plots_cgauss_byCoh(gfit,parsedData,mods,cohs,deltas,hdgs,conftask,RTtask)

%% now try fitting the fake data to recover the generative parameters
% 
% % options.errfun = 'dots3DMP_fit_2Dacc_err_nSims';
options.errfun = 'dots3DMP_fit_2Dacc_err_sepbounds_noSim';

% % options.nreps  = 100;
% % options.confModel = 'evidence+time';


% choose whether to run fit with interpolated headings

% this is sort of redundant  for now, because model fits are
% generated via Monte Carlo and are going to be too noisy for a nice
% interpolated fit

% SJ 10/2021, no longer doing model fits via Monte Carlo
options.runInterpFit = 1; 

options.fitMethod = 'fms'; %'fms','global','multi','pattern','bads'
% options.fitMethod = 'global';
% options.fitMethod = 'multi';
% options.fitMethod = 'pattern';
% options.fitMethod = 'bads';

% initial guess (or hand-tuned params)
kmult   = 30;
kvis    = kmult.*cohs';
kves    = mean(kvis);
BVes    = 0.6;
BVis    = 1.2;
BComb   = 0.8;
Tnd     = 300;
Ttc     = 0; % time to confidence!

fixed   = [0 1 1 1 1 1 1 1];
guess   = [kves kvis(1:2) BVes BVis BComb Tnd Ttc];

if conftask==2 % PDW
    theta = 0.6;

    fixed   = [0 1 1 1 1 1 1 1 1];

    guess   = [guess theta];
end

% ************************************
% set all fixed to 1 for hand-tuning, or 0 for full fit
fixed(:)=1;
% ************************************

% plot error trajectory (prob doesn't work with parallel fit methods)
options.ploterr  = 1;
options.dummyRun = 0;
options.RTtask   = RTtask;
options.conftask = conftask; % 1 - sacc endpoint, 2 - PDW

if options.ploterr, options.fh = 400; end

% remove 
% removethese = data.RT == max(data.RT);
% fnames = fieldnames(data);
% 
% for f=1:length(fnames)
%     data.(fnames{f})(removethese) = [];
% end

[X, err_final, fit, fitInterp] = dots3DMP_fitDDM(data,options,guess,fixed);

% plot it!
dots3DMP_plots_fit_byCoh(data,fitInterp,conftask,RTtask,0)
