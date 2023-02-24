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
confModel = 'evidence+time'; % 'evidence+time','evidence_only','time_only'
useVelAcc = 0; % toggle time-varying )vel/acc) vs. constant drift rate

plotExampleTrials = 0;

nreps = 1000; % number of repetitions of each unique trial type
              % start small to verify it's working, then increase
              % (ntrials depends on num unique trial types)

cohs = [0.4 0.8]; % visual coherence levels (these are really just labels, since k's are set manually)
hdgs = [-12 -6 -3 -1.5 -eps eps 1.5 3 6 12];
deltas = [-5 0 5]; % conflict angle; positive means vis to the right
mods = [1 2 3]; % stimulus modalities: ves, vis, comb
duration = 3000; % stimulus duration (ms)

lose_flag = 1;
plot_flag = 1;

%% PARAMS

kmult = 30; % try to reduce number of params
kvis  = kmult*cohs; % [20 40];
kves  = mean(kvis); % for now, assume straddling
% knoise = [0.07 0.07]; %??
sigmaVes = 0.05;
sigmaVis = [sigmaVes sigmaVes];
BVes     = 0.6;
BVis     = 1.2; % fixed across cohs
BComb    = 0.8;
muTnd    = 300; % fixed across mods SJ 10/11/2021 [unlike Drugo]
sdTnd = 0; % fixed SD

if conftask==2
    timeToConf = 0; % additional processing time for confidence
    theta = 2; % threshold for high bet in logOdds, ignored if conftask==1
else
    timeToConf = 0;
end
duration = duration + timeToConf;

origParams.kmult = kmult;
origParams.kvis = kvis;
origParams.kves = kves;
origParams.sigmaVes = sigmaVes;
origParams.sigmaVis = sigmaVis;
origParams.BVes = BVes;
origParams.BVis = BVis;
origParams.BComb = BComb;
origParams.muTnd = muTnd;
origParams.sdTnd = sdTnd;
origParams.ttc   = timeToConf;
if conftask==2
    origParams.theta = theta; % PDW only
end


%% calculate log odds corr maps using Wolpert MoI method

% assume the mapping is based on an equal amount of experience with the 
% *three* levels of reliability (ves, vis-low, vis-high) hence k and sigma
% are their averages, for the purpose of expected logOddsCorr
k = mean([kves kvis]);

R.t = 0.001:0.001:duration/1000;
R.Bup = B;
R.drift = k * sind(hdgs(hdgs>=0)); % takes only unsigned drift rates
R.lose_flag = lose_flag;
R.plotflag = plot_flag; % 1 = plot, 2 = plot and export_fig
P =  images_dtb_2d(R);

% RVes.t = 0.001:0.001:duration/1000;
% RVes.Bup = Bs(1);
% RVes.drift = k * sind(hdgs(hdgs>=0)); % takes only unsigned drift rates
% RVes.lose_flag = lose_flag;
% RVes.plotflag = plot_flag; % 1 = plot, 2 = plot and export_fig
% PVes =  images_dtb_2d(RVes);
% 
% RVis.t = 0.001:0.001:duration/1000;
% RVis.Bup = Bs(2);
% RVis.drift = k * sind(hdgs(hdgs>=0)); % takes only unsigned drift rates
% RVis.lose_flag = lose_flag;
% RVis.plotflag = plot_flag; % 1 = plot, 2 = plot and export_fig
% PVis =  images_dtb_2d(RVis);
% 
% RComb.t = 0.001:0.001:duration/1000;
% RComb.Bup = Bs(3);
% RComb.drift = k * sind(hdgs(hdgs>=0)); % takes only unsigned drift rates
% RComb.lose_flag = lose_flag;
% RComb.plotflag = plot_flag; % 1 = plot, 2 = plot and export_fig
% PComb =  images_dtb_2d(RComb);

RComb.t = 0.001:0.001:duration/1000;
RComb.Bup = BComb;
RComb.drift = k * sind(hdgs(hdgs>=0)); % takes only unsigned drift rates
RComb.lose_flag = lose_flag;
RComb.plotflag = plot_flag; % 1 = plot, 2 = plot and export_fig
PComb =  images_dtb_2d(RComb);

% try a single conf map
R.t = 0.001:0.001:duration/1000;
R.Bup = B;
R.drift = k * sind(hdgs(hdgs>=0)); % takes only unsigned drift rates
R.lose_flag = 1;
R.plotflag = 0; % 1 = plot, 2 = plot and export_fig
P =  images_dtb_2d(R);


%% create acceleration and velocity profiles (arbitrary for now)
% SJ 04/2020
% Hou et al. 2019, peak vel = 0.37m/s, SD = 210ms
% 07/2020 lab settings...160cm in 1.3s, sigma=0.14
ampl = 0.16; % movement in metres
pos = normcdf(1:duration,duration/2,0.14*duration)*ampl;
vel = gradient(pos)*1000; % metres/s
acc = gradient(vel);

% normalize (by max or by mean?) and take abs of acc 
vel = vel/mean(vel);
acc = abs(acc)/mean(abs(acc));
% vel = vel/max(vel);
% acc = abs(acc)/max(abs(acc));

if useVelAcc==0
%     vel = ones(size(vel))*mean(vel);
    vel = ones(size(vel));
    acc = vel;
end


%% build trial list

[hdg, modality, coh, delta, ntrials] = dots3DMP_create_trial_list(hdgs,mods,cohs,deltas,nreps,0); % don't shuffle

% make constant dur, assuming RT task, or (equivalently, as far as
% the simulation is concerned) fixed duration with internal bounded
% process (Kiani et al. 2008)
dur = ones(ntrials,1) * duration;

Tnds = muTnd + randn(ntrials,1).*sdTnd; % obs


%% bounded evidence accumulation

% assume momentary evidence is proportional to sin(heading),
% as in drugowitsch et al 2014

dv_all = cell(ntrials,1); % shouldn't need to store every trial's DV, but if you want to, it's here

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
    
    Tnd = Tnds(n) / 1000; % Tnd for nth trial in seconds

    switch modality(n)
        
        case 1
            %B = Bs(1); P = PVes; R = RVes;
            mu = acc .* kves * sind(hdg(n)) / 1000; % mean of momentary evidence
                % (I'm guessing drift rate in images_dtb is per second, hence div by 1000)
            s = [sigmaVes sigmaVes]; % standard deviaton vector (see below)
        case 2
%             B = BVis; P = PVis; R = RVis;
            mu = vel .* kvis(cohs==coh(n)) * sind(hdg(n)) / 1000;
            s = [sigmaVis(cohs==coh(n)) sigmaVis(cohs==coh(n))];
        case 3
%             B = BComb; P = PComb; R = RComb;
            % positive delta defined as ves to the left, vis to the right
            muVes = acc .* kves               * sind(hdg(n)-delta(n)/2) / 1000;
            muVis = vel .* kvis(cohs==coh(n)) * sind(hdg(n)+delta(n)/2) / 1000;
            % optimal weights (Drugo et al.)
            wVes = sqrt( kves^2 / (kvis(cohs==coh(n))^2 + kves^2) );
            wVis = sqrt( kvis(cohs==coh(n))^2 / (kvis(cohs==coh(n))^2 + kves^2) );
            
            % corrupt optimal weights with noise
%             wVes = wVes + randn*knoise(1);
%             wVis = wVis + randn*knoise(2);

%             wVes = rand; wVis = 1 - wVes;

            mu = wVes.*muVes + wVis.*muVis;
            
            % the DV is a sample from a dist with mean = weighted sum of
            % means. thus the variance is the weighted sum of variances
            % (error propagation formula):
            sigmaComb = sqrt(wVes.^2 .* sigmaVes^2 + wVis.^2 .* sigmaVis(cohs==coh(n))^2); % assume zero covariance
            s = [sigmaComb sigmaComb];
    end

    Mu = [mu; -mu]'; % mean vector for 2D DV

    % convert correlation to covariance matrix
    V = diag(s)*S*diag(s);
%     dv = [0 0; cumsum(mvnrnd(Mu,V,dur(n)-1))]; % bivariate normrnd
    dv = [0 0; cumsum(mvnrnd(Mu,V))]; % dv is now scaled by physical signal vel/acc

    % because Mu is signed according to heading (positive=right),
    % dv(:,1) corresponds to evidence favoring rightward, not evidence
    % favoring the correct decision (as in Kiani eqn. 3 and images_dtb)

    dv_all{n} = dv;
    % decision outcome
    cRT1 = find(dv(1:dur(n)-timeToConf,1)>=B, 1);
    cRT2 = find(dv(1:dur(n)-timeToConf,2)>=B, 1);
    
    % the options are:
    % (1) only right accumulator hits bound,
    if ~isempty(cRT1) && isempty(cRT2)
        RT(n) = cRT1/1000;
        finalV(n) = dv(cRT1+timeToConf,2); % only 1 hit, so 2 is the loser
        hitBound(n) = 1;
        choice(n) = 1;
    % (2) only left accumulator hits bound,
    elseif isempty(cRT1) && ~isempty(cRT2)
        RT(n) = cRT2/1000;
        finalV(n) = dv(cRT2+timeToConf,1); % only 2 hit, so 1 is the loser
        hitBound(n) = 1;
        choice(n) = -1;
    % (3) neither hits bound,
    elseif isempty(cRT1) && isempty(cRT2)
        RT(n) = (dur(n)-timeToConf)/1000;
        
            % which DV matters for confidence if neither hits bound? 
            % SJ 07/2020 logOddsCorrMap is fixed, so just shift finalV up
            % so that 'winner' did hit bound,
        whichWon = dv(dur(n)-timeToConf,:)==max(dv(dur(n)-timeToConf,:));
        finalV(n) = dv(end,~whichWon) + B-dv(end,whichWon);
        % ^  shifting the losing dv up by whatever the
        % difference is between the bound and the winning dv

        hitBound(n) = 0;
        a = [1 -1];
        choice(n) = a(whichWon);
    % (4) or both do
    else
        RT(n) = min([cRT1 cRT2])/1000;
        whichWon = [cRT1<=cRT2 cRT1>cRT2];
        finalV(n) = dv(min([cRT1 cRT2]+timeToConf),~whichWon); % the not-whichWon is the loser
        hitBound(n) = 1;
        a = [1 -1];
        choice(n) = a(whichWon);
    end
    
    diffV = abs((P.y+B)-finalV(n));
    diffT = abs(R.t-(RT(n)+timeToConf));
            
    switch confModel
        case 'evidence+time'
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
        case 'evidence_only'
            if conftask==1
                conf(n) = max(diffV) ./ range(P.y);
            elseif conftask==2
                conf(n) = max(diffV) > theta;
            end
        case 'time_only'
            if contask==1
                conf(n) = 1 - RT(n) ./ range(P.t);
            elseif conftask==2
                conf(n) = RT(n) < theta;
            end
    end
    
    if rand<confLapse(modality(n)) && conftask==2
        conf(n) = true;
    end
                         
    if isnan(conf(n)), conf(n)=0; end % if dvs are almost overlapping, force conf to zero as it can sometimes come out as NaN
    RT(n) = RT(n) + Tnd;

    if plotExampleTrials

        if modality(n)==1 && hdg(n)==3 && choice(n)==1 && doneWith1==0 % make a better plot for talk/poster;
                              % must not shuffle trial list for this to go in correct in order

            figure(1000); set(gcf, 'Color', [1 1 1], 'Position', [100 100 350 375], 'PaperPositionMode', 'auto');
            hold on; box off;
            xlabel('Time (ms)');
            ylabel('Accum. evidence');
            ylim([-1.25 1.5].*B); xlim([0 duration]);
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

choice(choice==0) = sign(randn(sum(choice==0),1)); % not needed under usual circumstances

% sanity check:
correct = (choice==1 & hdg>0) | (choice==-1 & hdg<0) | ...
    (rand<0.5 & (hdg==0 | abs(hdg)<abs(delta)));
pCorrect_total = sum(correct) / ntrials

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

% data.correct = (data.choice==2 & data.heading>0) | (data.choice==1 & data.heading<0) | ...
%     (rand<0.5 & (data.heading==0 | abs(data.heading)<abs(data.delta)));

if conftask==2
    data.PDW=data.conf;
end
subject = 'simul';


%% save it
% cd('/Users/chris/Documents/MATLAB')
% cd('/Users/stevenjerjian/Desktop/FetschLab/Analysis')
save(sprintf('2DAccSim_conftask%d_%dtrs.mat',conftask,ntrials),'data','cohs','deltas','hdgs','mods','origParams','RTtask','conftask','subject')



%% plots
if 1
    mods   = unique(data.modality);
    cohs   = unique(data.coherence);
    deltas = unique(data.delta);
    hdgs   = unique(data.heading);
    
    % means per condition, logistic fits
    parsedData = dots3DMP_parseData(data,mods,cohs,deltas,hdgs,conftask,RTtask);
    
    % gaussian fits
    gfit = dots3DMP_fit_cgauss(data,mods,cohs,deltas,conftask,RTtask);
    
    % plot it
    dots3DMP_plots(parsedData,mods,cohs,deltas,hdgs,conftask,RTtask)
%     dots3DMP_plots_cgauss_byCoh(gfit,parsedData,mods,cohs,deltas,hdgs,conftask,RTtask)
    
end

