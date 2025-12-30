function gfit = dots3DMP_fit_cgauss(data,mods,cohs,deltas,conftask,RTtask,D)
% calculate mean and standard deviations of gaussian fits to choice,
% confidence and RT data
%
% choice (PRight) is fit with a cumulative gaussian, returning mean (PSE)
% and sd (threshold)
% 
% confidence/PDW is fit with an inverted gaussian, with mean + sd params,
% as well as a baseline and amplitude. Baseline dictates an offset down
% from 1 (i.e. a base rate of low confidence)
% error function differs depending on conftask, because PDW is a binary
% outcome, whereas SEP is continuous
%
% RT is fit with a regular (right-side-up) gaussian, so baseline in this
% case is an offset up from 0
%
% D specifies what to do in the combined condition 
% (use just delta==0, or pool all deltas)
%
%
% Major updates:
%
% SJ 2021 at some point, fixed fits for PDW and RT
% SJ 07-2021 converted to fcn for cleaner workspace and variable handling
% SJ 03-2023 added fmincon option (currently set within code, should be a
% function argument)

if ~isfield(data,'oneTargConf'); data.oneTargConf=zeros(size(data.choice)); end

% define anonymous functions for fitting:


%% define anonymous functions for fitting:

% CHOICES - cumulative gaussian
cgauss = @(b,hdg) 1/2 * ( 1 + erf( (hdg-b(1))./(b(2)*sqrt(2)) ) );
% for probabilities, error is negative log likelihood of observing the data, which is
% [ log(Pright(hdg)) + log(1-(~Pright(hdg))) ]
cgauss_err = @(param,choice,hdg) -(sum(log(cgauss(param,hdg(choice))))+sum(log(1-cgauss(param,hdg(~choice)))));

% b(2) and b(3) are mu and sigma
% b(1) is 'amplitude' above baseline, b(4) is baseline level

% CONFIDENCE - 'flipped' Gaussian
% 'baseline' for flippedGauss is the highBet side, because of flip
if conftask==1 % continuous, sacc endpoint
    % for continuous values, error is sum squared error
    flippedGauss = @(b,hdg) 1 - ( b(1) .* exp(-(hdg-b(2)).^2 ./ (2*b(3).^2)) + b(4));
    flippedGauss_err = @(param,SEP,hdg) sum((flippedGauss(param,hdg)-SEP).^2);
elseif conftask==2 % PDW, probabilities
    
    % this was an attempt to force min/max 0-1 but can result in
    % non-invertible Hessians (typically for smaller datasets)
    % if you want to constrain parameters, then use fmincon
%     flippedGauss = @(b,hdg) 1 - ( min(max(b(1),0),1) .* exp(-(hdg-b(2)).^2 ./ (2*b(3).^2)) + min(max(b(4),0),1));
    flippedGauss = @(b,hdg) 1 - ( b(1) .* exp(-(hdg-b(2)).^2 ./ (2*b(3).^2)) + b(4));

    % error term minimizes negative log likelihood of observing PDW data
    % log prob of observing high bet on all trials where subj bet high, + log prob of low bet on trials where subj bet low
    % equivalent to cgauss error function but with PDW in place of choice
    flippedGauss_err = @(param,pdw,hdg) -( sum(log(flippedGauss(param,hdg(pdw)))) + sum(log(1-flippedGauss(param,hdg(~pdw)))) );
end

% RT - Gaussian, error is sum squared because RT is cont variable
gauss = @(b,hdg) b(1) .* exp(-(hdg-b(2)).^2 ./ (2*b(3).^2)) + b(4);
gauss_err = @(param,RT,hdg) sum((gauss(param,hdg)-RT).^2);

unc = 0; % saves biases from fminunc instead of fminsearch (SEs always are fminunc, and plots are always fminsearch)
useCon = 'none'; % 'all', 'none', 'delta'

%% set some fitting options + initial parameters

% parameter initial guesses
guess_cgauss = [0 2];
guess_fgauss = [0.05 0 2 0.5];
guess_gauss  = [0.5 0 2 1];

% bounds for fmincon
cgauss_LB = [-5 0];
cgauss_UB = [5 2.5];

% in fgauss case - bsln
% amp mu sigma bsln
fgauss_LB = [0 -5 0 0];
fgauss_UB = [1 5 4.5 1];

gauss_LB = [0 -5 0 0];
gauss_UB = [1 5 4.5 1];

fitOptions = optimset('display','final','MaxFunEvals',1e6,'MaxIter',1e6,'TolFun',1e-20);

%% first fit for all trials, irrespective of delta

if nargin<7 || isempty(D)
    D = length(deltas)+1; % (the extra column we made for pooling across deltas)
    % OR select just delta=0:
    % D = find(deltas==0);
end

% initialize vars for storing param fits
% deal func looks nicer, but is slow for some reason...
n = nan(length(mods),length(cohs),D);

muChoice = n; muChoiceSE = n;
sigmaChoice = n; sigmaChoiceSE = n;
fvalChoice = n;

amplConf = n; amplConfSE = n;
muConf = n; muConfSE = n;
sigmaConf = n; sigmaConfSE = n;
baselineConf = n; baselineConfSE = n;
fvalConf = n;

amplRT = n; amplRTse = n;
muRT = n; muRTse = n;
sigmaRT = n; sigmaRTse = n;
baselineRT = n; baselineRTse = n;
fvalRT = n;

for c = 1:length(cohs)
    % choice
    for m = 1:length(mods)     % m c d h
        if mods(m)==1
            I = data.modality==mods(m);
        else
            if D==length(deltas)+1
                I = data.modality==mods(m) & data.coherence==cohs(c); % all trials irrespective of delta
            else
                I = data.modality==mods(m) & data.coherence==cohs(c) & data.delta==deltas(D);
            end
        end
        I = I & ~data.oneTargConf;


        if strcmp(useCon,'all')
            [beta,fval] = fmincon(@(x) cgauss_err(x,data.choice(I)==2,data.heading(I)), guess_cgauss,[],[],[],[],...
                cgauss_LB,cgauss_UB,[],fitOptions);
        else
            [beta,fval] = fminsearch(@(x) cgauss_err(x,data.choice(I)==2,data.heading(I)), guess_cgauss, fitOptions);
        end

        [betaUnc,fvalunc,flag,~,~,hessian] = fminunc(@(x) cgauss_err(x,data.choice(I)==2,data.heading(I)), guess_cgauss,fitOptions);
        
        SE = sqrt(diag(inv(hessian)));
        muChoiceSE(m,c,D) = SE(1);
        sigmaChoiceSE(m,c,D) = SE(2);
        flagChoice(m,c,D) = flag;
        if unc
            muChoice(m,c,D) = betaUnc(1);
            sigmaChoice(m,c,D) = betaUnc(2);
            fvalChoice(m,c,D) = fvalunc;
        else
            muChoice(m,c,D) = beta(1);
            sigmaChoice(m,c,D) = beta(2);
            fvalChoice(m,c,D) = fval;
        end
    end
    
    % conf
    if conftask
        for m = 1:length(mods)
            if mods(m)==1
                I = data.modality==mods(m);
            else
                if D==length(deltas)+1
                    I = data.modality==mods(m) & data.coherence==cohs(c); % all trials irrespective of delta
                else
                    I = data.modality==mods(m) & data.coherence==cohs(c) & data.delta==deltas(D);
                end
            end
            I = I & ~data.oneTargConf;
           
            
%             fprintf('Fitting PDW for mod %d\n',mods(m))
            if conftask==1 % sacc endpoint
                Y = data.conf(I);
            elseif conftask==2 % PDW
                Y = data.PDW(I)==1;
            end

            if strcmp(useCon,'all')
                [beta,fval] = fmincon(@(x) flippedGauss_err(x,Y,data.heading(I)), guess_fgauss, [],[],[],[],...
                    fgauss_LB,fgauss_UB, [], fitOptions);
            else
                [beta,fval] = fminsearch(@(x) flippedGauss_err(x,Y,data.heading(I)), guess_fgauss, fitOptions);
            end
            [betaUnc,fvalunc,flag,~,~,hessian] = fminunc(@(x) flippedGauss_err(x,Y,data.heading(I)), guess_fgauss, fitOptions);


            SE = sqrt(diag(inv(hessian)));
            amplConfSE(m,c,D) = SE(1);
            muConfSE(m,c,D) = SE(2);
            sigmaConfSE(m,c,D) = SE(3);
            baselineConfSE(m,c,D) = SE(4);
            flagConf(m,c,D) = flag;

            if unc
                amplConf(m,c,D) = betaUnc(1);
                muConf(m,c,D) = betaUnc(2);
                sigmaConf(m,c,D) = betaUnc(3);
                baselineConf(m,c,D) = betaUnc(4);
                fvalConf(m,c,D) = fvalunc;
            else
                amplConf(m,c,D) = beta(1);
                muConf(m,c,D) = beta(2);
                sigmaConf(m,c,D) = beta(3);
                baselineConf(m,c,D) = beta(4);
                fvalConf(m,c,D) = fval;
            end
        end
    end
    
    % RT
    if RTtask
        for m = 1:length(mods)
            if mods(m)==1
                I = data.modality==mods(m);
            else
                if D==length(deltas)+1
                    I = data.modality==mods(m) & data.coherence==cohs(c); % all trials irrespective of delta
                else
                    I = data.modality==mods(m) & data.coherence==cohs(c) & data.delta==deltas(D);
                end
            end
            I = I & ~data.oneTargConf;

            if strcmp(useCon,'all')
                [beta,fval] = fmincon(@(x) gauss_err(x,data.RT(I),data.heading(I)), guess_gauss, [],[],[],[],...
                    gauss_LB,gauss_UB,[],fitOptions);
            else
                [beta,fval] = fminsearch(@(x) gauss_err(x,data.RT(I),data.heading(I)), guess_gauss,fitOptions);
            end

            [betaUnc,fvalunc,flag,~,~,hessian] = fminunc(@(x) gauss_err(x,data.RT(I),data.heading(I)), guess_gauss,fitOptions);
            SE = sqrt(diag(inv(hessian)));
            amplRTse(m,c,D) = SE(1);
            muRTse(m,c,D) = SE(2);
            sigmaRTse(m,c,D) = SE(3);
            baselineRTse(m,c,D) = SE(4);
            flagRT(m,c,D) = flag;
            if unc
                amplRT(m,c,D) = betaUnc(1);
                muRT(m,c,D) = betaUnc(2);
                sigmaRT(m,c,D) = betaUnc(3);
                baselineRT(m,c,D) = betaUnc(4);
                fvalRT(m,c,D) = fvalunc;
            else
                amplRT(m,c,D) = beta(1);
                muRT(m,c,D) = beta(2);
                sigmaRT(m,c,D) = beta(3);
                baselineRT(m,c,D) = beta(4);
                fvalRT(m,c,D) = fval;
            end
        end
    end
    
end


%% now separate by delta

for c = 1:length(cohs)
    % choice
    for d = 1:length(deltas)     % m c d h
        I = data.modality==3 & data.coherence==cohs(c) & data.delta==deltas(d);
         I = I & ~data.oneTargConf;

         if any(strcmp(useCon,{'all','delta'}))
             [beta,fval] = fmincon(@(x) cgauss_err(x,data.choice(I)==2,data.heading(I)), guess_cgauss,[],[],[],[],...
                 cgauss_LB,cgauss_UB,[],fitOptions);
         else
             [beta,fval] = fminsearch(@(x) cgauss_err(x,data.choice(I)==2,data.heading(I)), guess_cgauss, fitOptions);
         end

         [betaUnc,fvalunc,flag,~,~,hessian] = fminunc(@(x) cgauss_err(x,data.choice(I)==2,data.heading(I)), guess_cgauss,fitOptions);
         SE = sqrt(diag(inv(hessian)));
        muChoiceSE(3,c,d) = SE(1);
        sigmaChoiceSE(3,c,d) = SE(2);
        flagChoice(3,c,d) = flag;
        if unc
            muChoice(3,c,d) = betaUnc(1);
            sigmaChoice(3,c,d) = betaUnc(2);
            fvalChoice(3,c,d) = fvalunc;
        else
            muChoice(3,c,d) = beta(1);
            sigmaChoice(3,c,d) = beta(2);
            fvalChoice(3,c,d) = fval;
        end
    end
    
    % conf
    if conftask
        for d = 1:length(deltas)
            I = data.modality==3 & data.coherence==cohs(c) & data.delta==deltas(d);
            I = I & ~data.oneTargConf;

            if conftask==1 % sacc endpoint
                Y = data.conf(I);
            elseif conftask==2 % PDW
                Y = data.PDW(I)==1;
            end

            if any(strcmp(useCon,{'all','delta'}))
                [beta,fval] = fmincon(@(x) flippedGauss_err(x,Y,data.heading(I)), guess_fgauss, [],[],[],[],...
                    fgauss_LB,fgauss_UB, [], fitOptions);
            else
                [beta,fval] = fminsearch(@(x) flippedGauss_err(x,Y,data.heading(I)), guess_fgauss,fitOptions);
            end
            [betaUnc,fvalunc,flag,~,~,hessian] = fminunc(@(x) flippedGauss_err(x,Y,data.heading(I)), guess_fgauss,fitOptions);

            SE = sqrt(diag(inv(hessian)));
            amplConfSE(3,c,d) = SE(1);
            muConfSE(3,c,d) = SE(2);
            sigmaConfSE(3,c,d) = SE(3);
            baselineConfSE(3,c,d) = SE(4);
            flagConf(3,c,d) = flag;

            if unc
                amplConf(3,c,d) = betaUnc(1);
                muConf(3,c,d) = betaUnc(2);
                sigmaConf(3,c,d) = betaUnc(3);
                baselineConf(3,c,d) = betaUnc(4);
                fvalConf(3,c,d) = fvalunc;
            else
                amplConf(3,c,d) = beta(1);
                muConf(3,c,d) = beta(2);
                sigmaConf(3,c,d) = beta(3);
                baselineConf(3,c,d) = beta(4);
                fvalConf(3,c,d) = fval;
            end
        end
    end
    
    % RT
    if RTtask
        for d = 1:length(deltas)
            I = data.modality==3 & data.coherence==cohs(c) & data.delta==deltas(d);
             I = I & ~data.oneTargConf;
             
             if any(strcmp(useCon,{'all','delta'}))
                 [beta,fval] = fmincon(@(x) gauss_err(x,data.RT(I),data.heading(I)), guess_gauss, [],[],[],[],...
                     gauss_LB,gauss_UB,[],fitOptions);
             else
                 [beta,fval] = fminsearch(@(x) gauss_err(x,data.RT(I),data.heading(I)), guess_gauss,fitOptions);
            end
            [betaUnc,fvalunc,flag,~,~,hessian] = fminunc(@(x) gauss_err(x,data.RT(I),data.heading(I)), guess_gauss,fitOptions);

            SE = sqrt(diag(inv(hessian)));
            amplRTse(3,c,d) = SE(1);
            muRTse(3,c,d) = SE(2);
            sigmaRTse(3,c,d) = SE(3);
            baselineRTse(3,c,d) = SE(4);
            flagRT(3,c,d) = flag;
            if unc
                amplRT(3,c,d) = betaUnc(1);
                muRT(3,c,d) = betaUnc(2);
                sigmaRT(3,c,d) = betaUnc(3);
                baselineRT(3,c,d) = betaUnc(4);
                fvalRT(3,c,d) = fvalunc;
            else
                amplRT(3,c,d) = beta(1);
                muRT(3,c,d) = beta(2);
                sigmaRT(3,c,d) = beta(3);
                baselineRT(3,c,d) = beta(4);
                fvalRT(3,c,d) = fval;
            end
        end
    end
    
end

%% save outputs to struct, for clean output

gfit = struct();

gfit.choice.mu = muChoice;
gfit.choice.muSE = muChoiceSE;
gfit.choice.sigma = sigmaChoice;
gfit.choice.sigmaSE = sigmaChoiceSE;
gfit.choice.flag = flagChoice;

gfit.choice.func = cgauss;
gfit.choice.err = cgauss_err;
gfit.choice.guess = guess_cgauss;

gfit.choice.fval = fvalChoice;

if conftask
    gfit.conf.ampl = amplConf;
    gfit.conf.mu   = muConf;
    gfit.conf.sigma = sigmaConf;
    gfit.conf.bsln = baselineConf;
    gfit.conf.flag = flagConf;

    gfit.conf.amplSE = amplConfSE;
    gfit.conf.muSE  = muConfSE;
    gfit.conf.sigmaSE = sigmaConfSE;
    gfit.conf.bslnSE = baselineConfSE;

    gfit.conf.func = flippedGauss;
    gfit.conf.err = flippedGauss_err;
    gfit.conf.guess = guess_fgauss;

    gfit.conf.fval = fvalConf;
end

if RTtask
    gfit.RT.ampl = amplRT;
    gfit.RT.mu   = muRT;
    gfit.RT.sigma = sigmaRT;
    gfit.RT.bsln = baselineRT;
    gfit.RT.flag = flagRT;

    gfit.RT.amplSE = amplRTse;
    gfit.RT.muSE  = muRTse;
    gfit.RT.sigmaSE = sigmaRTse;
    gfit.RT.bslnSE = baselineRTse;

    gfit.RT.func = gauss;
    gfit.RT.err  = gauss_err;
    gfit.RT.guess = guess_gauss;
    
    gfit.RT.fval = fvalRT;
end

gfit.D = D;
