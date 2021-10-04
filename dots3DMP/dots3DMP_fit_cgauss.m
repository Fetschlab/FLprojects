function gfit = dots3DMP_fit_cgauss(data,mods,cohs,deltas,conftask,RTtask)
% SJ 07-2021 converted to function for cleaner workspace

% define anonymous functions for fitting:

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
    flippedGauss_err = @(param,SEP,hdg) nansum((flippedGauss(param,hdg)-SEP).^2);
elseif conftask==2 % PDW, probabilities
    
    % dont force min/max 0,1... resulting in non-invertible Hessians
%     flippedGauss = @(b,hdg) 1 - ( min(max(b(1),0),1) .* exp(-(hdg-b(2)).^2 ./ (2*b(3).^2)) + min(max(b(4),0),1));
    flippedGauss = @(b,hdg) 1 - ( b(1) .* exp(-(hdg-b(2)).^2 ./ (2*b(3).^2)) + b(4));

    % negative log likelihood of observing PDW data
    % log prob of observing high bet on all trials where subj bet
    % high, + log prob of low bet on trials where subj bet low
    % equivalent to cgauss error function
    flippedGauss_err = @(param,pdw,hdg) -( sum(log(flippedGauss(param,hdg(pdw)))) + sum(log(1-flippedGauss(param,hdg(~pdw)))) );
end

% RT - Gaussian, error is sum squared because RT is cont variable
gauss = @(b,hdg) b(1) .* exp(-(hdg-b(2)).^2 ./ (2*b(3).^2)) + b(4);
gauss_err = @(param,RT,hdg) sum((gauss(param,hdg)-RT).^2);

unc = 0; % saves biases from fminunc instead of fminsearch (SEs always are fminunc, and plots are always fminsearch)

% parameter initial guesses
guess_cgauss = [0 3];
guess_fgauss = [0.2 0 6 0.5];
guess_gauss  = [0.5 0 6 1];

fitOptions = optimset('display','none','MaxFunEvals',1000,'MaxIter',1000);

%% first, for all trials irrespective of delta
% D = length(deltas)+1; % (the extra column we made for pooling across deltas)
% OR select just delta=0:
D = find(deltas==0);

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
        if m==1
            I = data.modality==mods(m);
        else
            if D==length(deltas)+1
                I = data.modality==mods(m) & data.coherence==cohs(c); % all trials irrespective of delta
            else
                I = data.modality==mods(m) & data.coherence==cohs(c) & data.delta==deltas(D);
            end
        end
        
        [beta,fval] = fminsearch(@(x) cgauss_err(x,data.choice(I)==2,data.heading(I)), guess_cgauss,fitOptions);
        [betaUnc,fvalunc,flag,~,~,hessian] = fminunc(@(x) cgauss_err(x,data.choice(I)==2,data.heading(I)), guess_cgauss,fitOptions);
        SE = sqrt(diag(inv(hessian)));
        muChoiceSE(m,c,D) = SE(1);
        sigmaChoiceSE(m,c,D) = SE(2);
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
            if m==1
                I = data.modality==mods(m);
            else
                if D==length(deltas)+1
                    I = data.modality==mods(m) & data.coherence==cohs(c); % all trials irrespective of delta
                else
                    I = data.modality==mods(m) & data.coherence==cohs(c) & data.delta==deltas(D);
                end
            end
            
            if conftask==1 % sacc endpoint
                [beta,fval] = fminsearch(@(x) flippedGauss_err(x,data.conf(I),data.heading(I)), guess_fgauss,fitOptions);
                [betaUnc,fvalunc,flag,~,~,hessian] = fminunc(@(x) flippedGauss_err(x,data.conf(I),data.heading(I)), guess_fgauss,fitOptions);
            elseif conftask==2 % PDW
                [beta,fval] = fminsearch(@(x) flippedGauss_err(x,data.PDW(I)==1,data.heading(I)), guess_fgauss,fitOptions);
                [betaUnc,fvalunc,flag,~,~,hessian] = fminunc(@(x) flippedGauss_err(x,data.PDW(I)==1,data.heading(I)), guess_fgauss,fitOptions);
            end
            
            SE = sqrt(diag(inv(hessian)));
            amplConfSE(m,c,D) = SE(1);
            muConfSE(m,c,D) = SE(2);
            sigmaConfSE(m,c,D) = SE(3);
            baselineConfSE(m,c,D) = SE(4);
            fvalConf(m,c,D) = fval;
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
            if m==1
                I = data.modality==mods(m);
            else
                if D==length(deltas)+1
                    I = data.modality==mods(m) & data.coherence==cohs(c); % all trials irrespective of delta
                else
                    I = data.modality==mods(m) & data.coherence==cohs(c) & data.delta==deltas(D);
                end
            end
            [beta,fval] = fminsearch(@(x) gauss_err(x,data.RT(I),data.heading(I)), guess_gauss,fitOptions);
            [betaUnc,fvalunc,flag,~,~,hessian] = fminunc(@(x) gauss_err(x,data.RT(I),data.heading(I)), guess_gauss,fitOptions);
            SE = sqrt(diag(inv(hessian)));
            amplRTse(m,c,D) = SE(1);
            muRTse(m,c,D) = SE(2);
            sigmaRTse(m,c,D) = SE(3);
            baselineRTse(m,c,D) = SE(4);
            fvalRT(m,c,D) = fval;
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
        
        [beta,fval] = fminsearch(@(x) cgauss_err(x,data.choice(I)==2,data.heading(I)), guess_cgauss,fitOptions);
        [betaUnc,fvalunc,flag,~,~,hessian] = fminunc(@(x) cgauss_err(x,data.choice(I)==2,data.heading(I)), guess_cgauss,fitOptions);
        SE = sqrt(diag(inv(hessian)));
        muChoiceSE(3,c,d) = SE(1);
        sigmaChoiceSE(3,c,d) = SE(2);
        fvalChoice(3,c,d) = fval;
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
            
            if conftask==1 % sacc endpoint
                [beta,fval] = fminsearch(@(x) flippedGauss_err(x,data.conf(I),data.heading(I)), guess_fgauss,fitOptions);
                [betaUnc,fvalunc,flag,~,~,hessian] = fminunc(@(x) flippedGauss_err(x,data.conf(I),data.heading(I)), guess_fgauss,fitOptions);
            elseif conftask==2 % PDW
                [beta,fval] = fminsearch(@(x) flippedGauss_err(x,data.PDW(I)==1,data.heading(I)), guess_fgauss,fitOptions);
                [betaUnc,fvalunc,flag,~,~,hessian] = fminunc(@(x) flippedGauss_err(x,data.PDW(I)==1,data.heading(I)), guess_fgauss,fitOptions);
            end
            
            SE = sqrt(diag(inv(hessian)));
            amplConfSE(3,c,d) = SE(1);
            muConfSE(3,c,d) = SE(2);
            sigmaConfSE(3,c,d) = SE(3);
            baselineConfSE(3,c,d) = SE(4);
            fvalConf(3,c,d) = fval;

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
            [beta,fval] = fminsearch(@(x) gauss_err(x,data.RT(I),data.heading(I)), guess_gauss,fitOptions);
            [betaUnc,fvalunc,flag,~,~,hessian] = fminunc(@(x) gauss_err(x,data.RT(I),data.heading(I)), guess_gauss,fitOptions);
            SE = sqrt(diag(inv(hessian)));
            amplRTse(3,c,d) = SE(1);
            muRTse(3,c,d) = SE(2);
            sigmaRTse(3,c,d) = SE(3);
            baselineRTse(3,c,d) = SE(4);
            fvalRT(3,c,d) = fval;
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

% save outputs
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
