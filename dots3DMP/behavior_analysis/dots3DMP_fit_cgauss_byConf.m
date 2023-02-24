function gfit = dots3DMP_fit_cgauss_byConf(data,mods,cohs,deltas,conftask,RTtask,D)

% same as dots3DMP_fit_cgauss, but for splits by high and low confidence
% SJ 07-2021 converted to function for cleanear workspace

% define anonymous functions for fitting:

cgauss = @(b,hdg) 1/2 * ( 1 + erf( (hdg-b(1))./(b(2)*sqrt(2)) ) );
    % for probabilities, error is negative log likelihood of observing the data, which is
    % [ log(Pright(hdg)) + log(1-(~Pright(hdg))) ]
cgauss_err = @(param,choice,hdg) -(sum(log(cgauss(param,hdg(choice))))+sum(log(1-cgauss(param,hdg(~choice))))); 

% CONFIDENCE
% skip, obviously no sense fitting here

% RT
gauss = @(b,hdg) b(1) .* exp(-(hdg-b(2)).^2 ./ (2*b(3).^2)) + b(4);
gauss_err = @(param,SEP,hdg) sum((gauss(param,hdg)-SEP).^2);

unc = 0; % saves biases from fminunc instead of fminsearch (SEs always are fminunc, and plots are always fminsearch)


% parameter initial guesses
guess_cgauss = [0 3];
guess_gauss  = [1 0 6 1];

% split confidence into high and low 
% easy for PDW
% for sacc endpoint, for now just split by median
% eventually...split for each subject separately?
% or just make sure conf is normalized within subject before this stage
if conftask==1
    hiConf = data.conf >= median(data.conf);
elseif conftask==2 
    hiConf = data.PDW;
end


%% first, for all trials irrespective of delta

% SJ add D as option, 11/21/2022
if isempty(D) || nargin < 7
    D = length(deltas)+1; % (the extra column we made for pooling across deltas)
    % OR select just delta=0:
%     D = find(deltas==0);
end

% initialize vars for storing param fits
% deal func looks nicer, but is slow for some reason...
n = nan(length(mods)+1,length(cohs),length(deltas),2);
                                               % ^2 for high and low confidence fits

muPMF = n; muPMFse = n;
sigmaPMF = n; sigmaPMFse = n;

amplRT = n; amplRTse = n;
muRT = n; muRTse = n;
sigmaRT = n; sigmaRTse = n;
baselineRT = n; baselineRTse = n;

for c = 1:length(cohs)
    % choice
    for m = 1:length(mods)+1     % m c d h
        
        if m==length(mods)+1
            I = true(size(data.modality));
        elseif D==length(deltas)+1
            I = data.modality==mods(m) & data.coherence==cohs(c); % all trials irrespective of delta
        else
            I = data.modality==mods(m) & data.coherence==cohs(c) & data.delta==deltas(D);
        end
        
        if sum(I)==0, continue, end
        
        I_hi = I & hiConf & data.oneTargConf==0;
        I_lo = I & ~hiConf & data.oneTargConf==0;
        
        % for high conf
        beta = fminsearch(@(x) cgauss_err(x,data.choice(I_hi)==2,data.heading(I_hi)), [0 3]);
        [betaUnc,~,~,~,~,hessian] = fminunc(@(x) cgauss_err(x,data.choice(I_hi)==2,data.heading(I_hi)), [0 3]);
        SE = sqrt(diag(inv(hessian)));
        muPMFse(m,c,D,1) = SE(1);
        sigmaPMFse(m,c,D,1) = SE(2);
        if unc
            muPMF(m,c,D,1) = betaUnc(1);
            sigmaPMF(m,c,D,1) = betaUnc(2);
        else
            muPMF(m,c,D,1) = beta(1);
            sigmaPMF(m,c,D,1) = beta(2);
        end
        
        % repeat for low conf
        beta = fminsearch(@(x) cgauss_err(x,data.choice(I_lo)==2,data.heading(I_lo)), guess_cgauss);
        [betaUnc,~,~,~,~,hessian] = fminunc(@(x) cgauss_err(x,data.choice(I_lo)==2,data.heading(I_lo)), guess_cgauss);
        SE = sqrt(diag(inv(hessian)));
        muPMFse(m,c,D,2) = SE(1);
        sigmaPMFse(m,c,D,2) = SE(2);
        if unc
            muPMF(m,c,D,2) = betaUnc(1);
            sigmaPMF(m,c,D,2) = betaUnc(2);
        else
            muPMF(m,c,D,2) = beta(1);
            sigmaPMF(m,c,D,2) = beta(2);
        end
    end
    
    % RT
    if RTtask
        for m = 1:length(mods)+1      
            if m==length(mods)+1
                I = true(size(data.modality));
            elseif D==length(deltas)+1
                I = data.modality==mods(m) & data.coherence==cohs(c); % all trials irrespective of delta
            else
                I = data.modality==mods(m) & data.coherence==cohs(c) & data.delta==deltas(D);
            end
            if sum(I)==0, continue, end
            
            I_hi = I & hiConf & data.oneTargConf==0;
            I_lo = I & ~hiConf & data.oneTargConf==0;
            
            % high conf
            beta = fminsearch(@(x) gauss_err(x,data.RT(I_hi),data.heading(I_hi)), guess_gauss);
            [betaUnc,~,~,~,~,hessian] = fminunc(@(x) gauss_err(x,data.RT(I_hi),data.heading(I_hi)), guess_gauss);
            SE = sqrt(diag(inv(hessian)));
            amplRTse(m,c,D,1) = SE(1);
            muRTse(m,c,D,1) = SE(2);
            sigmaRTse(m,c,D,1) = SE(3);
            baselineRTse(m,c,D,1) = SE(4);
            if unc
                amplRT(m,c,D,1) = betaUnc(1);
                muRT(m,c,D,1) = betaUnc(2);
                sigmaRT(m,c,D,1) = betaUnc(3);
                baselineRT(m,c,D,1) = betaUnc(4);
            else
                amplRT(m,c,D,1) = beta(1);
                muRT(m,c,D,1) = beta(2);
                sigmaRT(m,c,D,1) = beta(3);
                baselineRT(m,c,D,1) = beta(4);
            end
            
            % repeat for low conf
            beta = fminsearch(@(x) gauss_err(x,data.RT(I_lo),data.heading(I_lo)), guess_gauss);
            [betaUnc,~,~,~,~,hessian] = fminunc(@(x) gauss_err(x,data.RT(I_lo),data.heading(I_lo)), guess_gauss);
            SE = sqrt(diag(inv(hessian)));
            amplRTse(m,c,D,2) = SE(1);
            muRTse(m,c,D,2) = SE(2);
            sigmaRTse(m,c,D,2) = SE(3);
            baselineRTse(m,c,D,2) = SE(4);
            if unc
                amplRT(m,c,D,2) = betaUnc(1);
                muRT(m,c,D,2) = betaUnc(2);
                sigmaRT(m,c,D,2) = betaUnc(3);
                baselineRT(m,c,D,2) = betaUnc(4);
            else
                amplRT(m,c,D,2) = beta(1);
                muRT(m,c,D,2) = beta(2);
                sigmaRT(m,c,D,2) = beta(3);
                baselineRT(m,c,D,2) = beta(4);
            end
        end
    end
    
end


%% now separate by delta

% this needs work to fix, if we want to bother showing it...

%{
for c = 1:length(cohs)
    % choice
    for d = 1:length(deltas)     % m c d h
        I = data.modality==3 & data.coherence==cohs(c) & data.delta==deltas(d);
        I_hi = I & hiConf;
        I_lo = I & ~hiConf;
            
        % high conf
        beta = fminsearch(@(x) cgauss_err(x,data.choice(I_hi)==2,data.heading(I_hi)), guess_cgauss);
        [betaUnc,~,~,~,~,hessian] = fminunc(@(x) cgauss_err(x,data.choice(I_hi)==2,data.heading(I_hi)), guess_cgauss);
        SE = sqrt(diag(inv(hessian)));
        muPMFse(3,c,d,1) = SE(1);
        sigmaPMFse(3,c,d,1) = SE(2);        
        if unc
            muPMF(3,c,d,1) = betaUnc(1);
            sigmaPMF(3,c,d,1) = betaUnc(2);
        else
            muPMF(3,c,d,1) = beta(1);
            sigmaPMF(3,c,d,1) = beta(2);
        end
        
        % repeat for low conf
        beta = fminsearch(@(x) cgauss_err(x,data.choice(I_lo)==2,data.heading(I_lo)), [0 3]);
        [betaUnc,~,~,~,~,hessian] = fminunc(@(x) cgauss_err(x,data.choice(I_lo)==2,data.heading(I_lo)), [0 3]);
        SE = sqrt(diag(inv(hessian)));
        muPMFse(3,c,d,2) = SE(1);
        sigmaPMFse(3,c,d,2) = SE(2);        
        if unc
            muPMF(3,c,d,2) = betaUnc(1);
            sigmaPMF(3,c,d,2) = betaUnc(2);
        else
            muPMF(3,c,d,2) = beta(1);
            sigmaPMF(3,c,d,2) = beta(2);
        end
    end
    
    % RT
    if RTtask
        for d = 1:length(deltas)
            I = data.modality==3 & data.coherence==cohs(c) & data.delta==deltas(d);
            I_hi = I & hiConf;
            I_lo = I & ~hiConf;
            
            % high conf
            beta = fminsearch(@(x) gauss_err(x,data.RT(I_hi),data.heading(I_hi)), guess_gauss);
            [betaUnc,~,~,~,~,hessian] = fminunc(@(x) gauss_err(x,data.RT(I_hi),data.heading(I_hi)), guess_gauss);
            SE = sqrt(diag(inv(hessian)));
            amplRTse(3,c,d,1) = SE(1);
            muRTse(3,c,d,1) = SE(2);
            sigmaRTse(3,c,d,1) = SE(3);
            baselineRTse(3,c,d,1) = SE(4);
           if unc
                amplRT(3,c,d,1) = betaUnc(1);
                muRT(3,c,d,1) = betaUnc(2);
                sigmaRT(3,c,d,1) = betaUnc(3);
                baselineRT(3,c,d,1) = betaUnc(4);
            else
                amplRT(3,c,d,1) = beta(1);
                muRT(3,c,d,1) = beta(2);
                sigmaRT(3,c,d,1) = beta(3);
                baselineRT(3,c,d,1) = beta(4);
           end
            
           % low conf
            beta = fminsearch(@(x) gauss_err(x,data.RT(I_lo),data.heading(I_lo)), guess_gauss);
            [betaUnc,~,~,~,~,hessian] = fminunc(@(x) gauss_err(x,data.RT(I_lo),data.heading(I_lo)), guess_gauss);
            SE = sqrt(diag(inv(hessian)));
            amplRTse(3,c,d,2) = SE(1);
            muRTse(3,c,d,2) = SE(2);
            sigmaRTse(3,c,d,2) = SE(3);
            baselineRTse(3,c,d,2) = SE(4);
           if unc
                amplRT(3,c,d,2) = betaUnc(1);
                muRT(3,c,d,2) = betaUnc(2);
                sigmaRT(3,c,d,2) = betaUnc(3);
                baselineRT(3,c,d,2) = betaUnc(4);
            else
                amplRT(3,c,d,2) = beta(1);
                muRT(3,c,d,2) = beta(2);
                sigmaRT(3,c,d,2) = beta(3);
                baselineRT(3,c,d,2) = beta(4);
            end
        end
    end
    
end
%}

% save outputs
gfit = struct();

gfit.choice.mu = muPMF;
gfit.choice.muSE = muPMFse;
gfit.choice.sigma = sigmaPMF;
gfit.choice.sigmaSE = sigmaPMFse;

gfit.choice.func = cgauss;
gfit.choice.err = cgauss_err;
gfit.choice.guess = guess_cgauss;

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

end