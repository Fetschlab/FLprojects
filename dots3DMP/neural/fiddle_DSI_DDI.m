% test some tuning indices

% CF 2024


%% simulate a tuning curve

clear

% anonymous function for von Mises (circular Gaussian)
tuning_vonMises = @(b,dir) b(1) * exp(b(2)*cosd(dir-b(3))) / (2*pi*besseli(0,b(2))) + b(4);

dirAxis = -180:45:180;

ampl = 30;
conc = 2; % concentration, aka inverse width, analogous to 1/sigma^2
peak = 0;
baseline = 10;
b = [ampl conc peak baseline];
    
frMean = tuning_vonMises(b,dirAxis);
figure;plot(dirAxis,frMean);


%% simulate an experiment

ntrials = 100;
dirs = randsample(dirAxis,ntrials,true)'; % stimulus direction on each trial
R = nan(size(dirs)); % response (spike rate) on each trial
for n = 1:ntrials
    
    R(n) = poissrnd(tuning_vonMises(b,dirs(n))); % poisson noise
    % OR
%     R(n) = tuning_vonMises(b,dirs(n)) + normrnd(0,sqrt(tuning_vonMises(b,dirs(n)))); % Gaussian noise with low variance
    % OR    
%     R(n) = tuning_vonMises(b,dirs(n)) + normrnd(0,8*sqrt(tuning_vonMises(b,dirs(n)))); % Gaussian noise with high variance, notice how DDI changes

    if R(n)<0; R(n)=0; end
end

% plot tuning
rMean = nan(size(dirAxis));
rSE = nan(size(dirAxis)); % standard error
rVar = nan(size(dirAxis)); % variance, for later
for n = 1:length(dirAxis)
    i = dirs==dirAxis(n); % select all trials for a given dir
    rMean(n) = nanmean(R(i));
    rVar(n) = nanstd(R(i))^2;
    rSE(n) = nanstd(R(i))/sqrt(nansum(i));
end
figure; errorbar(dirAxis, rMean, rSE);


%% compare indices

% a typical direction selectivity index
DSI = (max(rMean)-min(rMean)) / (max(rMean)+min(rMean)) % or something like that


% Ldir = % add equation for normalized vector sum

% direction discrimination index ; Nguyenkim & DeAngelis 2003,
% after Prince et al. 2002; measures ability of ideal observer to
% discriminate between min and max resp dirs ('coarse' discrim)
DDI = calcDDI(dirAxis,dirs,rMean,R)


% Fisher information: measures ability of ideal observer to discriminate
% between nearby directions ('fine' discrimi); can be estimated as
% the squared slope of the tuning curve at a reference value, divided by
% the variance at that value:

% First interpolate:
dirAxInterp = -180:180;
rMeanSmoothed = interp1(dirAxis,rMean,dirAxInterp,'spline');
rVarSmoothed  = interp1(dirAxis,rVar,dirAxInterp,'spline');

% get max slope (absolute value), and estimated variance at that same value
absSlope = abs(gradient(rMeanSmoothed));
FIslope = max(absSlope);
FIvar = rVarSmoothed(absSlope==max(absSlope));

% Fisher Info
FI = FIslope^2 / FIvar
% units are kinda arbitrary; if you did this for a bunch of neurons you
% could just normalize the values to a range of 0-1 so it is comparable to
% the other indices





%%

function DDI = calcDDI(dirAxis,dirs,rMean,R)
% direction discrimination index, Nguyenkim & DeAngelis 2003, after Prince
% et al. 2002; measures ability of ideal observer to discriminate between
% min and max resp dirs, i.e. takes into account response variability

    Rmax = max(rMean);
    Rmin = min(rMean);

    % sum squared error around the mean for each trial
    err = nan(size(dirs));
    for n = 1:length(dirs)
        thisMean = rMean(dirAxis==dirs(n));
        err(n) = R(n)-thisMean;
    end
    SSE = nansum(err.^2);

    N = length(R);
    M = length(dirAxis);

    DDI = (Rmax - Rmin) / (Rmax - Rmin + 2*sqrt(SSE/(N-M)));

end

