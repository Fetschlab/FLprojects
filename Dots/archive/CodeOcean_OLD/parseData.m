function parsedData = parseData(data,conftask,RTtask,RTCorrOnly)


if nargin < 4, RTCorrOnly = 0; end
if nargin < 3, RTtask = 0; end
if nargin < 2, conftask = 0; end

if RTtask==0
    data.RT = nan(size(data.choice));
end
switch conftask
    case 0
        data.wager = nan(size(data.choice));
        data.conf = nan(size(data.choice));
    case 1
        data.wager = nan(size(data.choice));
    case 2
        data.conf= nan(size(data.choice));
end

% before alpha offset, for choice/RT splits 
% [this is weird; see note in simDDM_postProcessing]
if ~isfield(data,'PDW_preAlpha') % TEMP: try this
    data.wager_preAlpha = data.wager;
end

cohs = unique(data.cohs);
if cohs(end)==1; cohs(end)=[]; end % one or more files had a few stray 100% coh trials in there by mistake

n = nan(length(cohs),1);

n_all = n; % all
n_high = n; % high
n_low = n; % low
pRight = n;
pCorrect = n;
pRightHigh = n;
pRightLow = n;

nRT_all = n; % all
nRT_high = n; % high
nRT_low = n; % low
nRT_corr = n; % correct
nRT_err = n; % error
RTmean = n;
RTse = n;
RTmeanHigh = n;
RTseHigh = n;
RTmeanLow = n;
RTseLow = n;
RTmeanCorr = n;
RTseCorr = n;
RTmeanErr = n;
RTseErr = n;


nPDW_all = n; % all
nPDW_corr = n; % corr
nPDW_err = n; % err
pHigh = n;
pHighCorr = n;
pHighErr = n;

nConf_all = n; % all
nConf_corr = n; % correct
nConf_err = n; % error
confMean = n;
confSE = n;
confMeanCorr = n;
confSEcorr = n;
confMeanErr = n;
confSEerr = n;


% for logistic regression (or whatever)
xVals = linspace(cohs(1),cohs(end),100);

K = ~isnan(data.RT); % need to start tagging non-RT trials as NaN!
if RTCorrOnly; K = K & data.correct==1; end

for c = 1:length(cohs)
    % choice
    J = data.cohs==cohs(c);
    n_all(c) = sum(J);
    pRight(c) = sum(J & data.choice==1) / n_all(c); % 0 is left, 1 is right

    pCorrect(c) = sum(J & data.correct==1) / n_all(c); 
    
    if conftask==1
        Jhi = data.cohs==cohs(c) & data.conf>=median(data.conf);
        Jlo = data.cohs==cohs(c) & data.conf<median(data.conf);
    elseif conftask==2
        Jhi = data.cohs==cohs(c) & data.wager_preAlpha==1;
        Jlo = data.cohs==cohs(c) & data.wager_preAlpha==0;
    end
    if conftask
        n_high(c) = sum(Jhi);
        pRightHigh(c) = sum(Jhi & data.choice==1) / n_high(c); 
        n_low(c) = sum(Jlo);
        pRightLow(c) = sum(Jlo & data.choice==1) / n_low(c);        
    end
    
    % RT
    nRT_all(c) = sum(J & K);
    RTmean(c) = mean(data.RT(J & K));
    RTse(c) = std(data.RT(J & K))/sqrt(nRT_all(c));
    
    if conftask
        nRT_high(c) = sum(Jhi & K);
        RTmeanHigh(c) = mean(data.RT(Jhi & K));
        RTseHigh(c) = std(data.RT(Jhi & K))/sqrt(nRT_high(c));

        nRT_low(c) = sum(Jlo & K);
        RTmeanLow(c) = mean(data.RT(Jlo & K));
        RTseLow(c) = std(data.RT(Jlo & K))/sqrt(nRT_low(c));
    end
    
    nRT_corr(c) = sum(J & K & data.correct==1);
    RTmeanCorr(c) = mean(data.RT(J & K & data.correct==1));
    RTseCorr(c) = std(data.RT(J & K & data.correct==1))/sqrt(nRT_corr(c));
    
    nRT_err(c) = sum(J & K & data.correct==0);
    RTmeanErr(c) = mean(data.RT(J & K & data.correct==0));
    RTseErr(c) = std(data.RT(J & K & data.correct==0))/sqrt(nRT_err(c));
    
    % pdw
    if conftask==2
        L = ~isnan(data.wager);
        nPDW_all(c) = sum(J & L);
        pHigh(c) = sum(J & L & data.wager==1) / nPDW_all(c); % 1 is high-bet

        Lcor = ~isnan(data.wager) & data.correct==1;
        nPDW_corr(c) = sum(J & Lcor);
        pHighCorr(c) = sum(J & Lcor & data.wager==1) / nPDW_corr(c);

        Lerr = ~isnan(data.wager) & data.correct==0;
        nPDW_err(c) = sum(J & Lerr);
        if nPDW_err(c)>15
            pHighErr(c) = sum(J & Lerr & data.wager==1) / nPDW_err(c);
        else % ignore small N (e.g. high coh low bet)
            pHighErr(c) = NaN;
        end
    end    
    
    % conf
    if conftask==1
        M = ~isnan(data.conf);
        nConf_all(c) = sum(J & M);
        confMean(c) = mean(data.conf(J & M));
        confSE(c) = std(data.conf(J & M))/sqrt(nConf_all(c));

        MM = ~isnan(data.conf) & data.correct==1;
        nConf_corr(c) = sum(J & MM);
        confMeanCorr(c) = mean(data.conf(J & MM));
        confSEcorr(c) = std(data.conf(J & MM))/sqrt(nConf_corr(c));

        MMM = ~isnan(data.conf) & data.correct==0;
        nConf_err(c) = sum(J & MMM);
        confMeanErr(c) = mean(data.conf(J & MMM));
        confSEerr(c) = std(data.conf(J & MMM))/sqrt(nConf_err(c));
    end    
end

pRightSE = sqrt( (pRight.*(1-pRight)) ./ n_all );
pCorrectSE = sqrt( (pCorrect.*(1-pCorrect)) ./ n_all );
pRightSEhigh = sqrt( (pRightHigh.*(1-pRightHigh)) ./ n_high );
pRightSElow = sqrt( (pRightLow.*(1-pRightLow)) ./ n_low );

pHighSE = sqrt( (pHigh.*(1-pHigh)) ./ nPDW_all );
pHighSEcorr = sqrt( (pHighCorr.*(1-pHighCorr)) ./ nPDW_corr );
pHighSEerr = sqrt( (pHighErr.*(1-pHighErr)) ./ nPDW_err );


% fit logistic regression

    % all trials
X = data.cohs;
y = data.choice==1; % 1 is right
[B1, ~, stats1] = glmfit(X, y, 'binomial');
yVals1 = glmval(B1,xVals,'logit');

    % high conf/bet
if conftask==1
    I = data.conf>=median(data.conf);
elseif conftask==2
    I = data.wager_preAlpha==1;
else
    I = false(size(data.wager));
end
X = data.cohs(I);
y = data.choice(I)==1;
[B2, ~, stats2] = glmfit(X, y, 'binomial');
yVals2 = glmval(B2,xVals,'logit');

    % low conf/bet
if conftask==1
    I = data.conf<median(data.conf);
elseif conftask==2
    I = data.wager_preAlpha==0;
else
    I = false(size(data.wager));
end
X = data.cohs(I);
y = data.choice(I)==1;
[B3, ~, stats3] = glmfit(X, y, 'binomial');
yVals3 = glmval(B3,xVals,'logit');

    % all trials with indicator term for conf
X = [data.cohs data.wager.*data.cohs];
y = data.choice==1; % 1 is right
[B4, ~, stats4] = glmfit(X, y, 'binomial');
   

parsedData = struct();

parsedData.n_all = n_all;
parsedData.n_high = n_high;
parsedData.n_low = n_low;
parsedData.nRT_all = nRT_all;
parsedData.nRT_high = nRT_high;
parsedData.nRT_low = nRT_low;
parsedData.nPDW_all = nPDW_all;
parsedData.nPDW_corr = nPDW_corr;
parsedData.nPDW_err = nPDW_err;
parsedData.nConf_all = nConf_all;
parsedData.nConf_corr = nConf_corr;
parsedData.nConf_err = nConf_err;

parsedData.pRight = pRight;
parsedData.pRightHigh = pRightHigh;
parsedData.pRightLow = pRightLow;
parsedData.pRightSE = pRightSE;
parsedData.pRightSEhigh = pRightSEhigh;
parsedData.pRightSElow = pRightSElow;
parsedData.pCorrect = pCorrect;
parsedData.pCorrectSE = pCorrectSE;
parsedData.xVals = xVals;
parsedData.yVals1 = yVals1;
parsedData.yVals2 = yVals2;
parsedData.yVals3 = yVals3;
parsedData.B1 = B1;
parsedData.B2 = B2;
parsedData.B3 = B3;
parsedData.B4 = B4;
parsedData.stats1 = stats1;
parsedData.stats2 = stats2;
parsedData.stats3 = stats3;
parsedData.stats4 = stats4;
parsedData.yVals1 = yVals1;
parsedData.yVals2 = yVals2;
parsedData.yVals3 = yVals3;

if conftask==1
    parsedData.confMean = confMean;
    parsedData.confSE = confSE;
    parsedData.confMeanCorr = confMeanCorr;
    parsedData.confSEcorr = confSEcorr;
    parsedData.confMeanErr = confMeanErr;
    parsedData.confSEerr = confSEerr;
elseif conftask==2
    parsedData.pHigh = pHigh;
    parsedData.pHighSE = pHighSE;
    parsedData.pHighCorr = pHighCorr;
    parsedData.pHighSEcorr = pHighSEcorr;
    parsedData.pHighErr = pHighErr;
    parsedData.pHighSEerr = pHighSEerr;
end

if RTtask
    parsedData.RTmean = RTmean;
    parsedData.RTse = RTse;
    parsedData.RTmeanHigh = RTmeanHigh;
    parsedData.RTseHigh = RTseHigh;
    parsedData.RTmeanLow = RTmeanLow;
    parsedData.RTseLow = RTseLow;
    parsedData.RTmeanCorr = RTmeanCorr;
    parsedData.RTseCorr = RTseCorr;
    parsedData.RTmeanErr = RTmeanErr;
    parsedData.RTseErr = RTseErr;
end





