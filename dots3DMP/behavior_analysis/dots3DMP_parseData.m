function parsedData = dots3DMP_parseData(data,mods,cohs,deltas,hdgs,conftask,RTtask,useAbsHdg,RTCorrOnly)
% SJ 07-2021 converted to function, for cleaner workspace

if nargin < 9, RTCorrOnly = 0; end
if nargin < 8, useAbsHdg = 0; end
if nargin < 7, RTtask = 0; end
if nargin < 6, conftask = 0; end

%% parse data
% create and use matrices of summary data indexed by variables of interest

% switch back to 1:2 if 0:1 (e.g. from fitting code)
if max(data.choice)==1
    data.choice = data.choice+1;
end

if useAbsHdg
    hdgs = unique(abs(hdgs));
    data.heading = abs(data.heading);
end
if RTCorrOnly
    use_trs = data.correct | abs(data.heading) < 1e-6;
end

n = nan(length(mods),length(cohs),length(deltas)+1,length(hdgs));
                               % add extra column^ for pooling all trials irrespective of delta
pRight = n;
pCorrect = n;

RTmean = n; RTse = n;
confMean = n; confSE = n;
xVals = hdgs(1):0.1:hdgs(end);
yVals = nan(length(mods),length(cohs),length(deltas)+1,length(xVals));

B = cell(length(mods),length(cohs),length(deltas)+1);
stats = cell(length(mods),length(cohs),length(deltas)+1);
plotLogistic = nan(length(mods),length(cohs),length(deltas)+1);

for m = 1:length(mods)
for c = 1:length(cohs)
for d = 1:length(deltas)+1 % add extra column for all trials irrespective of delta
        
    for h = 1:length(hdgs)
        if d==length(deltas)+1
            J = data.modality==mods(m) & data.coherence==cohs(c) & data.heading==hdgs(h); % all trials irrespective of delta
        else
            J = data.modality==mods(m) & data.coherence==cohs(c) & data.heading==hdgs(h) & data.delta==deltas(d);
        end
        
        n(m,c,d,h) = nansum(J);
%         if sum(J), keyboard, end
        pRight(m,c,d,h) = nansum(J & data.choice==2) / n(m,c,d,h); % 2 is rightward!!       
        pCorrect(m,c,d,h) = nansum(J & data.correct) / n(m,c,d,h);
        
        if RTtask
            if RTCorrOnly
                RTmean(m,c,d,h) = nanmean(data.RT(J & use_trs));
                RTse(m,c,d,h) = nanstd(data.RT(J))/sqrt(nansum(J & use_trs));
            else
                RTmean(m,c,d,h) = nanmean(data.RT(J));
                RTse(m,c,d,h) = nanstd(data.RT(J))/sqrt(n(m,c,d,h));
            end
        else
            RTmean(m,c,d,h) = NaN;
            RTse(m,c,d,h) = NaN;
        end
        
        if conftask==1 % saccEndpoint
            
            if RTCorrOnly
                confMean(m,c,d,h) = nanmean(data.conf(J & use_trs));
                confSE(m,c,d,h) = nanstd(data.conf(J & use_trs))/sqrt(sum(J & use_trs));
            else
                confMean(m,c,d,h) = nanmean(data.conf(J));
                confSE(m,c,d,h) = nanstd(data.conf(J))/sqrt(n(m,c,d,h));
            end
            
        elseif conftask==2 % PDW
            % ignore 1-target trials here
            if isfield(data,'oneTargConf')
                %disp('ignoring one-target confidence trials')
                J = J & ~data.oneTargConf;
            end
            
            if RTCorrOnly
                confMean(m,c,d,h) = nansum(J & data.PDW==1 & use_trs) / sum(J & use_trs);
            else
                confMean(m,c,d,h) = nansum(J & data.PDW==1) / sum(J); % 1 is high
            end
            % SE gets calculated below
            
        else % no conf at all
            confMean(m,c,d,h) = NaN;
            confSE(m,c,d,h) = NaN;
        end
    end

    % fit logistic regression
    if d==length(deltas)+1
        K = data.modality==mods(m) & data.coherence==cohs(c); % all trials irrespective of delta
    else
        K = data.modality==mods(m) & data.coherence==cohs(c) & data.delta==deltas(d);
    end
    if nansum(~isnan(data.heading(K)))>=3*length(hdgs) && length(hdgs)>5
        X = data.heading(K);
        y = data.choice(K)==2; % 2 is rightward
        [B{m,c,d}, ~, stats{m,c,d}] = glmfit(X, y, 'binomial');
        yVals(m,c,d,:) = glmval(B{m,c,d},xVals,'logit');
        plotLogistic(m,c,d) = 1;
    else
        plotLogistic(m,c,d) = 0;
    end

end
end
end

% standard error of proportion
pRightSE = sqrt( (pRight.*(1-pRight)) ./ n );
pCorrectSE = sqrt( (pCorrect.*(1-pCorrect)) ./ n );

if conftask==2 % PDW
    confSE = sqrt( (confMean.*(1-confMean)) ./ n );
end

% copy vestib-only data to all coherences, to aid plotting
for c=1:length(cohs)
    n(1,c,:,:) = n(1,1,:,:);
    pRight(1,c,:,:) = pRight(1,1,:,:);
    pRightSE(1,c,:,:) = pRightSE(1,1,:,:);
    pCorrect(1,c,:,:) = pCorrect(1,1,:,:);
    pCorrectSE(1,c,:,:) = pCorrectSE(1,1,:,:);
    confMean(1,c,:,:) = confMean(1,1,:,:);
    confSE(1,c,:,:) = confSE(1,1,:,:);
    RTmean(1,c,:,:) = RTmean(1,1,:,:);
    RTse(1,c,:,:) = RTse(1,1,:,:);
    yVals(1,c,:,:) = yVals(1,1,:,:);
    plotLogistic(1,c,:) = plotLogistic(1,1,:);
    B(1,c,:) = B(1,1,:);
    stats(1,c,:) = stats(1,1,:);
end

parsedData = struct();
% parsedData.subject = data.filename{1}(1:5);
parsedData.n = n;
parsedData.pRight = pRight;
parsedData.pRightSE = pRightSE;
parsedData.pCorrect = pCorrect;
parsedData.pCorrectSE = pCorrectSE;
parsedData.xVals = xVals;
parsedData.yVals = yVals;
parsedData.plotLogistic = plotLogistic;
parsedData.B = B;
parsedData.stats = stats;

if conftask
    parsedData.confMean = confMean;
    parsedData.confSE = confSE;
end

if RTtask
    parsedData.RTmean = RTmean;
    parsedData.RTse = RTse;
end
