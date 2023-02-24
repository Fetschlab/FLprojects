function parsedData = dots3DMP_parseData_byConf(data,mods,cohs,deltas,hdgs,conftask,RTtask)
% SJ 07-2021 converted to function, for cleaner workspace
% does exactly the same as regular parseData, but splits all output
% variables for high and low confidence
% in SEP confidence task, high and low is determined by median split across
% all confidence ratings - use caution with this, particularly if conf
% ratings are unnormalized!

if nargin < 7, RTtask = 0; end
if conftask == 0, error('Cannot split data by confidence for non-conf task!'); end

% switch back to 1:2 if 0:1 (e.g. from fitting code)
if max(data.choice)==1
    data.choice = data.choice+1;
end

%% parse data
% create and use matrices of summary data indexed by variables of interest

n = nan(length(mods)+1,length(cohs),length(deltas)+1,length(hdgs),2);
                               % add extra column^ for pooling all trials irrespective of delta
[pRight, RTmean, RTse, confMean, confSE] = deal(n);

xVals = hdgs(1):0.1:hdgs(end);
yVals = nan(length(mods)+1,length(cohs),length(deltas)+1,length(xVals),2);

B = cell(length(mods)+1,length(cohs),length(deltas)+1,2);
stats = cell(length(mods)+1,length(cohs),length(deltas)+1,2);
plotLogistic = nan(length(mods)+1,length(cohs),length(deltas)+1,2);

if conftask==1
    % currently just a median split across subjs
    % should be fine if we have normalized within subj and within L/R,
    % otherwise need to be careful...
    hiConf = data.conf >= median(data.conf);
elseif conftask==2 
    hiConf = data.PDW;
end

for m = 1:length(mods)+1
for c = 1:length(cohs)
for d = 1:length(deltas)+1 % add extra column for all trials irrespective of delta

    for h = 1:length(hdgs)
        if m==length(mods)+1
            J = data.heading==hdgs(h);
        elseif d==length(deltas)+1
            J = data.modality==mods(m) & data.coherence==cohs(c) & data.heading==hdgs(h); % all trials irrespective of delta
        else
            J = data.modality==mods(m) & data.coherence==cohs(c) & data.heading==hdgs(h) & data.delta==deltas(d);
        end

        Jhi = J & hiConf & data.oneTargConf==0;
        Jlo = J & ~hiConf & data.oneTargConf==0;
        
        n(m,c,d,h,1) = sum(Jhi);
        n(m,c,d,h,2) = sum(Jlo);

        pRight(m,c,d,h,1) = sum(Jhi & data.choice==2) / n(m,c,d,h,1); % 2 is rightward
        pRight(m,c,d,h,2) = sum(Jlo & data.choice==2) / n(m,c,d,h,2); % 2 is rightward

        if RTtask
            RTmean(m,c,d,h,1) = mean(data.RT(Jhi));
            RTmean(m,c,d,h,2) = mean(data.RT(Jlo));
            
            RTse(m,c,d,h,1) = std(data.RT(Jhi))/sqrt(n(m,c,d,h,1));
            RTse(m,c,d,h,2) = std(data.RT(Jlo))/sqrt(n(m,c,d,h,2));
        end
        
        if conftask==1
            confMean(m,c,d,h,1) = mean(data.conf(Jhi));
            confSE(m,c,d,h,1) = std(data.conf(Jhi))/sqrt(n(m,c,d,h,1));
            
            confMean(m,c,d,h,2) = mean(data.conf(Jlo));
            confSE(m,c,d,h,2) = std(data.conf(Jlo))/sqrt(n(m,c,d,h,2));
            
        else % PDW
            confMean(m,c,d,h,1) = mean(data.PDW(Jhi)); % 1 is high
            confMean(m,c,d,h,2) = mean(data.PDW(Jlo)); 

            % SE gets calculated below
        end            
    end

    % fit logistic regression
    if m==length(mods)+1
        K = true(size(data.modality));
    elseif d==length(deltas)+1
        K = data.modality==mods(m) & data.coherence==cohs(c); % all trials irrespective of delta
    else
        K = data.modality==mods(m) & data.coherence==cohs(c) & data.delta==deltas(d);
    end
    Khi = K & hiConf & data.oneTargConf==0;

    if sum(~isnan(data.heading(Khi)))>=3*length(hdgs) && length(hdgs)>5
        X = data.heading(Khi);
        y = data.choice(Khi)==2; % 2 is rightward
        [B{m,c,d,1}, ~, stats{m,c,d,1}] = glmfit(X, y, 'binomial');
        yVals(m,c,d,:,1) = glmval(B{m,c,d,1},xVals,'logit');
        plotLogistic(m,c,d,1) = 1;
    else
        plotLogistic(m,c,d,1) = 0;
    end
    
    Klo = K & ~hiConf & data.oneTargConf==0;
    if sum(~isnan(data.heading(Klo)))>=3*length(hdgs) && length(hdgs)>5
        X = data.heading(Klo);
        y = data.choice(Klo)==2; % 2 is rightward
        [B{m,c,d,2}, ~, stats{m,c,d,2}] = glmfit(X, y, 'binomial');
        yVals(m,c,d,:,2) = glmval(B{m,c,d,2},xVals,'logit');
        plotLogistic(m,c,d,2) = 1;
    else
        plotLogistic(m,c,d,2) = 0;
    end

end
end
end

pRightSE = sqrt( (pRight.*(1-pRight)) ./ n );
if conftask==2
    confSE = sqrt( (confMean.*(1-confMean)) ./ n );
end

% copy vestib-only data to all coherences, to aid plotting
for c=1:length(cohs)
    n(1,c,:,:,:) = n(1,1,:,:,:);
    pRight(1,c,:,:,:) = pRight(1,1,:,:,:);
    pRightSE(1,c,:,:,:) = pRightSE(1,1,:,:,:);
    confMean(1,c,:,:,:) = confMean(1,1,:,:,:);
    confSE(1,c,:,:,:) = confSE(1,1,:,:,:);
    RTmean(1,c,:,:,:) = RTmean(1,1,:,:,:);
    RTse(1,c,:,:,:) = RTse(1,1,:,:,:);
    confMean(1,c,:,:,:) = confMean(1,1,:,:,:);
    confSE(1,c,:,:,:) = confSE(1,1,:,:,:);
    yVals(1,c,:,:,:) = yVals(1,1,:,:,:);
    plotLogistic(1,c,:,:) = plotLogistic(1,1,:,:);
    B(1,c,:,:) = B(1,1,:,:);
    stats(1,c,:,:) = stats(1,1,:,:);
end


parsedData = struct();
%parsedData.subject = data.filename{1}(1:5);
parsedData.n = n;
parsedData.pRight = pRight;
parsedData.pRightSE = pRightSE;
parsedData.xVals = xVals;
parsedData.yVals = yVals;
parsedData.plotLogistic = plotLogistic;
parsedData.B = B;
parsedData.stats = stats;

parsedData.confMean = confMean;
parsedData.confSE = confSE;


if RTtask
    parsedData.RTmean = RTmean;
    parsedData.RTse = RTse;
end

