function parsedData = dots3DMP_parseData_multiConf(data,mods,cohs,deltas,hdgs,confGroup,conftask,RTtask,removeOneTarg,splitPDW)
% SJ 07-2021 converted to function, for cleaner workspace
% does exactly the same as regular parseData, but splits all output
% variables by confGroup
% confGroup must be pre-specified by user, should be a vector matching
% length of vectors in data
% most basic would be just use PDW (high/low). 
% can also factor in 1-target trials
% or assign groupings based on 

if nargin < 10, splitPDW = 1; end
if nargin < 9, removeOneTarg = 0; end
if nargin < 8, RTtask = 0; end
if conftask == 0, error('Cannot split data by confidence for non-conf task!'); end

% need to remove one-target trials!
if removeOneTarg && isfield(data,'oneTargConf')
    fprintf('removing 1-target confidence trials\n')
    fnames = fieldnames(data);
    inds = logical(data.oneTargConf);
    for f=1:length(fnames)
        data.(fnames{f})(inds) = [];
    end
    confGroup(inds)=[];
end

if conftask==1
    hiConf = data.conf >= median(data.conf);
elseif conftask==2
    hiConf = data.PDW;
end
if isempty(confGroup)
    confGroup = hiConf;
end
groupSplit = length(unique(confGroup));


if splitPDW
    if ~removeOneTarg && isfield(data,'oneTargConf')
        hiConf(logical(data.oneTargConf))=2;
    end
    confGroup = confGroup + (double(hiConf)*groupSplit);
end

[cgs,ia,ic] = unique(confGroup);
ncg = length(cgs);

%% parse data
% create and use matrices of summary data indexed by variables of interest

n = nan(length(mods)+1,length(cohs),length(deltas)+1,length(hdgs),ncg);
                               % add extra column^ for pooling all trials irrespective of delta
[pRight, RTmean, RTse, confMean, confSE] = deal(n);

xVals = hdgs(1):0.1:hdgs(end);
yVals = nan(length(mods)+1,length(cohs),length(deltas)+1,length(xVals),ncg);

B = cell(length(mods)+1,length(cohs),length(deltas)+1,ncg);
stats = cell(length(mods)+1,length(cohs),length(deltas)+1,ncg);
plotLogistic = nan(length(mods)+1,length(cohs),length(deltas)+1,ncg);

for nc=1:ncg
thisConf = confGroup==cgs(nc); %ic==nc;

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
        
        Jc = J & thisConf;
        
        n(m,c,d,h,nc) = sum(Jc);
        
        pRight(m,c,d,h,nc) = sum(Jc & data.choice==2) / n(m,c,d,h,nc); % 2 is rightward
        
        if RTtask
            RTmean(m,c,d,h,nc) = mean(data.RT(Jc));
            RTse(m,c,d,h,nc) = std(data.RT(Jc))/sqrt(n(m,c,d,h,nc));
        end
        
        if conftask==1
            confMean(m,c,d,h,nc) = mean(data.conf(Jc));
            confSE(m,c,d,h,nc) = std(data.conf(Jc))/sqrt(n(m,c,d,h,nc));
            
        else % PDW
            confMean(m,c,d,h,nc) = mean(data.PDW(Jc)); % 1 is high
            % SE gets calculated below
        end
    end

    % fit logistic regression for choice data
    if m==length(mods)+1
        K = true(size(data.modality));
    elseif d==length(deltas)+1
        K = data.modality==mods(m) & data.coherence==cohs(c); % all trials irrespective of delta
    else
        K = data.modality==mods(m) & data.coherence==cohs(c) & data.delta==deltas(d);
    end
    K = K & confGroup==cgs(nc);

    if sum(~isnan(data.heading(K)))>=3*length(hdgs) && length(hdgs)>5
        X = data.heading(K);
        y = data.choice(K)==2; % 2 is rightward
        [B{m,c,d,nc}, ~, stats{m,c,d,nc}] = glmfit(X, y, 'binomial');
        yVals(m,c,d,:,nc) = glmval(B{m,c,d,nc},xVals,'logit');
        plotLogistic(m,c,d,nc) = 1;
    else
        plotLogistic(m,c,d,nc) = 0;
    end
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
parsedData.confGroups = cgs;
parsedData.confGroupSplit = groupSplit;
parsedData.confMean = confMean;
parsedData.confSE = confSE;


if RTtask
    parsedData.RTmean = RTmean;
    parsedData.RTse = RTse;
end

