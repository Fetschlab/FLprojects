
% for tuning functions
tuning_vonMises = @(b,dir) b(1) * exp(b(2)*cosd(dir-b(3))) / (2*pi*besseli(0,b(2))) + b(4);
tuning_vonMises_err = @(b,dir,FR) nansum((tuning_vonMises(b,dir)-FR).^2); % sum squared error

% check max trial len, nunits, ntrials (for preallocation)
nU = nan(1,length(dataCell));
nTrials_tuning = nan(1,length(dataCell));
nTrials_exp = nan(1,length(dataCell));
pct99_tuning = nan(1,length(dataCell));
pct99_exp = nan(1,length(dataCell));
maxlen_tuning = nan(1,length(dataCell));
maxlen_exp = nan(1,length(dataCell));
for n = 1:length(dataCell)
    
%     % tmp    
%     any(isnan(dataCell{1}.data.dots3DMP.events.PDW) ~= isnan(dataCell{1}.data.dots3DMP.events.choice))
%     invalidTr = isnan(dataCell{1}.data.dots3DMP.events.choice);
        % should be covered by .goodtrial
    
    if max(dataCell{n}.data.dots3DMP.events.choice)==2 % convert choices to 0:1
        dataCell{n}.data.dots3DMP.events.choice = dataCell{n}.data.dots3DMP.events.choice-1;
    end
    
%     goodtr = dataCell{n}.data.dots3DMPtuning.events.goodtrial; % field goodtrial seems not valid for tuning
    goodtr = isnan(dataCell{n}.data.dots3DMPtuning.events.breakfix);  
    mStart = dataCell{n}.data.dots3DMPtuning.events.stimOn(goodtr)';
    mEnd = dataCell{n}.data.dots3DMPtuning.events.stimOff(goodtr)';
    dur = round((mEnd-mStart)*1000);
    pct99_tuning(n) = ceil(prctile(dur,99)); % find plausible max dur
    maxlen_tuning(n) = max(dur); % or actual max dur
    
    goodtr = dataCell{n}.data.dots3DMP.events.goodtrial;
    mStart = dataCell{n}.data.dots3DMP.events.stimOn(goodtr)';
    mEnd = dataCell{n}.data.dots3DMP.events.stimOff(goodtr)';
    dur = round((mEnd-mStart)*1000);
    pct99_exp(n) = ceil(prctile(dur,99)); % find plausible max dur
    maxlen_exp(n) = max(dur); % or actual max dur

    % sanity checks
    if any(isnan(dataCell{n}.data.dots3DMP.events.choice(goodtr)))
        keyboard
    end
    if any(isnan(dataCell{n}.data.dots3DMP.events.PDW(goodtr)))
        keyboard
    end    
    if length(dataCell{1}.data.dots3DMPtuning.spiketimes) ~= length(dataCell{1}.data.dots3DMP.spiketimes)
        keyboard
    end    
    
    nU(n) = length(dataCell{n}.data.dots3DMPtuning.spiketimes);
    nTrials_tuning(n) = sum(dataCell{n}.data.dots3DMPtuning.events.goodtrial);
    nTrials_exp(n) = sum(dataCell{n}.data.dots3DMP.events.goodtrial);
end
rasterLen_tuning = max(pct99_tuning)+1; % +1 for a cushion, some weird edge effects
% rasterLen_exp = max(pct99_exp)+1;
% rasterLen_tuning = max(maxlen_tuning)+1;
rasterLen_exp = max(maxlen_exp)+1;

nUnit = sum(nU) - size(excludes,1);
sessID = nan(nUnit,1); % keep track of which session a given unit was in

% preallocate

% raster: [unit,trials,time]
raster_tuning =    nan(nUnit,max(nTrials_tuning),rasterLen_tuning+sum(extRaster));
raster_dotsOn = nan(nUnit,max(nTrials_exp),rasterLen_exp+sum(extRaster));
raster_RT =     nan(nUnit,max(nTrials_exp),rasterLen_exp+sum(extRaster));

% aggregated, normalized psth [unit,conditions(x4),time]:
psthNorm_dotsOn_one =   nan(nUnit,rasterLen_exp+sum(extRaster));
psthNorm_dotsOn_two =   nan(nUnit,rasterLen_exp+sum(extRaster));
psthNorm_dotsOn_three = nan(nUnit,rasterLen_exp+sum(extRaster));
psthNorm_dotsOn_four =  nan(nUnit,rasterLen_exp+sum(extRaster));
psthNorm_RT_one =       nan(nUnit,rasterLen_exp+sum(extRaster));
psthNorm_RT_two =       nan(nUnit,rasterLen_exp+sum(extRaster));
psthNorm_RT_three =     nan(nUnit,rasterLen_exp+sum(extRaster));
psthNorm_RT_four =      nan(nUnit,rasterLen_exp+sum(extRaster));

% spike rates/counts: [unit,trials]
% % % spRate_tuning =  nan(nUnit,max(nTrials_tuning));
spRate_exp =  nan(nUnit,max(nTrials_exp));
spCount_exp = nan(nUnit,max(nTrials_exp));

% misc
prefDir_fit = nan(nUnit,1);
CP_high =     nan(nUnit,1);
CP_low =      nan(nUnit,1);
CP_all =      nan(nUnit,1);
ConfP_pref =  nan(nUnit,1);
ConfP_null =  nan(nUnit,1);
ConfP_all =   nan(nUnit,1);

