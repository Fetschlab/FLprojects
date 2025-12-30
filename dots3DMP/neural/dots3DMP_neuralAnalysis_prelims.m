% PRELIMS

% for tuning functions
tuning_vonMises = @(b,dir) b(1) * exp(b(2)*cosd(dir-b(3))) / (2*pi*besseli(0,b(2))) + b(4);
tuning_vonMises_err = @(b,dir,FR) nansum((tuning_vonMises(b,dir)-FR).^2); % sum squared error

tuning_line = @(b,dir) b(1) + b(2)*dir;
tuning_line_err = @(b,dir,FR) nansum((tuning_line(b,dir)-FR).^2); % sum squared error


% check max trial len, nunits, ntrials (for preallocation)
for n = 1:length(newDataStruct)

    ev_tun = newDataStruct(n).data.dots3DMPtuning.events;
    ev_exp = newDataStruct(n).data.dots3DMP.events;
    un = newDataStruct(n).data.dots3DMP.units;
    
    mStart = ev_tun.motionOn';
    mEnd = ev_tun.stimOff';
    dur = round((mEnd-mStart)*1000);
%     maxlen_tun(n) = round(max(dur));
    pct99_tun(n) = ceil(prctile(dur,99)); % use this instead, due to outliers
    
    mStart = ev_exp.motionOn';
    mEnd = ev_exp.postTargHold';
    dur = round((mEnd-mStart)*1000);
%     maxlen_exp(n) = round(max(dur));
    pct99_exp(n) = ceil(prctile(dur,99)); % use this instead, due to outliers
    
    nU(n) = length(un.spiketimes);
    nTrials_tun(n) = length(ev_tun.motionOn(ev_tun.goodtrial));
    nTrials_exp(n) = length(ev_exp.motionOn(ev_exp.goodtrial));
end
rasterLen_tun = max(pct99_tun)+1;
rasterLen_exp = max(pct99_exp)+1; % +1 for a cushion, some weird edge bugs

nUnit = sum(nU) - size(excludes,1);
sessID = nan(nUnit,1); % keep track of which session a given unit was in


% preallocate

% raster: [unit,trials,time]
raster_tun =    nan(nUnit,max(nTrials_tun),rasterLen_tun+sum(extRaster));
raster_stimOn = nan(nUnit,max(nTrials_exp),rasterLen_exp+sum(extRaster));
raster_RT =     nan(nUnit,max(nTrials_exp),rasterLen_exp+sum(extRaster));

% aggregated, normalized psth {modality}[unit,conditions(x4),time]:
for m=1:3
    psthNorm_stimOn_one{m} =   nan(nUnit,rasterLen_exp+sum(extRaster));
    psthNorm_stimOn_two{m} =   nan(nUnit,rasterLen_exp+sum(extRaster));
    psthNorm_stimOn_three{m} = nan(nUnit,rasterLen_exp+sum(extRaster));
    psthNorm_stimOn_four{m} =  nan(nUnit,rasterLen_exp+sum(extRaster));
    psthNorm_RT_one{m} =       nan(nUnit,rasterLen_exp+sum(extRaster));
    psthNorm_RT_two{m} =       nan(nUnit,rasterLen_exp+sum(extRaster));
    psthNorm_RT_three{m} =     nan(nUnit,rasterLen_exp+sum(extRaster));
    psthNorm_RT_four{m} =      nan(nUnit,rasterLen_exp+sum(extRaster));
end

% spike rates/counts: [unit,trials]
% % % spRate_tun =  nan(nUnit,max(nTrials_tun));
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

