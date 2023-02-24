function [dotX_3D,dotY_3D,dotZ_3D,dotSize,dotX,dotY] = generateDots3D_offline(PDS)

% generate 3D motion dots offline, using random seed from trials

% requires CreateUniformDotsIn3DFrustum and RandLim, and access to PDS
% struct fields

%  SJ 01-2022 

% is passing in the entire PDS struct the most useful way to do it?

% presumably this would be an intermediate step to running some kind of
% analysis on the dots e.g. ME, which would then get saved to parallel the
% cleaned data Struct.
% we don't want to save this back into the PDS struct on the server,
% because that defeats the point of removing the dots in the first place.

% make the worker function do just one trial, and then wrapper loops over trials,
% skipping vestib ones and storing dots data in file corresponding to main PDS one which can then be
% merged into dataStruct later if desired?

% then test that all trials from recording are recoverable
% test with mixed mods (including handling vestib trials)


%  dotsInPlaneOfScreen=0;

viewdist = PDS.initialParameters{3}.display.viewdist;
dist = viewdist*10; % convert to mm


if isempty(PDS.initialParameters{3}.stimulus.eyeHeight) %|| p.trial.stimulus.eyeHeight==0
    tempEyeHeight = [];
else
    tempEyeHeight = 0;
    %     tempEyeHeight = p.trial.stimulus.eyeHeight;
end

%%
for t = 1:length(PDS.data)
    
    if PDS.conditions{t}.stimulus.modality==1
        % allocate memory for dot positions and sizes
        dotX_3D = nan;
        dotY_3D = nan;
        dotZ_3D = nan;
        dotSize = nan;
    else
    
    traj = PDS.data{t}.stimulus.visTraj;

    % allocate memory for dot positions and sizes
    dotX_3D{t} = nan(size(traj,1),PDS.initialParameters{3}.stimulus.ndots);
    dotY_3D{t} = nan(size(traj,1),PDS.initialParameters{3}.stimulus.ndots);
    dotZ_3D{t} = nan(size(traj,1),PDS.initialParameters{3}.stimulus.ndots);
    dotSize{t} = nan(size(traj,1),PDS.initialParameters{3}.stimulus.ndots);
    
    % set the trialRNG back using the trialSeed!
    reset(PDS.data{t}.stimulus.rngs.trialRNG,PDS.data{t}.stimulus.rngs.trialSeed)
    
    % Generate the first frame
    
    [x,y,z] = CreateUniformDotsIn3DFrustum(PDS.initialParameters{3}.stimulus.ndots, PDS.initialParameters{3}.stimulus.fovy, 1/PDS.initialParameters{4}.stimulus.ar, ...
        PDS.initialParameters{3}.stimulus.clipNear*10, PDS.initialParameters{3}.stimulus.clipFar*10, tempEyeHeight,PDS.data{t}.stimulus.rngs.trialRNG);
    % distances in mm
    
    z = z + dist;
    
    % the first frame is x,y,z, regardless of coherence
    dotX_3D{t}(1,:) = x;
    dotY_3D{t}(1,:) = y;
    dotZ_3D{t}(1,:) = z;
    
    tx = -traj(1,4); ty = traj(1,3); tz = -traj(1,5); % translation
    
    % update dot size for every frame, based on Z distance to eye plane
    %     distToEye = sqrt((x+tx).^2 + (y+ty).^2 + ((z+tz)-p.trial.display.viewdist*10).^2); % mm
    distToEyePlane = abs((z+tz)-dist); % mm
    dotSize{t}(1,:) = viewdist./(distToEyePlane/10) *  PDS.initialParameters{3}.stimulus.dotSizeInWorld * PDS.initialParameters{4}.display.ppcm;
    
    velProfile = PDS.data{t}.stimulus.velProfile;
    if ~isempty(velProfile)
        
        %***********************************************************
        % scale dot lifetime by inverse velocity, normalized to
        % the range specified by min/maxDotLifetime
        %***********************************************************
        
        coh = PDS.conditions{t}.stimulus.coherence;
        
        velProfileNorm = (velProfile-min(velProfile)) / max(velProfile-min(velProfile));
        
        velProfileNorm(velProfileNorm<0.05) = 0.05; % avoid div by 0
        flippedVel = 1./velProfileNorm;
        % OR
        % flippedVel = 1-velProfileNorm;
        
        flippedVelNorm = (flippedVel-min(flippedVel)) / max(flippedVel-min(flippedVel));
        
        lifetime = flippedVelNorm * (PDS.initialParameters{3}.stimulus.maxDotLifetime-PDS.initialParameters{3}.stimulus.minDotLifetime) + PDS.initialParameters{3}.stimulus.minDotLifetime;
        lifetime(end+1:length(traj)) = lifetime(end); % pad the end, like traj
        
        % kluge long lifetimes so there are no conspicuous dot jumps early/late
        % here = length(lifetime)-find(flipud(lifetime)<max(flipud(lifetime)),1,'first') + 2;
        % lifetime(here:end) = 999;
        
        % kluge shift lifetime curve by ~15 frames
        % lifetime2 = lifetime(15:end);
        % lifetime2(end+1:end+15) = lifetime(end);
        % lifetime = lifetime2;
        
        % OR instead of shifting, stretch lifetime profile using resample()
        lifetime2 = resample(lifetime,round(length(lifetime)*1.2),length(lifetime));
        clip = length(lifetime2)-length(lifetime);
        lifetime2(1 : (round(clip/2)-1)) = [];
        lifetime2(end-(round(clip/2)-1) : end) = [];
        % *** WARNING: MUST BE CHECKED WHENEVER MOTION PROFILE CHANGES
        lifetime = lifetime2;
        
        % figure; plot(traj(:,5)+250); hold on; plot(velProfile,'g'); plot(lifetime,'r');
        
    else % otherwise, use infinite lifetime (becomes 100% coh)
        lifetime = ones(size(traj,1),1)*200;
        coh = 1;
    end
    
    % now implement coherence for the remaining frames
    % this step occurs in world coordinates with no 'camera' movement at all;
    % it's simply replotting a random subset of dots on each frame, or every
    % N frames, dictated by lifetime.
    ndots = PDS.initialParameters{3}.stimulus.ndots;
    
    for f = 2:size(traj,1)
        %     disp(num2str(mod(f-1,round(lifetime(f)))));
        
        % SJ 10/2020 shouldn't this line be here to update dotSize every frame?
        tx = -traj(f,4); ty = traj(f,3); tz = -traj(f,5); % translation
        
        if mod(f-1,round(lifetime(f)))==0
            
            %         disp([num2str(f) '  ' num2str(lifetime(f))]);
            
            %find the index of 'noise' dots in this frame (or set of N=lifetime frames)
            L = rand(PDS.data{t}.stimulus.rngs.trialRNG,ndots,1) > abs(coh);
            if sum(L)>0
                % make new random locations for the noise dots
                [xnoise,ynoise,znoise] = CreateUniformDotsIn3DFrustum(sum(L), PDS.initialParameters{3}.stimulus.fovy, 1/PDS.initialParameters{4}.stimulus.ar, ...
                    PDS.initialParameters{3}.stimulus.clipNear*10, PDS.initialParameters{3}.stimulus.clipFar*10, tempEyeHeight,PDS.data{t}.stimulus.rngs.trialRNG);
                % distances in mm
                znoise = znoise + dist;
                
                
                x(L) = xnoise;
                y(L) = ynoise;
                z(L) = znoise;
            end
            
        end
        
        % save dot positions (again, world-fixed, not on screen or even
        % relative to moving observer -- those can be calculated from these)
        dotX_3D{t}(f,:) = x;
        dotY_3D{t}(f,:) = y;
        dotZ_3D{t}(f,:) = z;
        
        % update dot size for every frame, based on *Z* distance to eye plane
        % (not euclidian distance?)
        %     distToEyePlane = sqrt((x+tx).^2 + (y+ty).^2 + ((z+tz)-p.trial.display.viewdist*10).^2); % mm
        distToEyePlane = abs((z+tz)-dist); % mm
        dotSize{t}(f,:) = viewdist./(distToEyePlane/10) *  PDS.initialParameters{3}.stimulus.dotSizeInWorld * PDS.initialParameters{4}.display.ppcm;
        
    end
    
    % screen coords of dots SJ 06-2022
    
    X = bsxfun(@plus,dotX_3D, -traj(:,4));
    Y = bsxfun(@plus,dotY_3D, traj(:,3));
    Z = bsxfun(@plus,dotZ_3D, -traj(:,5));
    
    % Z minus dist so that z is relative to the eye again!
    dotX{t} = dist * (X./-(Z-dist)) * PDS.initialParameters{4}.display.ppcm/10;
    dotY{t} = dist * (Y./-(Z-dist)) * PDS.initialParameters{4}.display.ppcm/10;
    
    % add display center to get to raw screen pixels?
    dotX{t} = dotX{t} + PDS.initialParameters{4}.display.ctr(1);
    dotY{t} = -dotY{t} + PDS.initialParameters{4}.display.ctr(2); % negate y to get to screen co-ords

    end
    
end
