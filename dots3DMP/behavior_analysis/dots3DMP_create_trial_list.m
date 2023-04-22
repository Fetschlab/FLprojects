function [hdg,modality,coh,delta,ntrials] = dots3DMP_create_trial_list(hdgs,mods,cohs,deltas,nreps,shuff)

if nargin<6, shuff = 1; end

%% build trial list
% (can't just randsample the condition vectors, because certain combinations of
% modality, coh, delta etc are invalid)
    
% repeat heading list once for ves, ncoh for vis, and ncoh x ndelta for comb
numHdgGroups = any(ismember(mods,1)) + ...
               any(ismember(mods,2)) * length(cohs) + ...
               any(ismember(mods,3)) * length(cohs)*length(deltas);
hdg = repmat(hdgs', numHdgGroups, 1);

lh = length(hdgs);
coh = nan(size(hdg));
delta = nan(size(hdg));
modality = nan(size(hdg));

% kluge for ves, call it the lowest coh by convention
if any(ismember(mods,1))
    coh(1:lh) = cohs(1); 
    delta(1:lh) = 0;
    modality(1:lh) = 1;
    last = lh;
else
    last = 0;    
end

if any(ismember(mods,2)) % loop over coh for vis
    for c = 1:length(cohs)
        these = last+(c-1)*lh+1 : last+(c-1)*lh+lh;
        coh(these) = cohs(c);
        delta(these) = 0;
        modality(these) = 2;
    end
    last = these(end);
end

if any(ismember(mods,3)) % loop over coh and delta for comb
    for c = 1:length(cohs)
        for d = 1:length(deltas)
            here = last + 1 + (c-1)*lh*length(deltas) + (d-1)*lh;
            these = here:here+lh-1;
            coh(these) = cohs(c);
            delta(these) = deltas(d);
            modality(these) = 3;
        end
    end
end


% now replicate times nreps and shuffle (or not):
condlist = [hdg modality coh delta];
% trialTable = repmat(condlist,nreps,1); 
if size(condlist,2)~=4
    error('conditions lists should be column vectors!')
end

trialTable = repmat(condlist,nreps,1);
if shuff
    % why shuffle? well, sometimes it's easier to find particular trial
    % types to test out when debugging
    trialTable = Shuffle(trialTable,2);
end

hdg      = trialTable(:,1);  
modality = trialTable(:,2);  
coh      = trialTable(:,3);  
delta    = trialTable(:,4);  
ntrials  = size(trialTable,1);



