function [tuning_curves, posterior, pop_lh, hdg_intp, simdata] = dots3DMP_likelihood_decoding(data,meanFRs,hdgs,mods,cohs,deltas,numtrs,step)

% numtrs - number of simulated trials
% step - tuning curve stepsize

lh_func = @(g,tc,r_obs) (exp(-g*tc) .* bsxfun(@power,g*tc,r_obs)) ./ factorial(r_obs); % vectorized

numunits = length(unique(data.unitnum));
hdg_intp = min(hdgs):step:max(hdgs);

[hdg,modality,coh,delta,ntrials] = dots3DMP_create_trial_list(hdgs',mods',cohs',deltas',numtrs,0);
simdata.heading = hdg;
simdata.modality = modality;
simdata.coherence = coh;
simdata.delta = delta;

% pre-allocate
tuning_curves = nan(length(mods),length(cohs),length(hdg_intp),numunits);
% r_obs         = nan(length(mods),length(cohs),length(deltas),length(hdgs),numunits,numtrs);
lh            = nan(length(simdata.heading),length(hdg_intp),numunits);

g = 1; % Poisson gain

for m=1:length(mods)
for c=1:length(cohs)
for u=1:numunits
    x = squeeze(meanFRs(m,c,deltas==0,:,u));   % use delta==0 only!
    if any(isnan(x)),continue,end % skip
    tuning_curves(m,c,:,u) = interp1(hdgs,x,hdg_intp,'linear');
 
    % linear regression
%     p = polyfit(hdgs,x,1);
%     y = polyval(p,xq);
%     tuning_curves(m,c,d,:,u) = y;
end
end
end

for m=1:length(mods)
for c=1:length(cohs)
for d=1:length(deltas)
for u=1:numunits

    % tuning curve for given stimulus
    tc = squeeze(tuning_curves(m,c,:,u)); 

    for h=1:length(hdgs)
        % random draws of single-trial firing for each condition
        % ... and compute likelihood

        J = data.modality==mods(m) & data.coherence==cohs(c) & data.heading==hdgs(h) & data.delta==deltas(d) & data.unitnum==u;
        
        % out of all actual trials sum(J) meet criteria. what was spikeRate
        % on this subset of trials?
        r_trials = round(data.spikeRates(J));
        if isempty(r_trials), continue, end
        shuf_ind = randi(length(r_trials),1,numtrs);
        clear r

        r(:,1) = r_trials(shuf_ind); % Spike rates assigned to 100 simTrials by sampling from actual observed trial spike rates stored in 'r_trials'
        %         r_obs(m,c,d,h,u,:) = r; % Seemingly unused..?

        if any(sum(bsxfun(@power,g*tc,r')) == Inf) || any(factorial(r) == Inf)
            %                 disp(['Inf problem! Cell ID #' num2str(n)]);
            tc = tc .* 0.8;  % This kluges one cell whose high FR makes the equation exceed
            r = round(r.*0.8); % the largest representable positive floating point number
        end
        if any(sum(bsxfun(@power,g*tc,r')) == Inf) || any(factorial(r) == Inf)
            %             disp(['Inf problem after scale-down! Cell ID #' num2str(n)]);
            tc = tc .* 0.5;  % One(?) cell needs even more of a scale-down
            r = round(r.*0.5);
        end
        if any(sum(bsxfun(@power,g*tc,r')) == Inf) || any(factorial(r) == Inf)
            disp(['Inf problem after double scale-down! Cell ID #' num2str(n)]);
            keyboard
        end
        if sum(tc) == 0
            disp(['Zero problem: Cell ID #' num2str(n)]);
            tc = ones(size(tc)); % cell does not contribute, but does not zero out the likelihood either
        end

        % boolean index for condition in simulated trials data struct 
        K = simdata.modality==mods(m) & simdata.coherence==cohs(c) & simdata.heading==hdgs(h) & simdata.delta==deltas(d);
        lh(K,:,u) = lh_func(g,tc,r')';

    end
end
end
end
end

lh(lh==0) = .01;               % kluge to avoid prod(0...)
pop_lh = squeeze(prod(lh,3, 'omitnan')); % population likelihood, product across units
% pop_lh = squeeze(prod(lh,6,'omitnan'));   % population likelihood, product across units

% normalize over interpolated headings to get posterior, assuming flat prior
% posterior = pop_lh ./ sum(pop_lh,5,'omitnan');
posterior = pop_lh ./ sum(pop_lh,2,'omitnan');


%%

% [~,zerohdg] = min(abs(hdg_intp));
zerohdg = find(hdg_intp==0);

negSum = nansum(posterior(:,1:zerohdg-1),2);
posSum = nansum(posterior(:,zerohdg+1:end),2);
simdata.choice = double(posSum>negSum)+1;

% [~,maxpos] = max(posterior,[],5);
% MAPhdg     = squeeze(hdg_intp(maxpos));
% simChoice  = sign(MAPhdg);
% simChoice(simChoice==0) = sign(randn);

