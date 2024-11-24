function [tuning_curves, posterior, pop_lh, hdg_intp, simdata] = dots3DMP_likelihood_decoding(data,meanFRs,hdgs,mods,cohs,deltas,numtrs,step)

% numtrs - number of simulated trials
% step - tuning curve stepsize, how finely heading values will be sampled, over which tc is calculated

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
lh            = nan(length(simdata.heading),length(hdg_intp),numunits);

g = 1; % Poisson gain

for m=1:length(mods)
for c=1:length(cohs)
for u=1:numunits
    x = squeeze(meanFRs(m,c,deltas==0,:,u));   % use delta==0 only!
    if any(isnan(x)),continue,end % skip
    tuning_curves(m,c,:,u) = interp1(hdgs,x,hdg_intp,'linear'); % TC for each unit on each stim condition (m x c x all h)
 
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
        
        r_trials = round(data.spikeRates(J));        % out of all actual trials 'sum(J)' meet criteria. what is spikeRate on each of these matching condition trials?
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

lh(lh==0) = .01;  % kluge to avoid prod(0...), % lh dimensions:
% each trial condition repeated 100 times, Each trial is assigned a
% likelihood distribution of length 'hdgs_intp', each unit has a unique
% tr x lh_value matrix

pop_lh = squeeze(prod(lh,3,'omitnan')); % population likelihood, product across units

% normalize over interpolated headings to get posterior, assuming flat prior
posterior = pop_lh ./ sum(pop_lh,2,'omitnan');

%% Simulating choice behavior from decoded posterior

% area under the curve
% zerohdg = find(hdg_intp==0);
% 
% negSum = sum(posterior(:,1:(zerohdg-1)),2,'omitnan');
% posSum = sum(posterior(:,(zerohdg+1):end),2,'omitnan');
% simdata.choice = double(posSum>negSum)+1;

% Max estimate
[~,maxpos] = max(posterior,[],2);
MAPhdg     = hdg_intp(maxpos);
simChoice  = sign(MAPhdg);
simChoice(simChoice==0) = sign(randn); % assign L/R choice at random, for decoded 0 hdg trials
simChoice = simChoice*2; simChoice(simChoice<0) = 1; % Make choices 1 & 2 for L & R, respectively
simdata.choice = simChoice';

end


