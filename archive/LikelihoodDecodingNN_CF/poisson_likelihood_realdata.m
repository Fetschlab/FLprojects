% Computes likelihood functions from resampled neuronal responses (found in
% all_cells.mat).  This function is called from within population_likelihood.m

function L = poisson_likelihood_realdata(s,j,k,i)


% load(['C:\Program Files\MATLAB\R2007a\work\neurons_indiv_trials\all_cells.mat']);
% load(['C:\Program Files\MATLAB\R2007a\work\neurons_indiv_trials\good_cells_only.mat']);

% FORMAT: resp_trial_all{n}{s,j,k,i} , resp_mean_all{n}{s,j,k}

N = length(resp_trial_all);
% hdg = -20:0.1:20;
hdg = -10:0.1:10;
heading = [-10 -3.5 -1.225 0 1.225 3.5 10];
g = 1;

i_index = find(round(hdg*10) == round(heading(i)*10)); % fix matlab rounding problem
i_index = i_index(1);

f = zeros(N,length(hdg));
r = zeros(N,1);
temp_L = zeros(N,length(hdg));

% for qq = 1:20
for n = 1:N
    
    % Jazayeri and Movshon
%     n_reps = length(resp_trial_all{n}{s,j,k,i});
%     n_shuf = randperm(n_reps);
%     r(n) = round(resp_trial_all{n}{s,j,k,i}(n_shuf(1)));
%     f(n) = mean(resp_trial_all{n}{s,j,k,i})
%     logL_n = r(n)*log(f(n)) - f(n) - log(factorial(r(n)));

    % Ma et al. (Seung and Somp, Sanger etc.)
    n_reps = length(resp_trial_all{n}{s,j,k,i});
    n_shuf = randperm(n_reps);
    r(n) = round(resp_trial_all{n}{s,j,k,i}(n_shuf(1)));
    % tuning curve is forced to be the zero-conflict condition.  This is important!
    f(n,:) = interp1(heading,resp_mean_all{n}{s,2,k},hdg,'linear');
    if sum((g*f(n,:)).^r(n)) == Inf
        f(n,:) = f(n,:) * 0.8; % This kluges one cell whose FR makes the equation exceed matlab's maximum number size
        r(n) = round(r(n)*0.8); % Check if that can be increased!
    end
%     figure(10); plot(hdg,f(n,:),'-','LineWidth',2,'Color',colorMap(n,:)); hold on;
	temp_L(n,:) = ( exp(-g*f(n,:)) .* (g*f(n,:)).^r(n) ) / factorial(r(n));

end

% ylabel('Firing Rate (sp/s)'); xlabel('Heading Angle (deg)');
% title('Vestibular');
% ylim([0 120]);


L = prod(temp_L);
% AAAAL = L'/sum(L');
% figure(11); plot(hdg,AAAAL); hold on;
% end %qq


% if sum(L) == Inf
%     asdf = 1;
% end

return