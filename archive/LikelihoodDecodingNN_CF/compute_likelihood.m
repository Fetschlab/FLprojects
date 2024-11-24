% Computes likelihood functions from resampled neuronal responses (found in
% all_cells.mat).  This function is called from within population_likelihood.m

function [r,L] = compute_likelihood(unique_heading,resp_trial,resp_mean_all_temp,s,j,k)

g=1;
hdg = unique_heading(1):0.1:unique_heading(end);
sizeR = size(resp_trial);
N = sizeR(1);
reps = sizeR(2);
r = resp_trial(:,round(rand*(reps-1))+1);

% tuning curve is either tailored for each condition (comb gets zero-conflict),
tunings = squeeze(resp_mean_all_temp(:,s,2,k,:));
% or, comb for all conditions (zero-conflict, vary coherence):
%     tunings = squeeze(resp_mean_all_temp(:,3,2,k,:));
% or, comb for all conditions (zero-conflict, fix coherence):
%     tunings = squeeze(resp_mean_all_temp(:,3,2,1,:));
% or, vestibular for all conditions
% tunings = squeeze(resp_mean_all_temp(:,1,2,1,:));

f = interp1(unique_heading',tunings',hdg,'linear');
f=f';

% % temp: shrink r and f to solve inf problem (happens for high fakedata weights >=0.9)
%f = f * 0.2;
%r = round(r*0.2);

flag_zero_problem = [];
for n = 1:N
%     if N ~=      % just skip the problematic cells and see if the result changes
        
    if sum((g*f(n,:)).^r(n)) == Inf || factorial(r(n)) == Inf
%           disp(['Inf problem! Cell ID #' num2str(n)]);
       f(n,:) = f(n,:) * 0.8; % This kluges one cell whose high FR makes the equation exceed
       r(n) = round(r(n)*0.8); % the largest representable positive floating point number
    end
    if sum((g*f(n,:)).^r(n)) == Inf || factorial(r(n)) == Inf
        disp(['Inf problem after scale-down! Cell ID #' num2str(n)]);
    end
    if sum(f(n,:)) == 0
%         disp(['Zero problem: Cell ID #' num2str(n) ', s,j,k = ' num2str([s j k])]);
        flag_zero_problem = [flag_zero_problem n];
        f(n,:) = 1;
    end
    
    % individual neuron likelihood (poisson):
    temp_L(n,:) = ( exp(-g*f(n,:)) .* (g*f(n,:)).^r(n) ) / factorial(r(n));
%     figure(11); plot(hdg,temp_L(n,:)); title(num2str(n)); pause;
    
end

if ~isempty(flag_zero_problem)
    disp(['Zero problem: cell #s ' num2str(flag_zero_problem)]);
end

L = prod(temp_L);

if sum(L) == 0
    L = prod(temp_L+10e-10);
end


return
