% Computes likelihood functions from resampled neuronal responses (found in
% all_cells.mat).  This function is called from within population_likelihood.m
                                                                                       %stim(trial),deltnum,cohnum,hdgnum,
function [r,L] = compute_likelihood_backup(heading,resp_trial_subset,resp_trial_all,resp_mean_all,s,j,k,i,use_logL,param1_interp,param2_interp,likemodel)

dbstop if error
N = length(resp_mean_all);
hdg = heading(1):0.1:heading(end);
g = 1;

% heading_times_ten = round(heading*10);
% hdg_times_ten = round(hdg*10);
% i_index = find(hdg_times_ten==heading_times_ten(i));

f = zeros(N,length(hdg));
r = zeros(N,1);
temp_L = zeros(N,length(hdg));
term1 = zeros(N,length(hdg));

% n_rand = randperm(length(resp_trial_all{1}{s,j,k,i}));
remove = false(N,1); % to remove skipped cells later 
for n = 1:N % N is nNeurons

    %% for corr version, cannot shuffle trials! (cannot have different trial index for each cell)
    % r(n) = round(resp_trial_all{n}{s,j,k,i}(n_rand(1))); 

    n_shuf = randperm(length(resp_trial_subset{n}{s,j,k,i}));
    
    % skip cells where this particular condition yields too few trials,
    % which happens a lot when conditioning on previous trial
    if length(n_shuf)<4
        r(n) = NaN;
        f(n,:) = nan(size(f(n,:)));
        term1(n,:) = nan(size(term1(n,:)));
        temp_L(n,:) = nan(size(temp_L(n,:)));
        remove(n) = true;
        continue 
    end
    
    if strcmp(likemodel,'Poisson')
        r(n) = round(resp_trial_subset{n}{s,j,k,i}(n_shuf(1)));
    else
        r(n) = resp_trial_subset{n}{s,j,k,i}(n_shuf(1)); % no need to round for histogram or fitted likelihood! (in fact it screws it up)
    end
       
    % which tuning curve to use for the 'f' in likelihood eqn?
    % either tailored for each condition (comb gets zero-conflict),
    f(n,:) = interp1(heading,resp_mean_all{n}{s,2,k},hdg,'linear');
    % OR, comb for all conditions (zero-conflict, vary coherence):
    % f(n,:) = interp1(heading,resp_mean_all{n}{3,2,k},hdg,'linear');
    % OR, comb for all conditions (zero-conflict, fix coherence):
    % f(n,:) = interp1(heading,resp_mean_all{n}{3,2,1},hdg,'linear');
    % OR, vestibular for all conditions
    % f(n,:) = interp1(heading,resp_mean_all{n}{1,2,1},hdg,'linear');
    % Yong's suggestion: subtract sum of all conditions' tunings for 2nd term
    % F(n,:) = interp1(heading,resp_mean_all{n}{1,2,1},hdg,'linear') + interp1(heading,resp_mean_all{n}{2,2,2},hdg,'linear') + interp1(heading,resp_mean_all{n}{3,2,2},hdg,'linear');
    % or mean?
    % F(n,:) = mean([interp1(heading,resp_mean_all{n}{1,2,1},hdg,'linear') ; interp1(heading,resp_mean_all{n}{2,2,2},hdg,'linear') ; interp1(heading,resp_mean_all{n}{3,2,2},hdg,'linear')]);

    if use_logL && strcmp(likemodel,'Poisson') % log-likelihood (Dayan and Abbott, Jazayeri and Movshon):
    
        f(n,f(n,:)==0) = .01; % can't have zeros with log method!
        term1(n,:) = r(n)*log(f(n,:));
        % term2(n,:) = f(n,:);
        % term3(n,:) = log(factorial(r(n)));
            % term 2 is just the sum of the f's (should be uniform), and term 3 is independent of hdg, so can omit them

        % figure(7);clf; plot(hdg,f(n,:),0,r(n),'ro'); title(['n=' num2str(n) '  r=' num2str(r(n))]); pause;

    else

%-------------------------------------------------------------
%-------------------------------------------------------------
            % % simple Poisson likelihood (Foldiak, Sanger, Dayan & Abbott, Ma et al.):
        if strcmp(likemodel,'Poisson')
            
            if sum((g*f(n,:)).^r(n)) == Inf || factorial(r(n)) == Inf
                disp(['Inf problem! Cell ID #' num2str(n)]);
        %            r(n)
        %            f(n,:)
                f(n,:) = f(n,:) * 0.8;  % This kluges one cell whose high FR makes the equation exceed
                r(n) = round(r(n)*0.8); % the largest representable positive floating point number
            end
            if sum((g*f(n,:)).^r(n)) == Inf || factorial(r(n)) == Inf
                disp(['Inf problem after scale-down! Cell ID #' num2str(n)]);
                f(n,:) = f(n,:) * 0.5;  % One(?) cell needs even more of a scale-down
                r(n) = round(r(n)*0.5); 
            end
            if sum((g*f(n,:)).^r(n)) == Inf || factorial(r(n)) == Inf
                disp(['Inf problem after double scale-down! Cell ID #' num2str(n)]);
                keyboard
            end
            if sum(f(n,:)) == 0
                disp(['Zero problem: Cell ID #' num2str(n)]);
                f(n,:) = 1; % cell does not contribute, but does not zero out the likelihood either
            end

            temp_L(n,:) = ( exp(-g*f(n,:)) .* (g*f(n,:)).^r(n) ) / factorial(r(n));
%             figure(12); clf; plot(hdg,temp_L(n,:)); take a look at some, but beware uncommenting this for a full run

        else
%-------------------------------------------------------------
%-------------------------------------------------------------
            % ad-hoc likelihood (from empirical firing rate distribution):

            % First set bins to be the entire range of firing rates across all
            % headings; this way, a r(n) that is completely outside the actual 
            % rate distribution for a particular heading will assign a probability
            % of zero to that heading, rather than some intermediate value (as was
            % happening for flat but narrow (square-like) histos)

            all_trials = [];
            for h = 1:length(heading)
                all_trials = [all_trials resp_trial_all{n}{s,2,k,h}];
            end
            bin_sep = (max(all_trials)-min(all_trials)) / 20; % the denominator is arbitrary
            if min(all_trials)-bin_sep < 0
                bins = 0 : bin_sep : max(all_trials)+bin_sep;
            else
                bins = min(all_trials)-bin_sep : bin_sep : max(all_trials)+bin_sep;
            end

%-------------------------------------------------------------
            % histogram method:
            if strcmp(likemodel,'histogram')
%                 figure(1);plot(hdg,f(n,:));
        %         figure; set(gcf,'Position',[50,600,1500,300]);
                like = zeros(1,length(heading));
                for h = 1:length(heading)
        %             subplot(1,7,h);hist(resp_trial_all{n}{s,2,k,h},bins);
                    [counts,binctr] = hist(resp_trial_all{n}{s,2,k,h},bins);
        %             if sum(counts==0) > 0 % INTERPOLATE SPARSE HISTOS!
        %                 counts = interp1(binctr(counts>0),counts(counts>0),binctr,'cubic');
        % %                 figure; plot(binctr,counts);
        %             end
                    % no, don't interp, just set zeroes to near-zero
                    counts(counts==0) = 0.1;
                    countProb = counts/sum(counts);
                    findBestBin = find(abs(binctr-r(n)) == min(abs(binctr-r(n))));
                    like(h) = countProb(findBestBin(1));
                end
                
%-------------------------------------------------------------
%-------------------------------------------------------------
            else
                % instead of histo, try fitting gamma or normal dist
                if r(n)==0
                    r(n)=0.1;
                end
                for h = 1:length(heading)            
        %             phat = gamfit(resp_trial_all{n}{s,2,k,h});
        %             estPDF = gampdf(bins,phat(1),phat(2));

        %             [muhat,sigmahat] = normfit(resp_trial_all{n}{s,2,k,h});
        %             estPDF = normpdf(bins,muhat,sigmahat);

                    % don't need to fit again if passing in interpolated params
                    if strcmp(likemodel,'Gamma')
                        estPDF = gampdf(bins,param1_interp(n,h),param2_interp(n,h));
                    elseif strcmp(likemodel,'Gaussian')
                        estPDF = normpdf(bins,param1_interp(n,h),param2_interp(n,h));
                    end

        %             figure(291);subplot(1,7,h);plot(bins,estPDF);
                    estPDF(estPDF<=0) = min(estPDF(estPDF>0))/2; % kluge zero values to avoid zeroing out
                    estPDF(isinf(estPDF)) = max(estPDF(~isinf(estPDF)))*2; % kluge inf values to avoid NaNs later
                    findBestBin = find(abs(bins-r(n)) == min(abs(bins-r(n))));
                    like(h) = estPDF(findBestBin(1));

        %             % temp: to test variance/mean against Poisson:
        %             outfile = ['C:\Documents and Settings\Chris\My Documents\MATLAB\var_mean_vestib.dat'];
        %             fid = fopen(outfile, 'a'); buff = sprintf('%4.3f\t%4.3f\t', sigmahat^2, muhat);
        %             fprintf(fid, '%s', buff); fprintf(fid, '\r\n'); fclose(fid);
                end

%-------------------------------------------------------------
            end

            temp_L(n,:) = interp1(heading,like,hdg,'linear');
%             figure(1); clf; plot(hdg,temp_L(n,:)); title(num2str([s j k i])); pause;
            
            if sum(temp_L(n,:))==0 || isnan(sum(temp_L(n,:)))
                asdf;
            end
%-------------------------------------------------------------
%-------------------------------------------------------------
        end

    end

end

if s==3 && j~=2
    sum(remove==0) % for the relevant conditions, display how many cells were kept
end

temp_L(remove,:) = []; % remove skipped cells (nans)
                       % term1 and f do not need this because nansum works

if use_logL && strcmp(likemodel,'Poisson')
    % L = nansum(term1) - nansum(f) - nansum(term3);
    L = nansum(term1) - nansum(f);  % term 2 is just the sum of the f's (should be uniform), and term 3 is independent of hdg, so can omit them
else
    if N>1
        L = prod(temp_L);
    else
        L = temp_L; % for singlecells
    end
%     L2 = prod(temp_L2); % compare Poisson and nonPoisson
end


return