% population_likelihood.m 
% Trial-by-trial simulation of the behavioral experiment, computing
% likelihood functions by drawing spike counts from real neuronal data.

% THIS VERSION SIMULATES RESPONSES WITH A DESIRED NOISE CORRELATION LEVEL,
% ASSUMING GAUSSIAN NOISE WITH A VARIANCE = FANO*MEAN

close all
clear all
BASE_PATH = 'C:\Program Files\MATLAB\R2007a\work\';
% BASE_PATH = '/home/chris/work/'; % if running on linux cluster

% batch_name = 'all_cells_yos'
% batch_name = 'all_cells_web'
% batch_name = 'all_cells'

% batch_name = 'cong_yos_gt0'
% batch_name = 'cong_web_gt0'
% batch_name = 'cong_both_gt0'

% batch_name = 'cong_yos_gt0.4'
% batch_name = 'cong_web_gt0.4'
batch_name = 'cong_both_gt0.4'

% batch_name = 'rocable';
% batch_name = 'opp_both_0.4'

% batch_name = 'failed_F_test'
% batch_name = 'failed_AIC'

% batch_name = 'Fig7_black_bars';
% batch_name = 'weights_really_dont_change';

% batch_name = 'small_WdiffLow'
% batch_name = 'large_WdiffLow'
% batch_name = 'small_WdiffLow_v2'
% batch_name = 'large_WdiffLow_v2'

% these are using the weights constrained from 0 to 1.5
% batch_name = 'large_WdiffLow_v3'  % if I don't see over-weighting here, I quit.  (YES, thank goodness)
% batch_name = 'small_WdiffLow_v3'  % and this one nearly abolishes over-weighting (with under-weighting at high coh)

% 'final' versions here, using cong/blackbars -- should still work because negative weights are more rare in this sample
% batch_name = 'small_WdiffLow_v4'
% batch_name = 'large_WdiffLow_v4'  

% batch_name = 'singlecue_Tratio_favors_ves'
% batch_name = 'singlecue_Tratio_favors_vis'
% batch_name = 'singlecue_Tratio_balanced_v1'
% batch_name = 'singlecue_Tratio_balanced_v2'

% load([BASE_PATH 'condition_list_neuro.mat']);
load([BASE_PATH 'condition_list_interp.mat']);

if BASE_PATH(1)=='/'
    load([BASE_PATH 'neurons_indiv_trials/' batch_name '.mat']);
    files = textread([BASE_PATH 'neurons (temp)/' batch_name '.txt'],'%q');
else
    load([BASE_PATH 'neurons_indiv_trials\' batch_name '.mat']);
    files = textread([BASE_PATH 'neurons (temp)\' batch_name '.txt'],'%q');
end

num_reps = 60;
num_sim_repeats = 1;
save_data = 1;

use_logL = 0;
use_fake_data = 1;
    make_fake = 0;
subtract_baseline = 0;
    
noise_corr_max = 0;
% noise_corr_max = 0.1503; % trained
% noise_corr_max = 0.2683; % naive

% noisemodel = 'gausian'; fano = 1.5;
noisemodel = 'poisson'; fano = 1;

repetition = num_reps;
cell_num = length(resp_mean_all)
% method = 'spline';
method = 'linear';

wichmann = 0;
plotflag = 1; % for psychometrics
prior = 0;
Jan = 0;
Yong = 0;
FILE = ['SIM'];

stimtype = condition_list(:,1);
coherence = condition_list(:,2);
heading = condition_list(:,3);
delta = condition_list(:,4);

unique_stim_type = unique(stimtype);
unique_coherence = unique(coherence);
unique_heading = unique(heading);
unique_delta = unique(delta);
unique_amplitude = .30;
unique_azimuth = 90;
con_txt = 'delta';
con_txt_2 = 'coherence';

x = -10:0.1:10;
total_tr = length(stimtype) * num_reps;
stim = zeros(1,total_tr); coh = zeros(1,total_tr);
hdg = zeros(1,total_tr); delt= zeros(1,total_tr);
choice = zeros(1,total_tr);


%------------------------------------------------
% construct covariance matrix, separately for each stim type.
% method from Shadlen et al. 1996;  code modified from Yong's.

% first convert mean firing rates to Yong's format, and *interpolate*
% X = unique_heading';
X = [-10 -3.5 -1.3 0 1.3 3.5 10];
% XI = unique_heading(1):0.2:unique_heading(end);
% XI = [-10 -3.5:0.2:-1.3 -1.2:0.1:1.2 1.3:0.2:3.5 10];
XI = unique_heading';

if ~make_fake  % for now, only do all this if not making fake single-cue responses

for n = 1:cell_num 
    firing{1,2,1}(n,:) = interp1(X,resp_mean_all{n}{1,2,1},XI,method);
    firing{2,2,1}(n,:) = interp1(X,resp_mean_all{n}{2,2,1},XI,method);
    firing{2,2,2}(n,:) = interp1(X,resp_mean_all{n}{2,2,2},XI,method);
    firing{3,1,1}(n,:) = interp1(X,resp_mean_all{n}{3,1,1},XI,method);
    firing{3,2,1}(n,:) = interp1(X,resp_mean_all{n}{3,2,1},XI,method);
    firing{3,3,1}(n,:) = interp1(X,resp_mean_all{n}{3,3,1},XI,method);
    firing{3,1,2}(n,:) = interp1(X,resp_mean_all{n}{3,1,2},XI,method);
    firing{3,2,2}(n,:) = interp1(X,resp_mean_all{n}{3,2,2},XI,method);
    firing{3,3,2}(n,:) = interp1(X,resp_mean_all{n}{3,3,2},XI,method);
end

% % build covariance matrix r_noise to approximate a desired maximum noise
% % correlation, scaled by signal correlation
% for stim_type = 1:3   
%     if stim_type == 1
%         cohnum = 1;
%     else
%         cohnum = 2;
%     end
%     for i = 1:cell_num
%         for j = 1:cell_num
%             
%             %***** option 1: STIM-TYPE-SPECIFIC R_SIGNAL
% %             temp = corrcoef(resp_mean_all{i}{stim_type,2,cohnum}, resp_mean_all{j}{stim_type,2,cohnum});
%             %***** option 2: FIXED R_SIGNAL BASED ON (for example) VISUAL TUNING
% %             temp = corrcoef(resp_mean_all{i}{2,2,2}, resp_mean_all{j}{2,2,2});
%             %***** option 3: R_signal based on both vis and ves tuning, appended 
%             temp = corrcoef([resp_mean_all{i}{1,2,1} resp_mean_all{i}{2,2,2}], [resp_mean_all{j}{1,2,1} resp_mean_all{j}{2,2,2}]);            
%             r_signal{stim_type}(i,j) = temp(2,1);
%             if j==i
%                 r_noise{stim_type}(i,j) = 1; 
%             else
%                  r_noise{stim_type}(i,j) = noise_corr_max*exp( (r_signal{stim_type}(i,j)-1)/0.4785 );
%             end 
%         end
%     end
% end

for i = 1:cell_num
    for j = 1:cell_num
        %***** option 4: separate R_signal for vis and ves, but R_noise
        %(for all conditions) is based on mult. regress of both
        temp = corrcoef(resp_mean_all{i}{1,2,1}, resp_mean_all{j}{1,2,1});
        r_signal{1}(i,j) = temp(2,1);
        temp = corrcoef(resp_mean_all{i}{2,2,2}, resp_mean_all{j}{2,2,2});
        r_signal{2}(i,j) = temp(2,1);
        if j==i
            r_noise{1}(i,j) = 1;
            r_noise{2}(i,j) = 1;
            r_noise{3}(i,j) = 1;
        else
            r_noise{1}(i,j) = 0.01475 + 0.0855*r_signal{1}(i,j) + 0.07039*r_signal{2}(i,j); % these values are from Yong, probably trained animals -- DOUBLE CHECK IF USED AGAIN IN REVISION
            r_noise{2}(i,j) = r_noise{1}(i,j);
            r_noise{3}(i,j) = r_noise{1}(i,j);
            if noise_corr_max == 0
                r_noise{1}(i,j) = 0;
                r_noise{2}(i,j) = 0;
                r_noise{3}(i,j) = 0;
            end
        end 
    end
end

% now generate partially correlated responses according to r_noise
if noisemodel == 'gausian'  % (this method only works for Gaussian noise)
    
    s=1; j=2; k=1; % vestib
    rootQ = sqrtm(r_noise{s});
    for i=1:length(XI)
        for r = 1:repetition
            z = normrnd(0,1,cell_num,1);
            y = rootQ * z;
            y = y.*sqrt(fano*firing{s,j,k}(:,i)) + firing{s,j,k}(:,i); % variance is fano*mean, and SD is sqrt of that
            y(find(y<0)) = 0; % no negative reponses
            firing_trial{s,j,k,i}(:,r) = y;
        end
    end

    s=2; j=2; % visual
    for k = 1:2 % coherence
        rootQ = sqrtm(r_noise{s});
        for i=1:length(XI)
            for r = 1:repetition
                z = normrnd(0,1,cell_num,1);
                y = rootQ * z;
                y = y.*sqrt(fano*firing{s,j,k}(:,i)) + firing{s,j,k}(:,i); % variance is fano*mean, and SD is sqrt of that
                y(find(y<0)) = 0; % no negative reponses
                firing_trial{s,j,k,i}(:,r) = y;
            end
        end
    end

    s=3; % combined
    for j = 1:3 % conflict angle
        for k = 1:2 % coherence
            rootQ = sqrtm(r_noise{s});
            for i=1:length(XI)
                for r = 1:repetition
                    z = normrnd(0,1,cell_num,1);
                    y = rootQ * z;
                    y = y.*sqrt(fano*firing{s,j,k}(:,i)) + firing{s,j,k}(:,i); % variance is fano*mean, and SD is sqrt of that
                    y(find(y<0)) = 0; % no negative reponses
                    firing_trial{s,j,k,i}(:,r) = y;
                end
            end
        end
    end

else % for Poisson, just generate firing_trial (zero correlation) using poissrnd at the interpolated values
    
    for n = 1:cell_num 
        for s = 1:3; % stim type
            for j = 1:3 % conflict angle
                for k = 1:2 % coherence
                    if ~isempty(firing{s,j,k})
                        for i=1:length(XI)
                            firing_trial{s,j,k,i}(n,1:repetition) = poissrnd(firing{s,j,k}(n,i),1,repetition);
                        end
                    end
                end
            end
        end
    end
    
end

end


% ------------------------------------------------------------------
if use_fake_data
% Here we have the option to replace the real combined responses with 
% simulated ones (weighted sums of single-cue responses) to test hypotheses
% about what causes vestibular over-weighting in the model

% first need to interpolate single-cue tuning curves
% (even though we don't need the high resolution for this version of
% poplike -- we'll pick off the desired values later)

method = 'linear';
xx = round(unique_heading'*100)/100; % rounding to 1.23, but okay
xi = xx(1):0.01:xx(end);

if make_fake
%    [resp_mean_all, makefake_baselines] = make_fake_data_for_poplike;
%    [resp_mean_all, makefake_baselines] = make_fake_data_for_poplike_realistic;
	[resp_mean_all, makefake_baselines] = make_fake_data_for_poplike_Alex;
    cell_num = length(resp_mean_all);
    for n = 1:cell_num
        firing{1,2,1}(n,:) = interp1(X,resp_mean_all{n}{1,2,1},XI,method);
        firing{2,2,1}(n,:) = interp1(X,resp_mean_all{n}{2,2,1},XI,method);
        firing{2,2,2}(n,:) = interp1(X,resp_mean_all{n}{2,2,2},XI,method);
%         firing{3,1,1}(n,:) = interp1(X,resp_mean_all{n}{3,1,1},XI,method);
%         firing{3,2,1}(n,:) = interp1(X,resp_mean_all{n}{3,2,1},XI,method);
%         firing{3,3,1}(n,:) = interp1(X,resp_mean_all{n}{3,3,1},XI,method);
%         firing{3,1,2}(n,:) = interp1(X,resp_mean_all{n}{3,1,2},XI,method);
%         firing{3,2,2}(n,:) = interp1(X,resp_mean_all{n}{3,2,2},XI,method);
%         firing{3,3,2}(n,:) = interp1(X,resp_mean_all{n}{3,3,2},XI,method);
        for s = 1:2; % stim type
            for j = 2 % conflict angle
                for k = 1:2 % coherence
                    if ~isempty(firing{s,j,k})
                        for i=1:length(XI)
                            firing_trial{s,j,k,i}(n,1:repetition) = poissrnd(firing{s,j,k}(n,i),1,repetition);
                        end
                    end
                end
            end
        end
    end
end

cell_num = length(resp_mean_all);
repetition = num_reps;

% Wves_range = 0.32; % 0.63/0.32 % low/high coh
% Wvis_range = 0.86; % 0.57/0.86

load weights.mat  %  NOTE: TO RETURN TO THIS STRATEGY, MUST REMAKE WEIGHTS.MAT (current one is matched to congruent cell list, not the 59 "sig boths"
veslow = weights(:,1); vislow = weights(:,2);
veshigh = weights(:,3); vishigh = weights(:,4);
disp('generating fake data with specified weights');
for p = 1 %:length(Wves_range)
    p
    for q = 1 %:length(Wvis_range)
        q
        x_low = [0.63 0.57];
        x_high = [0.32 0.86];
        
        for n = 1:cell_num
    
            % instead of means, sample Wves and Wvis from Gaussian distributions with
            % same mean/SD -- SD's are: 0.56/0.61 (ves) and 0.76/0.46 (vis)

%             x_low(1) = randn*0.56+0.63;
%             x_high(1) = randn*0.61+0.32;
%             x_low(2) = randn*0.76+0.57;
%             x_high(2) = randn*0.46+0.86;
                        
%             no, need to sample with covariance structure too -- if wves is low/negative, wvis must be high/positive
%             define conditional probabilities? at least for binned
%             segments of the original dist...
            
            % forget about sampling, just recapitulate entire distribution, cell by cell
%             x_low(1) = veslow(n); x_low(2) = vislow(n);
%             x_high(1) = veshigh(n); x_high(2) = vishigh(n);
            % quick sanity check:
%             x_low(1) = 1; x_low(2) = 1;
%             x_high(1) = 1; x_high(2) = 1;
            
%             % sample vestib first, then vis conditionally
%             x_low(1) = randn*std(veslow) + mean(veslow);
%             if x_low(1)<0
%                 x_low(2) = randn*std(vislow(veslow<0)) + mean(vislow(veslow<0));
%             elseif x_low(1)>0 && x_low(1)<0.5
%                 x_low(2) = randn*std(vislow(veslow>0 & veslow<0.5)) + mean(vislow(veslow>0 & veslow<0.5));                
%             elseif x_low(1)>0.5 && x_low(1)<1
%                 x_low(2) = randn*std(vislow(veslow>0.5 & veslow<1)) + mean(vislow(veslow>0.5 & veslow<1));                
%             else % x_low(1)>1
%                 x_low(2) = randn*std(vislow(veslow>1)) + mean(vislow(veslow>1));
%             end
%             x_high(1) = randn*std(veshigh) + mean(veshigh);
%             if x_high(1)<0
%                 x_high(2) = randn*std(vishigh(veshigh<0)) + mean(vishigh(veshigh<0));
%             elseif x_high(1)>0 && x_high(1)<0.5
%                 x_high(2) = randn*std(vishigh(veshigh>0 & veshigh<0.5)) + mean(vishigh(veshigh>0 & veshigh<0.5));                
%             elseif x_high(1)>0.5 && x_high(1)<1
%                 x_high(2) = randn*std(vishigh(veshigh>0.5 & veshigh<1)) + mean(vishigh(veshigh>0.5 & veshigh<1));                
%             else % x_high(1)>1
%                 x_high(2) = randn*std(vishigh(veshigh>1)) + mean(vishigh(veshigh>1));
%             end

%             % sample vis first, then vestib conditionally
%             x_low(2) = randn*std(vislow) + mean(vislow);
%             if x_low(2)<0
%                 x_low(1) = randn*std(veslow(vislow<0)) + mean(veslow(vislow<0));
%             elseif x_low(2)>0 && x_low(2)<0.5
%                 x_low(1) = randn*std(veslow(vislow>0 & vislow<0.5)) + mean(veslow(vislow>0 & vislow<0.5));                
%             elseif x_low(2)>0.5 && x_low(2)<1
%                 x_low(1) = randn*std(veslow(vislow>0.5 & vislow<1)) + mean(veslow(vislow>0.5 & vislow<1));                
%             else % x_low(2)>1
%                 x_low(1) = randn*std(veslow(vislow>1)) + mean(veslow(vislow>1));
%             end
%             x_high(2) = randn*std(vishigh) + mean(vishigh);
%             if x_high(2)<0
%                 x_high(1) = randn*std(veshigh(vishigh<0)) + mean(veshigh(vishigh<0));
%             elseif x_high(2)>0 && x_high(2)<0.5
%                 x_high(1) = randn*std(veshigh(vishigh>0 & vishigh<0.5)) + mean(veshigh(vishigh>0 & vishigh<0.5));                
%             elseif x_high(2)>0.5 && x_high(2)<1
%                 x_high(1) = randn*std(veshigh(vishigh>0.5 & vishigh<1)) + mean(veshigh(vishigh>0.5 & vishigh<1));                
%             else % x_high(2)>1
%                 x_high(1) = randn*std(veshigh(vishigh>1)) + mean(veshigh(vishigh>1));
%             end

            if make_fake
                file = ['makefake' num2str(n)];
                baseline_rate = makefake_baselines(n); % from make_fake function
            else
                file = files{n};
                if BASE_PATH(1)=='/'
                    matfile = [BASE_PATH 'neurons_indiv_trials/' file '.mat'];
                else
                    matfile = [BASE_PATH 'neurons_indiv_trials\' file '.mat'];
                end
                load(matfile, 'baseline_rate');
            end
%             if n == 1 || n == 13 % temp-kluge a couple problem cells
%                 baseline_rate = baseline_rate/2;
%             end
                
            
            yi_vislow{n} = interp1(xx,firing{2,2,1}(n,:)-subtract_baseline*baseline_rate,xi,method);
            yi_vishigh{n} = interp1(xx,firing{2,2,2}(n,:)-subtract_baseline*baseline_rate,xi,method);
            yi_ves{n} = interp1(xx,firing{1,2,1}(n,:)-subtract_baseline*baseline_rate,xi,method);
            % figure;plot(xx,resp_mean_all{n}{2,2,1},'bo',xx,resp_mean_all{n}{2,2,2},'co',xx,resp_mean_all{n}{1,2,1},'ro');
            % hold on; plot(xi,yi_vislow{n},'b-',xi,yi_vishigh{n},'c-',xi,yi_ves{n},'r-');

            % also need to extrapolate to +/- delta/2 degrees outside largest heading
            x_chunk = xi(1 : find(xi==xx(2))-1);
            y_chunk = yi_vislow{n}(1 : find(xi==xx(2))-1);
            P = polyfit(x_chunk,y_chunk,1);
            xi_seg1 = xx(1)-4/2 : 0.01 : xx(1)-0.01;
            yi_seg1 = P(1)*xi_seg1 + P(2);
            x_chunk = xi(find(xi==xx(end-1)) : end);
            y_chunk = yi_vislow{n}(find(xi==xx(end-1)) : end);
            P = polyfit(x_chunk,y_chunk,1);
            xi_seg2 = xx(end)+0.01 : 0.01 : xx(end)+4/2;
            yi_seg2 = P(1)*xi_seg2 + P(2);
            yi_vislow{n} = [yi_seg1 yi_vislow{n} yi_seg2];

            x_chunk = xi(1 : find(xi==xx(2))-1);
            y_chunk = yi_vishigh{n}(1 : find(xi==xx(2))-1);
            P = polyfit(x_chunk,y_chunk,1);
            xi_seg1 = xx(1)-4/2 : 0.01 : xx(1)-0.01;
            yi_seg1 = P(1)*xi_seg1 + P(2);
            x_chunk = xi(find(xi==xx(end-1)) : end);
            y_chunk = yi_vishigh{n}(find(xi==xx(end-1)) : end);
            P = polyfit(x_chunk,y_chunk,1);
            xi_seg2 = xx(end)+0.01 : 0.01 : xx(end)+4/2;
            yi_seg2 = P(1)*xi_seg2 + P(2);
            yi_vishigh{n} = [yi_seg1 yi_vishigh{n} yi_seg2];

            x_chunk = xi(1 : find(xi==xx(2))-1);
            y_chunk = yi_ves{n}(1 : find(xi==xx(2))-1);
            P = polyfit(x_chunk,y_chunk,1);
            xi_seg1 = xx(1)-4/2 : 0.01 : xx(1)-0.01;
            yi_seg1 = P(1)*xi_seg1 + P(2);
            x_chunk = xi(find(xi==xx(end-1)) : end);
            y_chunk = yi_ves{n}(find(xi==xx(end-1)) : end);
            P = polyfit(x_chunk,y_chunk,1);
            xi_seg2 = xx(end)+0.01 : 0.01 : xx(end)+4/2;
            yi_seg2 = P(1)*xi_seg2 + P(2);
            yi_ves{n} = [yi_seg1 yi_ves{n} yi_seg2];

            xi_new = xx(1)-4/2 : 0.01 : xx(end)+4/2;
            % figure;plot(xx,resp_mean_all{n}{2,2,1},'bo',xx,resp_mean_all{n}{2,2,2},'co',xx,resp_mean_all{n}{1,2,1},'ro');
            % hold on; plot(xi_new,yi_vislow{n},'b-',xi_new,yi_vishigh{n},'c-',xi_new,yi_ves{n},'r-');

            % now simulate conflict conditions with the best-fitting linear weights
            delta_shift = 4/2/0.01;
            xi_new = round(xi_new*100)/100; % stupid rounding errors
            for h = 1:length(xx)
                comblowminus(h) = x_low(1)*yi_ves{n}(find(xi_new==xx(h))+delta_shift) + x_low(2)*yi_vislow{n}(find(xi_new==xx(h))-delta_shift) + subtract_baseline*baseline_rate;
                comblowzero(h) =  x_low(1)*yi_ves{n}(find(xi_new==xx(h))) + x_low(2)*yi_vislow{n}(find(xi_new==xx(h))) + subtract_baseline*baseline_rate;
                comblowplus(h) =  x_low(1)*yi_ves{n}(find(xi_new==xx(h))-delta_shift) + x_low(2)*yi_vislow{n}(find(xi_new==xx(h))+delta_shift) + subtract_baseline*baseline_rate;

                combhighminus(h) = x_high(1)*yi_ves{n}(find(xi_new==xx(h))+delta_shift) + x_high(2)*yi_vishigh{n}(find(xi_new==xx(h))-delta_shift) + subtract_baseline*baseline_rate;
                combhighzero(h) =  x_high(1)*yi_ves{n}(find(xi_new==xx(h))) + x_high(2)*yi_vishigh{n}(find(xi_new==xx(h))) + subtract_baseline*baseline_rate;
                combhighplus(h) =  x_high(1)*yi_ves{n}(find(xi_new==xx(h))-delta_shift) + x_high(2)*yi_vishigh{n}(find(xi_new==xx(h))+delta_shift) + subtract_baseline*baseline_rate;
            end

            % extrapolation can cause negative values, leading to NaNs from poissrnd below
            comblowminus(comblowminus<0) = 0;
            comblowzero(comblowzero<0) = 0;
            comblowplus(comblowplus<0) = 0;
            combhighminus(combhighminus<0) = 0;
            combhighzero(combhighzero<0) = 0;
            combhighplus(combhighplus<0) = 0;

%             figure(1); subplot(2,1,1); plot(X,resp_mean_all{n}{3,1,1},'bo-',X,resp_mean_all{n}{3,2,1},'g^-',X,resp_mean_all{n}{3,3,1},'rs-');
%             set(gca,'XTickMode','auto'); set(gca,'YTickMode','auto');
%             xlabel('Heading Angle'); ylabel('Firing Rate - Low Coh'); xlim([-15 15]); title(file);
%             yy = ylim;
%             
%             figure(1); subplot(2,1,2); plot(X,resp_mean_all{n}{3,1,2},'bo-',X,resp_mean_all{n}{3,2,2},'g^-',X,resp_mean_all{n}{3,3,2},'rs-');
%             set(gca,'XTickMode','auto'); set(gca,'YTickMode','auto');
%             xlabel('Heading Angle'); ylabel('Firing Rate - High Coh'); xlim([-15 15]); title('Real Data');
%             yy = ylim;
%             
%             figure(2); subplot(2,1,1); plot(xx,comblowminus,'bo-',xx,comblowzero,'g^-',xx,comblowplus,'rs-');
%             set(gca,'XTickMode','auto'); set(gca,'YTickMode','auto');
%             xlabel('Heading Angle'); ylabel('Firing Rate - Low Coh'); xlim([-15 15]); title(file);
%             yy = ylim;
%             text(-2,0.97*(yy(2)-yy(1))+yy(1),['Wves = ' num2str(x_low(1))]);
%             text(-2,0.93*(yy(2)-yy(1))+yy(1),['Wvis = ' num2str(x_low(2))]);
% %             text(-2,0.89*(yy(2)-yy(1))+yy(1),['DC = ' num2str(x_low(3))]);
%             
%             figure(2); subplot(2,1,2); plot(xx,combhighminus,'bo-',xx,combhighzero,'g^-',xx,combhighplus,'rs-');
%             set(gca,'XTickMode','auto'); set(gca,'YTickMode','auto');
%             xlabel('Heading Angle'); ylabel('Firing Rate - High Coh'); xlim([-15 15]); title('Linear Model');
%             yy = ylim;
%             text(-2,0.97*(yy(2)-yy(1))+yy(1),['Wves = ' num2str(x_high(1))]);
%             text(-2,0.93*(yy(2)-yy(1))+yy(1),['Wvis = ' num2str(x_high(2))]);
% %             text(-2,0.89*(yy(2)-yy(1))+yy(1),['DC = ' num2str(x_high(3))]);


            % lastly, generate single-trial responses for decoding
            % s=3; % combined
            % for j = 1:3 % conflict angle
            %     for k = 1:2 % coherence
            %         rootQ = sqrtm(r_noise{s});
            %         for i=1:length(XI)
            %             for r = 1:repetition
            %                 z = normrnd(0,1,cell_num,1);
            %                 y = rootQ * z;
            %                 y = y.*sqrt(fano*firing{s,j,k}(:,i)) + firing{s,j,k}(:,i); % variance is fano*mean, and SD is sqrt of that
            %                 y(find(y<0)) = 0; % no negative reponses
            %                 firing_trial{s,j,k,i}(:,r) = y;
            %             end
            %         end
            %     end
            % end

            
            % store fake data with desired weights (comb), and original 
            % singles, as resp_mean_weightmap
            resp_mean_weightmap{n,p,q}{1,2,1} = firing{1,2,1}(n,:);
            resp_mean_weightmap{n,p,q}{2,2,1} = firing{2,2,1}(n,:);
            resp_mean_weightmap{n,p,q}{2,2,2} = firing{2,2,2}(n,:);
    
            resp_mean_weightmap{n,p,q}{3,1,1} = comblowminus;
            resp_mean_weightmap{n,p,q}{3,2,1} = comblowzero;
            resp_mean_weightmap{n,p,q}{3,3,1} = comblowplus;

            resp_mean_weightmap{n,p,q}{3,1,2} = combhighminus;
            resp_mean_weightmap{n,p,q}{3,2,2} = combhighzero;
            resp_mean_weightmap{n,p,q}{3,3,2} = combhighplus;

            for i = 1:length(xx)  
                resp_trial_weightmap{n,p,q}{1,2,1,i} = firing_trial{1,2,1,i}(n,:);
                resp_trial_weightmap{n,p,q}{2,2,1,i} = firing_trial{2,2,1,i}(n,:);
                resp_trial_weightmap{n,p,q}{2,2,2,i} = firing_trial{2,2,2,i}(n,:);    
                if noisemodel == 'poisson'
                    resp_trial_weightmap{n,p,q}{3,1,1,i}(1:repetition) = poissrnd(comblowminus(i),1,repetition);
                    resp_trial_weightmap{n,p,q}{3,2,1,i}(1:repetition) = poissrnd(comblowzero(i),1,repetition);
                    resp_trial_weightmap{n,p,q}{3,3,1,i}(1:repetition) = poissrnd(comblowplus(i),1,repetition);

                    resp_trial_weightmap{n,p,q}{3,1,2,i}(1:repetition) = poissrnd(combhighminus(i),1,repetition);
                    resp_trial_weightmap{n,p,q}{3,2,2,i}(1:repetition) = poissrnd(combhighzero(i),1,repetition);
                    resp_trial_weightmap{n,p,q}{3,3,2,i}(1:repetition) = poissrnd(combhighplus(i),1,repetition);
                else % Gaussian
                    resp_trial_weightmap{n,p,q}{3,1,1,i}(1:repetition) = comblowminus(i) + randn(1,repetition)*sqrt(fano*comblowminus(i));
                    resp_trial_weightmap{n,p,q}{3,2,1,i}(1:repetition) = comblowzero(i) + randn(1,repetition)*sqrt(fano*comblowzero(i));
                    resp_trial_weightmap{n,p,q}{3,3,1,i}(1:repetition) = comblowplus(i) + randn(1,repetition)*sqrt(fano*comblowplus(i));

                    resp_trial_weightmap{n,p,q}{3,1,2,i}(1:repetition) = combhighminus(i) + randn(1,repetition)*sqrt(fano*combhighminus(i));
                    resp_trial_weightmap{n,p,q}{3,2,2,i}(1:repetition) = combhighzero(i) + randn(1,repetition)*sqrt(fano*combhighzero(i));
                    resp_trial_weightmap{n,p,q}{3,3,2,i}(1:repetition) = combhighplus(i) + randn(1,repetition)*sqrt(fano*combhighplus(i));
                end
            end
            
%             if sum([comblowminus comblowzero comblowplus combhighminus combhighzero combhighplus resp_mean_weightmap{n}{1,2,1} resp_mean_all{n}{2,2,1} resp_mean_all{n}{2,2,2}] < 0)
%                 disp(['ERROR -- negative mean FR in cellnum ' num2str(n)]);
%                 p
%                 q
%             end

%         pause
%         figure(1); clf;
%         figure(2); clf;
        end
        
    end
end    

% ------------------------------------------------------------------
else
    % now convert both means and single-trial responses back to my format
    clear resp_trial_all resp_mean_all  % not sure if I have to clear these
    for n = 1:cell_num % convert mean firing rates to Yong's format
        
        resp_mean_all{n}{1,2,1} = firing{1,2,1}(n,:);
        resp_mean_all{n}{2,2,1} = firing{2,2,1}(n,:);
        resp_mean_all{n}{2,2,2} = firing{2,2,2}(n,:);
        resp_mean_all{n}{3,1,1} = firing{3,1,1}(n,:);
        resp_mean_all{n}{3,2,1} = firing{3,2,1}(n,:);
        resp_mean_all{n}{3,3,1} = firing{3,3,1}(n,:);
        resp_mean_all{n}{3,1,2} = firing{3,1,2}(n,:);
        resp_mean_all{n}{3,2,2} = firing{3,2,2}(n,:);
        resp_mean_all{n}{3,3,2} = firing{3,3,2}(n,:);

        for i = 1:length(XI)
            resp_trial_all{n}{1,2,1,i} = firing_trial{1,2,1,i}(n,:);
            resp_trial_all{n}{2,2,1,i} = firing_trial{2,2,1,i}(n,:);
            resp_trial_all{n}{2,2,2,i} = firing_trial{2,2,2,i}(n,:);        
            resp_trial_all{n}{3,1,1,i} = firing_trial{3,1,1,i}(n,:);
            resp_trial_all{n}{3,2,1,i} = firing_trial{3,2,1,i}(n,:);
            resp_trial_all{n}{3,3,1,i} = firing_trial{3,3,1,i}(n,:);
            resp_trial_all{n}{3,1,2,i} = firing_trial{3,1,2,i}(n,:);
            resp_trial_all{n}{3,2,2,i} = firing_trial{3,2,2,i}(n,:);
            resp_trial_all{n}{3,3,2,i} = firing_trial{3,3,2,i}(n,:);
        end
        
    end
    
end % end use_fake_data

% % this is getting ridiculous.  reshape data for faster indexing later
% for p = 1 %:length(Wves_range)
%     for q = 1 %:length(Wvis_range)
%         for s=1:3
%             for j=1:3
%                 for k=1:2
%                     if ~isempty(resp_mean_weightmap{n,p,q}{s,j,k})
%                         for i=1:length(xx)
%                             for n = 1:cell_num
%                                 resp_mean_weightmap2{p,q}(n,s,j,k,i) = resp_mean_weightmap{n,p,q}{s,j,k}(i);
%                                 resp_trial_weightmap2{p,q}(n,s,j,k,i,1:repetition) = resp_trial_weightmap{n,p,q}{s,j,k,i};
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end

sig = 5;
Pr = normpdf(x,0,5);
Pr = Pr / sum(Pr);

%------------------------------------------------
if use_fake_data
% disp('computing likelihoods and generating choice distributions');
% for p = 1 %:length(Wves_range)
% for q = 1 %:length(Wvis_range)
%     disp(['p=' num2str(p) ', q=' num2str(q)]);
%     clear resp_mean_all_temp resp_trial_all_temp
% 	  resp_mean_all_temp = resp_mean_weightmap2{p,q};
%     resp_trial_all_temp = resp_trial_weightmap2{p,q};
for n = 1:cell_num
    resp_mean_all{n} = resp_mean_weightmap{n,1,1};
    resp_trial_all{n} = resp_trial_weightmap{n,1,1};
end
end

for repeat_sim = 1:num_sim_repeats
tic

%------------------------------------------------
for n = 1:num_reps
    repeat_sim
    n
    for t = 1:length(condition_list)
        trial = (n-1)*length(condition_list)+t;
        stim(trial) = stimtype(t);
        coh(trial) = coherence(t);
        hdg(trial) = heading(t);
        delt(trial) = delta(t);
        
%                 hdg(trial) = 1.2250; delt(trial) = 0;
        
        if unique_stim_type == 1
            deltnum = 2;
        else
            deltnum = find(unique(delta)==delt(trial));
        end
        cohnum = find(unique(coherence)==coh(trial));
        hdgnum = find(unique(heading)==hdg(trial));

% %         since debug mode is not working:
        s=stim(trial);j=deltnum;k=cohnum;i=hdgnum;
        
%         [r,L] = compute_likelihood_corr(unique_heading',resp_trial_all,resp_mean_all,stim(trial),deltnum,cohnum,hdgnum,use_logL);
        [r,L] = compute_likelihood_backup(unique_heading',resp_trial_all,resp_mean_all,stim(trial),deltnum,cohnum,hdgnum,use_logL);
%         choice(trial) = likelihood_ratio_test(unique_heading',resp_trial_all,resp_mean_all,stim(trial),deltnum,cohnum,hdgnum);

        %         % repeat the same stimulus a bunch of times and compute average pairwise correlation
%         for q = 1:10
%             [Rtemp(q,:),Ltemp] = compute_likelihood(resp_trial_all,resp_mean_all,stim(trial),deltnum,cohnum,hdgnum);
%         end
%         corr_mat = corrcoef(Rtemp);
%         nan_these = logical(diag(ones(length(corr_mat),1),0));
%         corr_mat(nan_these) = NaN;
%         mean_corr_actual(t) = mean(nanmean(corr_mat));
%         corr_mat = r_noise{s};
%         nan_these = logical(diag(ones(length(corr_mat),1),0));
%         corr_mat(nan_these) = NaN;
%         mean_corr_desired(t) = mean(nanmean(corr_mat));
        
        L = L/sum(L); % normalize
        Post = L; % no prior for now

        if sum(Post)>0
            Post = Post/sum(Post); % normalize
            if Jan % alternate decision rule: compare area under posterior <0 vs. >0
                choice(trial) = (sum(Post(find(x==0)+1:end)) > sum(Post(1:find(x==0)-1))) * 2 - 1;
            elseif Yong % Yong's idea is to compare L with that of a reference 'trial' at zero heading, instead of a noiseless reference of zero
                [rzero,Lzero] = compute_likelihood_backup(unique_heading',resp_trial_all,resp_mean_all,stim(trial),2,cohnum,find(unique(heading)==0),use_logL);
                Lzero = Lzero.*Pr; Lzero = Lzero/sum(Lzero);
                findmap = find(Post==max(Post));
                findref = find(Lzero==max(Lzero));
                choice(trial) = sign(x(findmap(1)) - x(findref(1)));
            else % MAP or MLE decision rule
                findmap = find(Post==max(Post));
                MAP = x(findmap(1));
                choice(trial) = sign(MAP);
            end
        else
            disp('Likelihood is all zero -- WTF?');
            choice(trial) = sign(randn);
        end
    end
%     nanmean(mean_corr_actual)
%     nanmean(mean_corr_desired)
end % end num_reps

stim_type = stim; heading = hdg; coherence = coh; delta = delt;

RIGHT = 1;
LEFT = -1;
right_pct = []; correct_pct = [];
for s = 1:length(unique_stim_type)
	for j = 1:length(unique_delta)
        for k = 1:length(unique_coherence)
            for i = 1:length(unique_heading)
                trials_select = logical( (stim_type == unique_stim_type(s)) & (heading == unique_heading(i)) & (delta == unique_delta(j)) & (coherence == unique_coherence(k)) );
                if sum(trials_select) == 0
                    right_pct{s,j,k}(i) = NaN;
                    correct_pct{s,j,k}(i) = NaN;
                elseif sum(trials_select) ~= num_reps
                    disp('ERROR: num_reps');
                    return;
                else
                    right_trials = (trials_select & (choice == RIGHT));
                    right_pct{s,j,k}(i) = sum(right_trials) / sum(trials_select);
                    if unique_heading(i) > 0
                        correct_pct{s,j,k}(i) = right_pct{s,j,k}(i);
                    else
                        correct_pct{s,j,k}(i) = 1 - right_pct{s,j,k}(i);
                    end
                end
                fit_data_psycho_cum{s,j,k,repeat_sim}(i,1) = unique_heading(i);
                fit_data_psycho_cum{s,j,k,repeat_sim}(i,2) = right_pct{s,j,k}(i);
                fit_data_psycho_cum{s,j,k,repeat_sim}(i,3) = sum(trials_select);
            end
        end
	end
end

%%%%%% use either pfit or cum_gaussfit_max1 to estimate threshold and bias
for s = 1:length(unique_stim_type)
    for j = 1:length(unique_delta)
        for k = 1:length(unique_coherence)
            if fit_data_psycho_cum{s,j,k,repeat_sim}(1,3) == 0 % identifies and skips invalid condition combinations (e.g., vestibular only with a nonzero conflict angle)
                Thresh_psy{s,j,k} = NaN;
                Bias_psy{s,j,k} = NaN;
                psy_perf{s,j,k} = [NaN , NaN];
            else
                if wichmann
                    wichman_psy = pfit(fit_data_psycho_cum{s,j,k,repeat_sim},'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'sens',0,'compute_stats','false','verbose','false');
                    Bias_psy{s,j,k} = wichman_psy.params.est(1);
                    Thresh_psy{s,j,k} = wichman_psy.params.est(2);
                    psy_perf{s,j,k} = [wichman_psy.params.est(1),wichman_psy.params.est(2)];
                else
                    [bbb,ttt] = cum_gaussfit_max1(fit_data_psycho_cum{s,j,k,repeat_sim});
                    Bias_psy{s,j,k} = bbb;
                    Thresh_psy{s,j,k} = ttt;
                    psy_perf{s,j,k} = [bbb,ttt];
                end                
            end
        end
    end
end

%--------------------------------------------------------------------------
% compute the predicted and actual weights/thresholds
if length(unique_stim_type) > 2
    for k = 1:length(unique_coherence)
        Wves_actual_minus(k) = ((Bias_psy{3,1,k} - Bias_psy{3,2,k}) - (-unique_delta(1)/2)) / unique_delta(1);
        Wves_actual_plus(k) = ((Bias_psy{3,3,k} - Bias_psy{3,2,k}) - (-unique_delta(end)/2)) / unique_delta(end);
        Wves_actual(k) = (Wves_actual_minus(k) + Wves_actual_plus(k)) / 2;
        Wves_predicted(k) = Thresh_psy{2,2,k}^2/(Thresh_psy{1,2,1}^2+Thresh_psy{2,2,k}^2);
        thresh_predicted(k) = sqrt((Thresh_psy{2,2,k}^2*Thresh_psy{1,2,1}^2)/(Thresh_psy{2,2,k}^2+Thresh_psy{1,2,1}^2));
        thresh_actual_d0(k) = Thresh_psy{3,2,k};
        thresh_actual_all(k) = (Thresh_psy{3,1,k}+Thresh_psy{3,2,k}+Thresh_psy{3,3,k}) / 3;
        thresh_actual_pm(k) = (Thresh_psy{3,1,k}+Thresh_psy{3,3,k}) / 2;
        bias_delta0(k) = Bias_psy{3,2,k};
    end
else
    for k = 1:length(unique_coherence)
        Wves_actual_minus(k) = NaN;
        Wves_actual_plus(k) = NaN;
        Wves_actual(k) = NaN;
        Wves_predicted(k) = NaN;
        thresh_predicted(k) = NaN;
        thresh_actual_d0(k) = NaN;
        thresh_actual_all(k) = NaN;
        thresh_actual_pm(k) = NaN;
        bias_delta0(k) = NaN;
    end
end
        
if plotflag
%--------------------------------------------------------------------------
% plot psychometric function here
H{1,1} = 'bo'; H{2,1} = 'b^'; H{3,1} = 'bs'; H{4,1} = 'b*'; F{1} = 'b-';
H{1,2} = 'go'; H{2,2} = 'g^'; H{3,2} = 'gs'; H{4,2} = 'g*'; F{2} = 'g-';
H{1,3} = 'ro'; H{2,3} = 'r^'; H{3,3} = 'rs'; H{4,3} = 'r*'; F{3} = 'r-';
H{1,4} = 'co'; H{2,4} = 'c^'; H{3,4} = 'cs'; H{4,4} = 'c*'; F{4} = 'c-';
H{1,5} = 'mo'; H{2,5} = 'm^'; H{3,5} = 'ms'; H{4,5} = 'm*'; F{5} = 'm-';
H{1,6} = 'yo'; H{2,6} = 'y^'; H{3,6} = 'ys'; H{4,6} = 'y*'; F{6} = 'y-';
H{1,7} = 'ko'; H{2,7} = 'k^'; H{3,7} = 'ks'; H{4,7} = 'k*'; F{7} = 'k-';
t=prior;
%first the single-cues
for s = 1:2
    figure(2+1+10*t);
    set(2+1+10*t,'Position', [200,50 700,600], 'Name', 'Heading Discrimination');
    if s == 1
        axes('position',[0.2,0.25, 0.6,0.5] );
    end
    % fit data with cumulative gaussian and plot both raw data and fitted curve
    legend_txt = [];

    % x = min(unique_heading) : 0.1 : max(unique_heading);
    % instead, force x range to be symmetric about zero (for staircase)
%     x = -max(abs(unique_heading)) : 0.1 : max(abs(unique_heading));

    for j = 1:length(unique_delta)    % <-- currently conflict_angle
        for k = 1:length(unique_coherence)  % <-- currently coherence
            figure(2+1+10*t);
            plot(unique_heading, fit_data_psycho_cum{s,j,k,repeat_sim}(:,2), H{k,s}, x, cum_gaussfit(psy_perf{s,j,k}, x), F{s});
            set(gca,'XTickMode','auto'); set(gca,'YTickMode','auto');
            xlabel('Heading Angle');   
            ylim([0,1]);
            ylabel('Percent Rightward Choices');
            hold on;
            legend_txt{j*2-1} = [num2str(unique_delta(j))];
            legend_txt{j*2} = [''];
        end
    end
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output some text of basic parameters in the figure
axes('position',[0.2,0.8, 0.6,0.15] );
xlim( [0,50] );
ylim( [2,10] );
text(0, 10, FILE);
text(15,11, ['amplitude = ' num2str(unique_amplitude)]);
text(15,10, ['stimtype = 1+2']);
text(30,10, ['azimuth = ' num2str(mean(unique_azimuth))]); % mean, because actual 'AZIMUTH' varies with conflict angle
text(45,10, ['repeats = ' num2str(num_reps)]);
text(0,8.3, 'stimtype');
text(10,8.3, con_txt_2);
text(20,8.3, 'bias');
text(30,8.3, 'thresh');
text(40,8.3, '%correct');
text(50,8.3, 'Wves-pred');
for s = 1:2
    if s == 1
        text(0, 6, num2str(s));     
        text(10,6,['N/A']);
        text(20,6, num2str(Bias_psy{s,find(unique_delta==0),1}) );
        text(30,6, num2str(Thresh_psy{s,find(unique_delta==0),1}) );
        text(40,6, num2str(mean(mean([1-fit_data_psycho_cum{s,find(unique_delta==0),1}(1:3,2) fit_data_psycho_cum{s,find(unique_delta==0),1}(5:7,2)]))) );
    else
        for k = 1:length(unique_coherence)  % <-- currently coherence
            text(0, -1.5*k+5.5, num2str(s));
            text(10,-1.5*k+5.5, num2str(unique_coherence(k)));
            text(20,-1.5*k+5.5, num2str(Bias_psy{s,find(unique_delta==0),k}) );
            text(30,-1.5*k+5.5, num2str(Thresh_psy{s,find(unique_delta==0),k}) );
            text(40,-1.5*k+5.5, num2str( num2str(mean(mean([1-fit_data_psycho_cum{s,find(unique_delta==0),k}(1:3,2) fit_data_psycho_cum{s,find(unique_delta==0),k}(5:7,2)]))) ));
            text(50,-1.5*k+5.5, num2str(Wves_predicted(k)) );
        end
    end
end % end single cues

axis off;
% print;
% close;

% lastly, the combined
s = 3;
for k = 1:length(unique_coherence)
    figure(s+k+10*t);
    set(s+k+10*t,'Position', [200,50 700,600], 'Name', 'Heading Discrimination');
    axes('position',[0.2,0.25, 0.6,0.5] );
    % fit data with cumulative gaussian and plot both raw data and fitted curve
    legend_txt = [];

    % x = min(unique_heading) : 0.1 : max(unique_heading);
    % instead, force x range to be symmetric about zero (for staircase)
%     x = -max(abs(unique_heading)) : 0.1 : max(abs(unique_heading));

    for j = 1:length(unique_delta)    % <-- currently conflict_angle
        figure(s+k+10*t);
        plot(unique_heading, fit_data_psycho_cum{s,j,k,repeat_sim}(:,2), H{k,j}, x, cum_gaussfit(psy_perf{s,j,k}, x),  F{j} );
        set(gca,'XTickMode','auto'); set(gca,'YTickMode','auto');
        xlabel('Heading Angle');   
        ylim([0,1]);
        ylabel('Percent Rightward Choices');
        hold on;
        legend_txt{j*2-1} = [num2str(unique_delta(j))];
        legend_txt{j*2} = [''];
    end

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % output some text of basic parameters in the figure
    axes('position',[0.2,0.8, 0.6,0.15] );
    xlim( [0,50] );
    ylim( [2,10] );
    text(0, 10, FILE);
    text(15,11, ['amplitude = ' num2str(unique_amplitude)]);
    text(15,10, ['stimtype = ' num2str(unique_stim_type(s))]);
    text(30,10, ['azimuth = ' num2str(mean(unique_azimuth))]); % mean, because actual 'AZIMUTH' varies with conflict angle
    text(45,10, ['repeats = ' num2str(num_reps)]);
    text(0,8.3, con_txt);
    text(8,8.3, con_txt_2);
    text(18,8.3, 'bias');
    text(28,8.3, 'thresh');
    text(37,8.3, '%correct');
    text(47,8.3, 'Wves-act');

    for j = 1:length(unique_delta)    % <-- currently conflict_angle
        text(0,8-j, num2str(unique_delta(j)));
        text(8,8-j,num2str(unique_coherence(k)));
        text(18,8-j,num2str(Bias_psy{s,j,k}) );
        text(28,8-j,num2str(Thresh_psy{s,j,k}) );
        text(37,8-j,num2str(num2str(mean(mean([1-fit_data_psycho_cum{s,j,k,repeat_sim}(1:3,2) fit_data_psycho_cum{s,j,k,repeat_sim}(5:7,2)])))));
        if j == 1
            text(47,8-j,num2str(Wves_actual_minus(k)) );
        elseif j == 3
            text(47,8-j,num2str(Wves_actual_plus(k)) );
        else
            text(48,8-j,num2str(Wves_actual(k)) );
        end
    end

    axis off;
%     print;
%     close;

end % end combined
% pause;
end % plotflag

% Tvislow(repeat_sim) = Thresh_psy{2,2,1};
% Tvishigh(repeat_sim) = Thresh_psy{2,2,2};
% Tves(repeat_sim) = Thresh_psy{1,2,1};
% Tcomblow(repeat_sim) = Thresh_psy{3,2,1};
% Tcombhigh(repeat_sim) = Thresh_psy{3,2,2};

W_ves_actual(:,repeat_sim) = Wves_actual;
W_ves_predicted(:,repeat_sim) = Wves_predicted;
% Thresh_actual(:,repeat_sim) = thresh_actual_all;
Thresh_actual(:,repeat_sim) = thresh_actual_d0;
Thresh_predicted(:,repeat_sim) = thresh_predicted;
if unique_stim_type == 1
    Thresh_ves(repeat_sim, :) = horzcat(Thresh_psy{1,1,:})
    Thresh_vis(repeat_sim, :) = zeros(repeat_sim, length(unique_coherence));
elseif unique_stim_type == 2
    Thresh_ves(repeat_sim, :) = zeros(repeat_sim, length(unique_coherence));
    Thresh_vis(repeat_sim, :) = horzcat(Thresh_psy{2,2,:})
else
    Thresh_ves(repeat_sim, :) = horzcat(Thresh_psy{1,2,:});
    Thresh_vis(repeat_sim, :) = horzcat(Thresh_psy{2,2,:});
end

if length(unique_stim_type) > 2
    e_cell_psy_raw{repeat_sim} = [unique_heading right_pct{2,2,1}' right_pct{2,2,2}' right_pct{1,2,1}' right_pct{3,1,1}' right_pct{3,2,1}' right_pct{3,3,1}' right_pct{3,1,2}' right_pct{3,2,2}' right_pct{3,3,2}'];
    xi = -10:0.1:10;
    e_cell_psy_fit{repeat_sim} = [xi' cum_gaussfit(psy_perf{2,2,1}, xi)' cum_gaussfit(psy_perf{2,2,2}, xi)' cum_gaussfit(psy_perf{1,2,1}, xi)' cum_gaussfit(psy_perf{3,1,1}, xi)' cum_gaussfit(psy_perf{3,2,1}, xi)' cum_gaussfit(psy_perf{3,3,1}, xi)' cum_gaussfit(psy_perf{3,1,2}, xi)' cum_gaussfit(psy_perf{3,2,2}, xi)' cum_gaussfit(psy_perf{3,3,2}, xi)'];
end
% 
% if save_data
% datestamp = date;
% % save(['likelihood_results_' batch_name '_CORR' num2str(noise_corr_max) '_' datestamp '.mat'])
% save(['likelihood_results_' batch_name '_CORR' num2str(noise_corr_max) '_' noisemodel '_' datestamp '_done.mat'])
% end

toc
end % end repeat_sim

Thresh_actual(Thresh_actual<0.1) = NaN
Thresh_predicted(Thresh_predicted<0.1) = NaN
Thresh_vis(Thresh_vis<0.1) = NaN
Thresh_ves(Thresh_ves<0.1) = NaN

% figure; plot(W_ves_predicted(1,:),W_ves_actual(1,:),'x',W_ves_predicted(2,:),W_ves_actual(2,:),'o',[0 1],[0 1],'k--');
% axis square; title('Weights');
% figure; plot(Thresh_predicted(1,:),Thresh_actual(1,:),'x',Thresh_predicted(2,:),Thresh_actual(2,:),'o',[0 2],[0 2],'k--');
% hold on; plot(Thresh_predicted(1,:),Thresh_vis(:,1),'rx',Thresh_predicted(2,:),Thresh_vis(:,2),'ro');
% plot(Thresh_predicted(1,:),Thresh_ves(:,1),'kx',Thresh_predicted(2,:),Thresh_ves(:,2),'ko');
% axis square; title('Thresholds');

if repeat_sim > 2
    figure; errorbar([1 2], nanmean(W_ves_predicted'), nanstd(W_ves_predicted')/sqrt(repeat_sim - sum(isnan(W_ves_predicted(2,:)))), 'bo-')
    hold on; errorbar([1 2], mean(W_ves_actual'), std(W_ves_actual')/sqrt(repeat_sim), 'ro-');
    xlim([0.5 2.5]); set(gca,'xtick', [1 2], 'xticklabel',[16 60]);
    title(['Weights, wich=' num2str(wichmann)]); legend('Predicted','Actual');
    
    figure; errorbar([1 2], nanmean(Thresh_predicted'), nanstd(Thresh_predicted')/sqrt(repeat_sim - sum(isnan(Thresh_predicted(2,:)))), 'co-')
    hold on; errorbar([1 2], nanmean(Thresh_actual'), nanstd(Thresh_actual')/sqrt(repeat_sim - sum(isnan(Thresh_actual(2,:)))), 'bo-');
    errorbar([1 2], nanmean(Thresh_vis), nanstd(Thresh_vis)/sqrt(repeat_sim - sum(isnan(Thresh_vis(:,2)))), 'rx-');
    errorbar([1 2], nanmean(Thresh_ves), nanstd(Thresh_ves)/sqrt(repeat_sim - sum(isnan(Thresh_ves(:,2)))), 'kx-');
    xlim([0.5 2.5]); set(gca,'xtick', [1 2], 'xticklabel',[16 60]);
    title(['Thresholds, wich=' num2str(wichmann)]); legend('Predicted','Actual','Visual','Vestib');
else
    figure; plot([1 2], W_ves_predicted, 'bo-')
    hold on; plot([1 2], W_ves_actual, 'ro-');
    xlim([0.5 2.5]); set(gca,'xtick', [1 2], 'xticklabel',[16 60]);
    title(['Weights, wich=' num2str(wichmann)]); legend('Predicted','Actual');
    
    figure; hold on;
    plot([1 2], Thresh_predicted, 'co-'); plot([1 2], Thresh_actual, 'bo-');
    plot([1 2], Thresh_vis', 'rx-'); plot([1 2], Thresh_ves', 'kx-');
    xlim([0.5 2.5]); set(gca,'xtick', [1 2], 'xticklabel',[16 60]);
    title(['Thresholds, wich=' num2str(wichmann)]); legend('Predicted','Actual','Visual','Vestib');
end

AAorigin1 = [nanmean(W_ves_predicted',1)' nanstd(W_ves_predicted',0,1)'/sqrt(repeat_sim - sum(isnan(W_ves_predicted(2,:)))) mean(W_ves_actual',1)' std(W_ves_actual',0,1)'/sqrt(repeat_sim)];
AAorigin2 = [nanmean(Thresh_predicted',1)' nanstd(Thresh_predicted',0,1)'/sqrt(repeat_sim - sum(isnan(Thresh_predicted(2,:)))) nanmean(Thresh_actual',1)' nanstd(Thresh_actual',0,1)'/sqrt(repeat_sim - sum(isnan(Thresh_actual(2,:)))) nanmean(Thresh_ves,1)' nanstd(Thresh_ves,0,1)'/sqrt(repeat_sim - sum(isnan(Thresh_ves(:,2)))) nanmean(Thresh_vis,1)' nanstd(Thresh_vis,0,1)'/sqrt(repeat_sim - sum(isnan(Thresh_vis(:,2))))];
AAorigin = [AAorigin1 AAorigin2];

e_cell_sum_fit = zeros(size(e_cell_psy_fit{repeat_sim}));
e_cell_sum_raw = zeros(size(e_cell_psy_raw{repeat_sim}));
for h = 1:repeat_sim
	e_cell_sum_fit = e_cell_sum_fit + e_cell_psy_fit{h};
    e_cell_sum_raw = e_cell_sum_raw + e_cell_psy_raw{h};
end
AAAe_cell_mean_fit = e_cell_sum_fit / repeat_sim;
AAAe_cell_mean_raw = e_cell_sum_raw / repeat_sim;

if save_data
datestamp = date;
% save(['likelihood_results_' batch_name '_CORR' num2str(noise_corr_max) '_' datestamp '.mat'])
save(['likelihood_results_' batch_name '_CORR' num2str(noise_corr_max) '_' noisemodel '_' datestamp '_sponsub' num2str(subtract_baseline) '.mat'])
end

clock