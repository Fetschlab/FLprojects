% population_likelihood.m 
% Trial-by-trial simulation of the behavioral experiment, computing
% likelihood functions by drawing spike counts from real neuronal data.

% THIS VERSION SIMULATES RESPONSES WITH A DESIRED NOISE CORRELATION LEVEL,
% ASSUMING GAUSSIAN NOISE WITH A VARIANCE = FANO*MEAN

%close all
clear all
% global cells_to_delete
dbstop if error

% BASE_PATH = 'C:\Program Files\MATLAB\R2007a\work\';
BASE_PATH = '/home/chris/backup/work/'; % if running on linux cluster

batch_name = 'cong_both_gt0.4'
% batch_name = 'cong_yos_gt0.4'
% batch_name = 'cong_web_gt0.4'
% batch_name = 'all_cells'
% batch_name = 'opp_both_0.4'

% finer timecourse steps (200 ms):
% timefolder = 'tc_400_600\';
% timefolder = 'tc_600_800\';
% timefolder = 'tc_800_1000\';
% timefolder = 'tc_1000_1200\';
% timefolder = 'tc_1200_1400\';
% timefolder = 'tc_1400_1600\';
% timefolder = 'tc_1600_1800\';
% timefolder = 'tc_1800_2000\';

% load([BASE_PATH 'condition_list_neuro.mat']);
load([BASE_PATH 'condition_list_interp.mat']);

if BASE_PATH(1)=='/' % mac/linux
    load([BASE_PATH 'neurons_indiv_trials/' batch_name '.mat']);
    files = textread([BASE_PATH 'neurons (temp)/' batch_name '.txt'],'%q');
else
    load([BASE_PATH 'neurons_indiv_trials\' batch_name '.mat']);
%     load([BASE_PATH 'neurons_indiv_trials\' timefolder batch_name '.mat']);
    files = textread([BASE_PATH 'neurons (temp)\' batch_name '.txt'],'%q');
end

num_stim_reps = 75;
num_resp_trials = 200;

save_data = 1;
run_resample = 0;
Jan = 1;

use_logL = 0;
use_fake_comb = 1;
    weights_to_apply = 'FishOpt'; include_DC = 0; mult_nonlin = 0; quad_nonlin = 0; % these params must be fixed at 0
        add_a_constant = 1;
    % weights_to_apply = 'Morgan';
%     weights_to_apply = 'NewMorgan'; include_DC = 1; mult_nonlin = 0; quad_nonlin = 0;
    % weights_to_apply = 'fixed';
make_fake_singles = 0;
keep_real_trials = 0; % i.e. coarse headings
subtract_baseline = 0;

noisemodel = 'Poisson'    % refers to generation of interpolated single-cue trials
% noisemodel = 'Gaussian'
% noisemodel = 'Gamma'

% likemodel = 'Poisson'
% likemodel = 'Gaussian'
% likemodel = 'Gamma'
likemodel = 'histogram'  % remember, no need to use histogram unless generating fake combs as weighted sums of single trials


cell_num = length(resp_mean_all);
% method = 'spline';
method = 'linear';

wichmann = 0;
plotflag = 1; % for psychometrics
FILE = ['SIM'];
repeat_sim = 1;

if keep_real_trials
    load('C:\Program Files\MATLAB\R2007a\work\condition_list_neuro.mat')
end

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
total_tr = length(stimtype) * num_stim_reps;
stim = zeros(1,total_tr); coh = zeros(1,total_tr);
hdg = zeros(1,total_tr); delt= zeros(1,total_tr);
choice = zeros(1,total_tr);

%-----------------------------------------------------------------------
% some preliminaries:
% X = unique_heading';
X = [-10 -3.5 -1.3 0 1.3 3.5 10];
% XI = unique_heading(1):0.2:unique_heading(end);
% XI = [-10 -3.5:0.2:-1.3 -1.2:0.1:1.2 1.3:0.2:3.5 10];
XI = unique_heading';

% replace resp with fake data if selected
if make_fake_singles
%    [resp_mean_all, makefake_baselines] = make_fake_data_for_poplike;
%    [resp_mean_all, makefake_baselines] = make_fake_data_for_poplike_realistic;
	[resp_mean_all, makefake_baselines, Wratio_opt] = make_fake_data_for_poplike_Alex;
    cell_num = length(resp_mean_all);
end

% % SUPER-KLUGE: temporarily erase some cells that are problematic for
% % fisher-optimal weights (ves & vis_low just barely opposite in slope: 15 16 21 22 36)
% super_kluge = [15 16 21 22 36];
super_kluge = [16 40 44]; % (NewMorgans: black bars only (sign. well-fit), excudes: 16, 40, 44)
if ~make_fake_singles && use_fake_comb && strcmp(batch_name,'cong_both_gt0.4') % && strcmp(weights_to_apply,'FishOpt')  % for now, remove these cells even for Morgan to make comparisons
    resp_mean_all(super_kluge) = []; resp_trial_all(super_kluge) = []; files(super_kluge) = [];
    cell_num = length(resp_mean_all)
end

% convert mean firing rates to Yong's format, and interpolate
for n = 1:cell_num
    for s = 1:3; % stim type
        for j = 1:3 % conflict angle
            for k = 1:2 % coherence
                if ~all(isnan(resp_mean_all{n}{s,j,k}))
                    firing{s,j,k}(n,:) = interp1(X,resp_mean_all{n}{s,j,k},XI,method);
                    param1_interp{s,j,k} = [];
                    param2_interp{s,j,k} = [];
                end
            end
        end
    end
end

xx = round(unique_heading'*100)/100; % rounding to 1.23, but okay
xi_step = 0.1;
xi = xx(1):xi_step:xx(end);
xi_short = -10:0.1:10;

% must replace real singles with linear fits in order to calculate optimal neural weights  
resp_mean_all_orig = resp_mean_all;
resp_trial_all_orig = resp_trial_all;
firing_orig = firing;

if use_fake_comb
if strcmp(weights_to_apply,'FishOpt')

%     clear firing resp_mean_all
    
    for n = 1:cell_num    
        
        % ugh, this below won't work (using real data instead of
        % fitting) because the weights vary in crazy ways with heading
        %             f1_low = interp1(X,resp_mean_all{n}{2,2,1},xi_short);
        %             f1_high = interp1(X,resp_mean_all{n}{2,2,2},xi_short);
        %             f2 = interp1(X,resp_mean_all{n}{1,2,1},xi_short);

        % instead just fit a line to each condition separately
        b1_low = regress(resp_mean_all_orig{n}{2,2,1}',[X' ones(length(X),1)]);
        f1_low = b1_low(1)*xi_short + b1_low(2);            
        b1_high = regress(resp_mean_all_orig{n}{2,2,2}',[X' ones(length(X),1)]);
        f1_high = b1_high(1)*xi_short + b1_high(2);
        b2 = regress(resp_mean_all_orig{n}{1,2,1}',[X' ones(length(X),1)]);
        f2 = b2(1)*(xi_short) + b2(2);

        % technically, weights are only correct if you replace singles with
        % these linear fits, but comment this out to try keeping orig reals
% %         firing{2,2,1}(n,:) = b1_low(1)*xx + b1_low(2);
% %         resp_mean_all{n}{2,2,1} = b1_low(1)*X + b1_low(2);            
% %         firing{2,2,2}(n,:) = b1_high(1)*xx + b1_high(2);
% %         resp_mean_all{n}{2,2,2} = b1_high(1)*X + b1_high(2);
% %         firing{1,2,1}(n,:) = b2(1)*xx + b2(2);
% %         resp_mean_all{n}{1,2,1} = b2(1)*X + b2(2);

        % also need comb now
        b3_low = regress(resp_mean_all_orig{n}{3,2,1}',[X' ones(length(X),1)]);
        f3_low = b3_low(1)*xi_short + b3_low(2);            
        b3_high = regress(resp_mean_all_orig{n}{3,2,2}',[X' ones(length(X),1)]);
        f3_high = b3_high(1)*xi_short + b3_high(2);
        firing{3,2,1}(n,:) = b3_low(1)*xx + b3_low(2);
        resp_mean_all{n}{3,2,1} = b3_low(1)*X + b3_low(2);            
        firing{3,2,2}(n,:) = b3_high(1)*xx + b3_high(2);
        resp_mean_all{n}{3,2,2} = b3_high(1)*X + b3_high(2);

% %         % relax assumption of linear singles (try spline at zero-heading)
% %         f1_low = interp1(X,resp_mean_all_orig{n}{2,2,1},xi_short,'linear');
% %         f1_high = interp1(X,resp_mean_all_orig{n}{2,2,2},xi_short,'linear');
% %         f2 = interp1(X,resp_mean_all_orig{n}{1,2,1},xi_short,'linear');
        
        f1prime_low = diff(f1_low);
        f1prime_high = diff(f1_high);
        f2prime = diff(f2);
        
        % kluge a couple cells with opposite slope for low/high coh, 
        % which are basically flat at low coh ("super-kluge")
        if sign(f1prime_low(100)) ~= sign(f1prime_high(100))
            f1prime_low(:) = .005 * sign(f1prime_high(100));
        end
        
        % Alex's Eq. 1.2
        Wratio_opt_low(n,:) = (f1_low(2:end).*f2prime) ./ (f1prime_low.*f2(2:end));
        Wratio_opt_high(n,:) = (f1_high(2:end).*f2prime) ./ (f1prime_high.*f2(2:end));
                
        % still may need to try Wratio that can vary with hdg, but for
        % now take the mean of some range (or just eval at zero)
%         a = rand+0.5;
        a = 1;
        x_low{n} = a*[mean(Wratio_opt_low(n,99:101)) 1]; % 65:135 is -3.5 to 3.5 hdg, 100 is zero
%        a = 1;  % CAREFUL! not sure it's legal for a to vary with coh
        x_high{n} = a*[mean(Wratio_opt_high(n,99:101)) 1];

        Wratio_opt{n} = [mean(Wratio_opt_low(n,:)) mean(Wratio_opt_high(n,:))];
            
%             % here's where we try to vary weights with heading
%             a = 0.85;
%             x_low{n} = ones(length(xi_short),2); x_high{n} = ones(length(xi_short),2);
%             x_low{n}(2:end,:) = a*[Wratio_opt_low(n,:)' ones(length(xi_short)-1,1)];
%             x_high{n}(2:end,:) = a*[Wratio_opt_high(n,:)' ones(length(xi_short)-1,1)];
%             x_low{n}(1,:) = x_low{n}(2,:); x_high{n}(1,:) = x_high{n}(2,:);

    end
    
    xl = cell2mat(x_low); xl=xl';
    xh = cell2mat(x_high); xh=xh';
%     xl(isinf(xl))=NaN;
%     xh(isinf(xh))=NaN;
    data1 = [xl(1:2:end-1) xh(1:2:end-1)];
%     asdf
%     figure; hist(xl(1:2:end-1));
%     figure; hist(xh(1:2:end-1)); pause;
     
elseif strcmp(weights_to_apply,'Morgan');

% % % %     clear firing resp_mean_all

    for n = 1:cell_num

%         % I guess also try this for Morgan case..  something is fishy with
%         % assuming linear tuning, so this is a place to start...
%         b1_low = regress(resp_mean_all_orig{n}{2,2,1}',[X' ones(length(X),1)]);
%         f1_low = b1_low(1)*xi_short + b1_low(2);            
%         b1_high = regress(resp_mean_all_orig{n}{2,2,2}',[X' ones(length(X),1)]);
%         f1_high = b1_high(1)*xi_short + b1_high(2);
%         b2 = regress(resp_mean_all_orig{n}{1,2,1}',[X' ones(length(X),1)]);
%         f2 = b2(1)*(xi_short) + b2(2);
% 
%         firing{2,2,1}(n,:) = b1_low(1)*xx + b1_low(2);
%         resp_mean_all{n}{2,2,1} = b1_low(1)*X + b1_low(2);            
%         firing{2,2,2}(n,:) = b1_high(1)*xx + b1_high(2);
%         resp_mean_all{n}{2,2,2} = b1_high(1)*X + b1_high(2);
%         firing{1,2,1}(n,:) = b2(1)*xx + b2(2);
%         resp_mean_all{n}{1,2,1} = b2(1)*X + b2(2);
% 
%         % also need comb now
%         b3_low = regress(resp_mean_all_orig{n}{3,2,1}',[X' ones(length(X),1)]);
%         f3_low = b3_low(1)*xi_short + b3_low(2);            
%         b3_high = regress(resp_mean_all_orig{n}{3,2,2}',[X' ones(length(X),1)]);
%         f3_high = b3_high(1)*xi_short + b3_high(2);
%         firing{3,2,1}(n,:) = b3_low(1)*xx + b3_low(2);
%         resp_mean_all{n}{3,2,1} = b3_low(1)*X + b3_low(2);            
%         firing{3,2,2}(n,:) = b3_high(1)*xx + b3_high(2);
%         resp_mean_all{n}{3,2,2} = b3_high(1)*X + b3_high(2);
        
        % LOW COHERENCE
        xdata = [resp_mean_all{n}{1,2,1}-subtract_baseline*baseline_rate ; resp_mean_all{n}{2,2,1}-subtract_baseline*baseline_rate];
        ydata = resp_mean_all{n}{3,2,1}-subtract_baseline*baseline_rate;
        XX = xdata';
        YY = ydata';
        [b,bint,r,rint,stats] = regress(YY,XX);
        x_low{n} = [b(1) b(2)];
        [r,p] = corrcoef(x_low{n}(1)*xdata(1,:)+x_low{n}(2)*xdata(2,:), ydata);
        R_corr_low(n) = r(1,2);
        P_corr_low(n) = p(1,2);

        % HIGH COHERENCE
        xdata = [resp_mean_all{n}{1,2,1}-subtract_baseline*baseline_rate ; resp_mean_all{n}{2,2,2}-subtract_baseline*baseline_rate];
        ydata = resp_mean_all{n}{3,2,2}-subtract_baseline*baseline_rate;
        XX = xdata';
        YY = ydata';
        [b,bint,r,rint,stats] = regress(YY,XX);
        x_high{n} = [b(1) b(2)];
        [r,p] = corrcoef(x_high{n}(1)*xdata(1,:)+x_high{n}(2)*xdata(2,:), ydata);
        R_corr_high(n) = r(1,2);
        P_corr_high(n) = p(1,2);

%         % visually inspect fit quality
%         figure(12); clf;
%         fake_low = x_low{n}(1)*resp_mean_all{n}{1,2,1} + x_low{n}(2)*resp_mean_all{n}{2,2,1};
%         fake_high = x_high{n}(1)*resp_mean_all{n}{1,2,1} + x_high{n}(2)*resp_mean_all{n}{2,2,2};
%         plot(X,resp_mean_all{n}{3,2,1},'c-o',X,fake_low,'c--s'); hold on;
%         plot(X,resp_mean_all{n}{3,2,2},'b-o',X,fake_high,'b--s');
%         pause;
        
    end

end

% lows = [R_corr_low' P_corr_low']
% highs = [R_corr_high' P_corr_high']
% pause

% % go ahead and make combined means here (e.g., for fitting param1+2), even though
% % they will typically be replaced by weighted sums of spike counts later
% for n = 1:cell_num
%     for j = 1:3 % conflict angle
%         for k = 1:2 % coherence
%             if ~isempty(firing{s,j,k})
%                 for i=1:length(XI)
%                     firing_trial{s,j,k,i}(n,1:num_resp_trials) = poissrnd(firing{s,j,k}(n,i),1,num_resp_trials);
%                 end
%                 param1_interp{s,j,k} = [];
%                 param2_interp{s,j,k} = [];
%             end
%         end
%     end
% end

end

% now generate trials for interpolated headings, using specified noise model
% *NOTE: combined trials made here will be replaced with weighted sums of single-
% cue trials, if 'use_fake_comb' option is specified
if strcmp(noisemodel,'Poisson')

    for n = 1:cell_num
        n
        for s = 1:3; % stim type
            for j = 1:3 % conflict angle
                for k = 1:2 % coherence
                    if ~isempty(firing{s,j,k})
                        for i=1:length(XI)
                            firing_trial{s,j,k,i}(n,1:num_resp_trials) = poissrnd(firing{s,j,k}(n,i),1,num_resp_trials);
                            % if replacing singles with linear fits, can get negative
                            % mean values and hence NaNs that should be zeroes
                            firing_trial{s,j,k,i}(n,isnan(firing_trial{s,j,k,i}(n,:))) = 0;
                        end
                        param1_interp{s,j,k} = [];
                        param2_interp{s,j,k} = [];
                    end
                end
            end
        end
    end
    
elseif strcmp(noisemodel,'Gaussian') % **fit the real data to get params for Gaussian or Gamma

    for n = 1:cell_num
        n
        for s = 1:3; % stim type
            for j = 1:3 % conflict angle
                for k = 1:2 % coherence
                    if ~isempty(firing{s,j,k})
                        for h = 1:length(X)
                            [muhat(h),sigmahat(h)] = normfit(resp_trial_all_orig{n}{s,j,k,h});
                        end
                        muhat_interp(n,:) = interp1(X,muhat,XI,'linear');
                        sigmahat_interp(n,:) = interp1(X,sigmahat,XI,'linear');
                        sigmahat_interp(n,sigmahat_interp(n,:)==0) = 0.1; % sigma can't be zero or normpdf will give nans
                        if keep_real_trials
                            % copy correct params to a variable to pass into compute_likelihood
                            param1_interp{s,j,k}(n,:) = muhat_interp(n,:);
                            param2_interp{s,j,k}(n,:) = sigmahat_interp(n,:);
                        else
                            for i=1:length(XI)
                                firing_trial{s,j,k,i}(n,1:num_resp_trials) = normrnd(muhat_interp(n,i),sigmahat_interp(n,i),1,num_resp_trials);
                                firing_trial{s,j,k,i}(n,firing_trial{s,j,k,i}(n,:)<0) = 0;
                                % because we had to clip responses at zero, must now re-fit 
                                % before coping params to a variable to pass into compute_likelihood
                                [muhat2(i),sigmahat2(i)] = normfit(firing_trial{s,j,k,i}(n,:));        
                            end
                            param1_interp{s,j,k}(n,:) = muhat2;
                            param2_interp{s,j,k}(n,:) = sigmahat2;
                        end
                    end
                end
            end
        end
    end
    
elseif strcmp(noisemodel,'Gamma')

    for n = 1:cell_num
        n
        for s = 1:3; % stim type
            for j = 1:3 % conflict angle
                for k = 1:2 % coherence
                    if ~isempty(firing{s,j,k})
                        clear phat
                        for h = 1:length(X)
                            % kluge: if all resps are equal, param1 becomes inf; so, perturb one trial by 5% in that case
                            if std(resp_trial_all{n}{s,j,k,h})<1e-3
                                resp_trial_all{n}{s,j,k,h}(1) = resp_trial_all{n}{s,j,k,h}(1) + 0.05*resp_trial_all{n}{s,j,k,h}(1);
                            end
                            phat(h,:) = gamfit(resp_trial_all{n}{s,j,k,h});
                            phat(h,isnan(phat(h,:))) = 0; % protect against NaNs (happens when all resp are zero)                                                        
                        end
                        shape_interp(n,:) = interp1(X,phat(:,1),XI,'linear');
                        shape_interp(n,shape_interp(n,:)<=0) = 0.1; % protect against zeros and negatives (they give nans)
                        scale_interp(n,:) = interp1(X,phat(:,2),XI,'linear');
                        scale_interp(n,scale_interp(n,:)<=0) = 0.1;
                        if keep_real_trials
                            % copy correct params to a variable to pass into compute_likelihood
                            param1_interp{s,j,k}(n,:) = shape_interp(n,:);
                            param2_interp{s,j,k}(n,:) = scale_interp(n,:);
                        else                            
                            for i=1:length(XI)
                                firing_trial{s,j,k,i}(n,1:num_resp_trials) = gamrnd(shape_interp(n,i),scale_interp(n,i),1,num_resp_trials);
    %                             debatable whether I need to (or should) re-fit
                                phat(i,:) = gamfit(firing_trial{s,j,k,i}(n,:));
                            end
                            % copy correct params to a variable to pass into compute_likelihood
                            temp = phat';
                            param1_interp{s,j,k}(n,:) = temp(1,:);
                            param2_interp{s,j,k}(n,:) = temp(2,:);
                        end
                    end
                end
            end
        end
    end
    
end

% ------------------------------------------------------------------
if use_fake_comb
% Here we have the option to replace the real combined responses with 
% simulated ones (weighted sums of single-cue responses) to test hypotheses
% about what causes vestibular over-weighting in the model

% first need to interpolate single-cue tuning curves
% (even though we don't need the high resolution for this version of
% poplike -- we'll pick off the desired values later)

cell_num = length(resp_mean_all);
disp('generating fake data with specified weights');

for v = 1 %:length(Wves_range)
    v
    for w = 1 %:length(Wvis_range)
        w
%         x_low{n} = [0.63 0.57];
%         x_high{n} = [0.32 0.86];

%         x_low{n} = [1 1];
%         x_high{n} = [1 1];

%         wtemp = cell2mat(Wratio_opt);
%         wtemp2(:,1) = wtemp(1:2:length(wtemp));
%         wtemp2(:,2) = wtemp(2:2:length(wtemp));
%         wtemp3 = mean(wtemp2);
%         a = 1;
%         x_low{n} = a*[wtemp3(1) 1];
%         x_high{n} = a*[wtemp3(2) 1];

        for n = 1:cell_num
            
            subtract_baseline = 0; % kluge for add_a_constant

            if strcmp(weights_to_apply,'fixed');
                % white crosses from Fig. 8a,c
                x_low{n} = [0.7 0.37];
                x_high{n} = [0.29 0.89];
            end

            if make_fake_singles
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
            
            if ~strcmp(noisemodel,'Poisson')
                param1_vislow{n} = interp1(xx,param1_interp{2,2,1}(n,:),xi,method);
                param1_vishigh{n} = interp1(xx,param1_interp{2,2,2}(n,:),xi,method);
                param1_ves{n} = interp1(xx,param1_interp{1,2,1}(n,:),xi,method);
                param2_vislow{n} = interp1(xx,param2_interp{2,2,1}(n,:),xi,method);
                param2_vishigh{n} = interp1(xx,param2_interp{2,2,2}(n,:),xi,method);
                param2_ves{n} = interp1(xx,param2_interp{1,2,1}(n,:),xi,method);
            end
            
            % figure;plot(xx,resp_mean_all{n}{2,2,1},'bo',xx,resp_mean_all{n}{2,2,2},'co',xx,resp_mean_all{n}{1,2,1},'ro');
            % hold on; plot(xi,yi_vislow{n},'b-',xi,yi_vishigh{n},'c-',xi,yi_ves{n},'r-');

            % also need to extrapolate to +/- delta/2 degrees outside largest heading
            x_chunk = xi(1 : find(xi==xx(2))-1);
            y_chunk = yi_vislow{n}(1 : find(xi==xx(2))-1);
            P = polyfit(x_chunk,y_chunk,1);
            xi_seg1 = xx(1)-4/2 : xi_step : xx(1)-xi_step;
            yi_seg1 = P(1)*xi_seg1 + P(2);
            x_chunk = xi(find(xi==xx(end-1)) : end);
            y_chunk = yi_vislow{n}(find(xi==xx(end-1)) : end);
            P = polyfit(x_chunk,y_chunk,1);
            xi_seg2 = xx(end)+xi_step : xi_step : xx(end)+4/2;
            yi_seg2 = P(1)*xi_seg2 + P(2);
            yi_vislow{n} = [yi_seg1 yi_vislow{n} yi_seg2];

            x_chunk = xi(1 : find(xi==xx(2))-1);
            y_chunk = yi_vishigh{n}(1 : find(xi==xx(2))-1);
            P = polyfit(x_chunk,y_chunk,1);
            xi_seg1 = xx(1)-4/2 : xi_step : xx(1)-xi_step;
            yi_seg1 = P(1)*xi_seg1 + P(2);
            x_chunk = xi(find(xi==xx(end-1)) : end);
            y_chunk = yi_vishigh{n}(find(xi==xx(end-1)) : end);
            P = polyfit(x_chunk,y_chunk,1);
            xi_seg2 = xx(end)+xi_step : xi_step : xx(end)+4/2;
            yi_seg2 = P(1)*xi_seg2 + P(2);
            yi_vishigh{n} = [yi_seg1 yi_vishigh{n} yi_seg2];

            x_chunk = xi(1 : find(xi==xx(2))-1);
            y_chunk = yi_ves{n}(1 : find(xi==xx(2))-1);
            P = polyfit(x_chunk,y_chunk,1);
            xi_seg1 = xx(1)-4/2 : xi_step : xx(1)-xi_step;
            yi_seg1 = P(1)*xi_seg1 + P(2);
            x_chunk = xi(find(xi==xx(end-1)) : end);
            y_chunk = yi_ves{n}(find(xi==xx(end-1)) : end);
            P = polyfit(x_chunk,y_chunk,1);
            xi_seg2 = xx(end)+xi_step : xi_step : xx(end)+4/2;
            yi_seg2 = P(1)*xi_seg2 + P(2);
            yi_ves{n} = [yi_seg1 yi_ves{n} yi_seg2];
            
            % ugly as hell, but need to repeat all three above blocks of
            % code 2 more times for param1 and param2 interpolation
                if ~strcmp(noisemodel,'Poisson')            
                x_chunk = xi(1 : find(xi==xx(2))-1);
                y_chunk = param1_vislow{n}(1 : find(xi==xx(2))-1);
                P = polyfit(x_chunk,y_chunk,1);
                xi_seg1 = xx(1)-4/2 : xi_step : xx(1)-xi_step;
                param1_seg1 = P(1)*xi_seg1 + P(2);
                x_chunk = xi(find(xi==xx(end-1)) : end);
                y_chunk = param1_vislow{n}(find(xi==xx(end-1)) : end);
                P = polyfit(x_chunk,y_chunk,1);
                xi_seg2 = xx(end)+xi_step : xi_step : xx(end)+4/2;
                param1_seg2 = P(1)*xi_seg2 + P(2);
                param1_vislow{n} = [param1_seg1 param1_vislow{n} param1_seg2];

                x_chunk = xi(1 : find(xi==xx(2))-1);
                y_chunk = param1_vishigh{n}(1 : find(xi==xx(2))-1);
                P = polyfit(x_chunk,y_chunk,1);
                xi_seg1 = xx(1)-4/2 : xi_step : xx(1)-xi_step;
                param1_seg1 = P(1)*xi_seg1 + P(2);
                x_chunk = xi(find(xi==xx(end-1)) : end);
                y_chunk = param1_vishigh{n}(find(xi==xx(end-1)) : end);
                P = polyfit(x_chunk,y_chunk,1);
                xi_seg2 = xx(end)+xi_step : xi_step : xx(end)+4/2;
                param1_seg2 = P(1)*xi_seg2 + P(2);
                param1_vishigh{n} = [param1_seg1 param1_vishigh{n} param1_seg2];

                x_chunk = xi(1 : find(xi==xx(2))-1);
                y_chunk = param1_ves{n}(1 : find(xi==xx(2))-1);
                P = polyfit(x_chunk,y_chunk,1);
                xi_seg1 = xx(1)-4/2 : xi_step : xx(1)-xi_step;
                param1_seg1 = P(1)*xi_seg1 + P(2);
                x_chunk = xi(find(xi==xx(end-1)) : end);
                y_chunk = param1_ves{n}(find(xi==xx(end-1)) : end);
                P = polyfit(x_chunk,y_chunk,1);
                xi_seg2 = xx(end)+xi_step : xi_step : xx(end)+4/2;
                param1_seg2 = P(1)*xi_seg2 + P(2);
                param1_ves{n} = [param1_seg1 param1_ves{n} param1_seg2];            

                %
                %

                x_chunk = xi(1 : find(xi==xx(2))-1);
                y_chunk = param2_vislow{n}(1 : find(xi==xx(2))-1);
                P = polyfit(x_chunk,y_chunk,1);
                xi_seg1 = xx(1)-4/2 : xi_step : xx(1)-xi_step;
                param2_seg1 = P(1)*xi_seg1 + P(2);
                x_chunk = xi(find(xi==xx(end-1)) : end);
                y_chunk = param2_vislow{n}(find(xi==xx(end-1)) : end);
                P = polyfit(x_chunk,y_chunk,1);
                xi_seg2 = xx(end)+xi_step : xi_step : xx(end)+4/2;
                param2_seg2 = P(1)*xi_seg2 + P(2);
                param2_vislow{n} = [param2_seg1 param2_vislow{n} param2_seg2];

                x_chunk = xi(1 : find(xi==xx(2))-1);
                y_chunk = param2_vishigh{n}(1 : find(xi==xx(2))-1);
                P = polyfit(x_chunk,y_chunk,1);
                xi_seg1 = xx(1)-4/2 : xi_step : xx(1)-xi_step;
                param2_seg1 = P(1)*xi_seg1 + P(2);
                x_chunk = xi(find(xi==xx(end-1)) : end);
                y_chunk = param2_vishigh{n}(find(xi==xx(end-1)) : end);
                P = polyfit(x_chunk,y_chunk,1);
                xi_seg2 = xx(end)+xi_step : xi_step : xx(end)+4/2;
                param2_seg2 = P(1)*xi_seg2 + P(2);
                param2_vishigh{n} = [param2_seg1 param2_vishigh{n} param2_seg2];

                x_chunk = xi(1 : find(xi==xx(2))-1);
                y_chunk = param2_ves{n}(1 : find(xi==xx(2))-1);
                P = polyfit(x_chunk,y_chunk,1);
                xi_seg1 = xx(1)-4/2 : xi_step : xx(1)-xi_step;
                param2_seg1 = P(1)*xi_seg1 + P(2);
                x_chunk = xi(find(xi==xx(end-1)) : end);
                y_chunk = param2_ves{n}(find(xi==xx(end-1)) : end);
                P = polyfit(x_chunk,y_chunk,1);
                xi_seg2 = xx(end)+xi_step : xi_step : xx(end)+4/2;
                param2_seg2 = P(1)*xi_seg2 + P(2);
                param2_ves{n} = [param2_seg1 param2_ves{n} param2_seg2];
            end
            
            xi_new = xx(1)-4/2 : xi_step : xx(end)+4/2;
            % figure;plot(xx,resp_mean_all{n}{2,2,1},'bo',xx,resp_mean_all{n}{2,2,2},'co',xx,resp_mean_all{n}{1,2,1},'ro');
            % hold on; plot(xi_new,yi_vislow{n},'b-',xi_new,yi_vishigh{n},'c-',xi_new,yi_ves{n},'r-');

            % now simulate conflict conditions with the desired neural weights
            delta_shift = 4/2/xi_step;
            xi_new = round(xi_new*100)/100; % stupid rounding errors
            for h = 1:length(xx)
                vesminus(h) = yi_ves{n}(find(xi_new==xx(h))+delta_shift);
                veszero(h) = yi_ves{n}(find(xi_new==xx(h)));
                vesplus(h) = yi_ves{n}(find(xi_new==xx(h))-delta_shift);
                vislowminus(h) = yi_vislow{n}(find(xi_new==xx(h))-delta_shift);
                vislowzero(h) = yi_vislow{n}(find(xi_new==xx(h)));
                vislowplus(h) = yi_vislow{n}(find(xi_new==xx(h))+delta_shift);
                vishighminus(h) = yi_vishigh{n}(find(xi_new==xx(h))-delta_shift);
                vishighzero(h) = yi_vishigh{n}(find(xi_new==xx(h)));
                vishighplus(h) = yi_vishigh{n}(find(xi_new==xx(h))+delta_shift);
                
                if ~strcmp(noisemodel,'Poisson')
                    p1_vesminus(h) = param1_ves{n}(find(xi_new==xx(h))+delta_shift);
                    p1_veszero(h) = param1_ves{n}(find(xi_new==xx(h)));
                    p1_vesplus(h) = param1_ves{n}(find(xi_new==xx(h))-delta_shift);
                    p1_vislowminus(h) = param1_vislow{n}(find(xi_new==xx(h))-delta_shift);
                    p1_vislowzero(h) = param1_vislow{n}(find(xi_new==xx(h)));
                    p1_vislowplus(h) = param1_vislow{n}(find(xi_new==xx(h))+delta_shift);
                    p1_vishighminus(h) = param1_vishigh{n}(find(xi_new==xx(h))-delta_shift);
                    p1_vishighzero(h) = param1_vishigh{n}(find(xi_new==xx(h)));
                    p1_vishighplus(h) = param1_vishigh{n}(find(xi_new==xx(h))+delta_shift);

                    p2_vesminus(h) = param2_ves{n}(find(xi_new==xx(h))+delta_shift);
                    p2_veszero(h) = param2_ves{n}(find(xi_new==xx(h)));
                    p2_vesplus(h) = param2_ves{n}(find(xi_new==xx(h))-delta_shift);
                    p2_vislowminus(h) = param2_vislow{n}(find(xi_new==xx(h))-delta_shift);
                    p2_vislowzero(h) = param2_vislow{n}(find(xi_new==xx(h)));
                    p2_vislowplus(h) = param2_vislow{n}(find(xi_new==xx(h))+delta_shift);
                    p2_vishighminus(h) = param2_vishigh{n}(find(xi_new==xx(h))-delta_shift);
                    p2_vishighzero(h) = param2_vishigh{n}(find(xi_new==xx(h)));
                    p2_vishighplus(h) = param2_vishigh{n}(find(xi_new==xx(h))+delta_shift);                
                end
            end

            % extrapolation can cause negative values, leading to NaNs from poissrnd below
            vesminus(vesminus<0) = 0;
            veszero(veszero<0) = 0;
            vesplus(vesplus<0) = 0;
            vislowminus(vislowminus<0) = 0;
            vislowzero(vislowzero<0) = 0;
            vislowplus(vislowplus<0) = 0;
            vishighminus(vishighminus<0) = 0;
            vishighzero(vishighzero<0) = 0;
            vishighplus(vishighplus<0) = 0;

            
            if strcmp(weights_to_apply,'NewMorgan')
                % here is where we will attempt to fit new Morgan weights by
                % including the nonzero conflict conditions

                hdg_orig = [1 2 14 26 38 50 51];  % major kluge here -- must change if hdg axis (interpolation) changes!
                hdg_orig = [hdg_orig hdg_orig+51 hdg_orig+102];
                
                % Standard Model (4 params, 2 for each coh)
                % LOW COHERENCE
                ydata = [firing{3,1,1}(n,:) firing{3,2,1}(n,:) firing{3,3,1}(n,:)]-subtract_baseline*baseline_rate;
                if include_DC %(without, followed by with, multiplicative nonlinearity, followed by quadratic)
                    if mult_nonlin
                        xdata = [vesminus veszero vesplus ; vislowminus vislowzero vislowplus ; vesminus.*vislowminus veszero.*vislowzero vesplus.*vislowplus ; ones(1,length(ydata))];
                    elseif quad_nonlin
                        xdata = [vesminus veszero vesplus ; vislowminus vislowzero vislowplus ; vesminus.^2 veszero.^2 vesplus.^2 ; vislowminus.^2 vislowzero.^2 vislowplus.^2 ; ones(1,length(ydata))];
                    else
                        xdata = [vesminus veszero vesplus ; vislowminus vislowzero vislowplus ; ones(1,length(ydata))];
                    end
                else
                    if mult_nonlin
                        xdata = [vesminus veszero vesplus ; vislowminus vislowzero vislowplus ; vesminus.*vislowminus veszero.*vislowzero vesplus.*vislowplus];
                    elseif quad_nonlin
                        xdata = [vesminus veszero vesplus ; vislowminus vislowzero vislowplus ; vesminus.^2 veszero.^2 vesplus.^2 ; vislowminus.^2 vislowzero.^2 vislowplus.^2];
                    else
                        xdata = [vesminus veszero vesplus ; vislowminus vislowzero vislowplus];
                    end                    
                end
                xdata = xdata(:,hdg_orig); % this could be done after fitting, or before -- shouldn't change result
                ydata = ydata(hdg_orig);
                XX = xdata';
                YY = ydata';
                [b,bint,r,rint,stats] = regress(YY,XX);
                x_low{n} = b;
%                 xdata = xdata(:,hdg_orig); % this could be done after fitting, or before -- shouldn't change result
%                 ydata = ydata(hdg_orig);
                yfit = b'*xdata;
                
                [r,p] = corrcoef(yfit,ydata);
                R_corr_low(n) = r(1,2);
                P_corr_low(n) = p(1,2);
                % VAF = 1-SSE/SST (Morgan et al.)
                SSElow(n) = sum((yfit-ydata).^2);
                SST = sum((ydata - mean(ydata)).^2);
                VAF_low(n) = 1-SSElow(n)/SST;
                % adjusted R-squared
                MST = SST/(length(ydata)-1);
                MSE = SSElow(n)/(length(ydata)-2-1);
                R2adj_low(n) = 1-MSE/MST;
%                 R2adj_low(n) = 1-(SSE*(length(ydata)-1))/(SST*(length(ydata)-2));
                % save residual error at zero-conflict (to look for source of threshold difference)
                residualErrLow(n,1:7) = ydata(8:14) - yfit(8:14);
                tempcorr = corrcoef(ydata(8:14),X);
                if sign(tempcorr(1,2))==-1
                    residualErrLow(n,:) = fliplr(residualErrLow(n,:));
                end
                
%                 figure(1); clf; subplot(2,1,1); 
%                 plot(xx,vesminus,'r',xx,veszero,'b',xx,vesplus,'g');
%                 subplot(2,1,2); plot(xx,vislowminus,'r',xx,vislowzero,'b',xx,vislowplus,'g');
%                 figure(2); clf;
%                 plot(xx,ydata(1:51),'r',xx,ydata(52:102),'b',xx,ydata(103:153),'g');
%                 hold on; plot(xx,yfit(1:51),'r--',xx,yfit(52:102),'b--',xx,yfit(103:153),'g--');
%                 title([num2str(n) ' ' num2str(R_corr_low(n)^2)]);
%                 pause;
                
                % HIGH COHERENCE
                ydata = [firing{3,1,2}(n,:) firing{3,2,2}(n,:) firing{3,3,2}(n,:)]-subtract_baseline*baseline_rate;
                if include_DC
                    if mult_nonlin
                        xdata = [vesminus veszero vesplus ; vishighminus vishighzero vishighplus ; vesminus.*vishighminus veszero.*vishighzero vesplus.*vishighplus ; ones(1,length(ydata))];
                    elseif quad_nonlin
                        xdata = [vesminus veszero vesplus ; vishighminus vishighzero vishighplus ; vesminus.^2 veszero.^2 vesplus.^2 ; vishighminus.^2 vishighzero.^2 vishighplus.^2 ; ones(1,length(ydata))];
                    else
                        xdata = [vesminus veszero vesplus ; vishighminus vishighzero vishighplus ; ones(1,length(ydata))];
                    end
                else
                    if mult_nonlin
                        xdata = [vesminus veszero vesplus ; vishighminus vishighzero vishighplus ; vesminus.*vishighminus veszero.*vishighzero vesplus.*vishighplus];
                    elseif quad_nonlin
                        xdata = [vesminus veszero vesplus ; vishighminus vishighzero vishighplus ; vesminus.^2 veszero.^2 vesplus.^2 ; vishighminus.^2 vishighzero.^2 vishighplus.^2];
                    else
                        xdata = [vesminus veszero vesplus ; vishighminus vishighzero vishighplus];
                    end                    
                end
                xdata = xdata(:,hdg_orig); % this could be done after fitting, or before -- shouldn't change result
                ydata = ydata(hdg_orig);
                XX = xdata';
                YY = ydata';
                [b,bint,r,rint,stats] = regress(YY,XX);
                x_high{n} = b;
%                 xdata = xdata(:,hdg_orig); % this could be done after fitting, or before -- shouldn't change result
%                 ydata = ydata(hdg_orig);
                yfit = b'*xdata;

                [r,p] = corrcoef(yfit,ydata);
                R_corr_high(n) = r(1,2);
                P_corr_high(n) = p(1,2);
                % VAF = 1-SSE/SST (Morgan et al.)
                SSEhigh(n) = sum((yfit - ydata).^2);
                SST = sum((ydata - mean(ydata)).^2);
                VAF_high(n) = 1-SSEhigh(n)/SST;
                % adjusted R-squared
                MST = SST/(length(ydata)-1);
                MSE = SSEhigh(n)/(length(ydata)-2-1);
                R2adj_high(n) = 1-MSE/MST;
                % save residual error at zero-conflict (to look for source of threshold difference)
                residualErrHigh(n,1:7) = ydata(8:14) - yfit(8:14);
                tempcorr = corrcoef(ydata(8:14),X);
                if sign(tempcorr(1,2))==-1
                    residualErrHigh(n,:) = fliplr(residualErrHigh(n,:));
                end
                
%                 subplot(2,1,2); plot(ydata); hold on; plot(yfit,'g')
%                 title([num2str(n) ' ' num2str(R_corr_high(n)^2)]);                
%                 pause;

                RSS_IndependentWeights = SSElow(n) + SSEhigh(n);

                %*********************************************************
                % Yoked Weights model (2 params, 3 w/ dc)
                hdg_orig = [hdg_orig hdg_orig+153];
                ydata = [firing{3,1,1}(n,:) firing{3,2,1}(n,:) firing{3,3,1}(n,:) firing{3,1,2}(n,:) firing{3,2,2}(n,:) firing{3,3,2}(n,:)]-subtract_baseline*baseline_rate;
                if include_DC
                    xdata = [vesminus veszero vesplus vesminus veszero vesplus ; vislowminus vislowzero vislowplus vishighminus vishighzero vishighplus; ones(1,length(ydata))];
                    xdata = xdata(:,hdg_orig); % this could be done after fitting, or before -- shouldn't change result
                    ydata = ydata(hdg_orig);
                    XX = xdata';
                    YY = ydata';
                    [b,bint,r,rint,stats] = regress(YY,XX);
                    x_yoked{n} = [b(1) b(2) b(3)];
                    % Yoked-LOW:
                    yfitlow = x_yoked{n}(1)*xdata(1,1:21) + x_yoked{n}(2)*xdata(2,1:21) + x_yoked{n}(3);
                    % Yoked-HIGH:
                    yfithigh = x_yoked{n}(1)*xdata(1,22:42) + x_yoked{n}(2)*xdata(2,22:42) + x_yoked{n}(3);
                else
                    xdata = [vesminus veszero vesplus vesminus veszero vesplus ; vislowminus vislowzero vislowplus vishighminus vishighzero vishighplus];
                    xdata = xdata(:,hdg_orig); % this could be done after fitting, or before -- shouldn't change result
                    ydata = ydata(hdg_orig);
                    XX = xdata';
                    YY = ydata';
                    [b,bint,r,rint,stats] = regress(YY,XX);
                    x_yoked{n} = [b(1) b(2)];
                    % Yoked-LOW:
                    yfitlow = x_yoked{n}(1)*xdata(1,1:21) + x_yoked{n}(2)*xdata(2,1:21);
                    % Yoked-HIGH:
                    yfithigh = x_yoked{n}(1)*xdata(1,22:42) + x_yoked{n}(2)*xdata(2,22:42);
                end                    
                    
                [r,p] = corrcoef(yfitlow,ydata(1:21));
                R_corr_yoked_low(n) = r(1,2);
                P_corr_yoked_low(n) = p(1,2);
                % VAF = 1-SSE/SST (Morgan et al.)
                SSElowYok(n) = sum((yfitlow-ydata(1:21)).^2);
                SST = sum((ydata(1:21) - mean(ydata(1:21))).^2);
                VAF_yoked_low(n) = 1-SSElowYok(n)/SST;
                % adjusted R-squared
                MST = SST/(length(ydata(1:21))-1);
                MSE = SSElowYok(n)/(length(ydata(1:21))-2-1);
                R2adj_yoked_low(n) = 1-MSE/MST;

                [r,p] = corrcoef(yfithigh,ydata(22:42));
                R_corr_yoked_high(n) = r(1,2);
                P_corr_yoked_high(n) = p(1,2);
                % VAF = 1-SSE/SST (Morgan et al.)
                SSEhighYok(n) = sum((yfithigh-ydata(22:42)).^2);
                SST = sum((ydata(22:42) - mean(ydata(22:42))).^2);
                VAF_yoked_high(n) = 1-SSEhighYok(n)/SST;
                % adjusted R-squared
                MST = SST/(length(ydata(22:42))-1);
                MSE = SSEhighYok(n)/(length(ydata(22:42))-2-1);
                R2adj_yoked_high(n) = 1-MSE/MST;

                RSS_YokedWeights = SSElowYok(n) + SSEhighYok(n);
                
                %*********************************************************
                % Sequential F test
                if include_DC
                    DF_yoked = 3; % model 1
                    DF_indep = 6; % model 2
                else
                    DF_yoked = 2; % model 1
                    DF_indep = 4; % model 2
                end

                N = length(ydata);
                % N = length(horzcat(horzcat(resp{3,2,1,:},horzcat(resp{3,2,2,:}))));

                % % from Wikipedia (http://en.wikipedia.org/wiki/F-test)
                % F = (RSS_YokedWeights-RSS_IndependentWeights)/(DF_indep-DF_yoked) / (RSS_IndependentWeights/(N-DF_indep));
                % Fdist = fpdf(0:0.1:1000, DF_indep-DF_yoked, N-DF_indep);
                % p_Ftest = 1 - sum(Fdist(1:round(F*10)))/sum(Fdist); % one-tailed (see bookmark on F table)

                % % from GraphPad (http://www.graphpad.com/help/prism5/prism5help.html?howtheftestworks.htm):
                % F = ((RSS_YokedWeights-RSS_IndependentWeights)/RSS_IndependentWeights) / ((DF_indep-DF_yoked)/DF_indep);
                % Fdist = fpdf(0:0.1:1000, DF_indep-DF_yoked, DF_indep);
                % p_Ftest = 1 - sum(Fdist(1:round(F*10)))/sum(Fdist);

                % % from another site (colorado?)
                % F = ((SSER-SSE)/2) / (SSE/(N-7));

                % from Aihua (via Greg)
                F = [(RSS_YokedWeights-RSS_IndependentWeights)/(DF_indep-DF_yoked)]/[RSS_IndependentWeights/(N-DF_indep)];
                p_Ftest(n,1) = 1-fcdf(F,DF_indep-DF_yoked,N-DF_indep);

                % AIC, wikipedia:
                AIC_yoked = 2*DF_yoked + N*log(RSS_YokedWeights);
                AIC_indep = 2*DF_indep + N*log(RSS_IndependentWeights);
                % % corrected version for small N (AICc)
                AICc_yoked(n,1) = AIC_yoked + (2*DF_yoked*(DF_yoked+1)) / (N-DF_yoked-1);
                AICc_indep(n,1) = AIC_indep + (2*DF_indep*(DF_indep+1)) / (N-DF_indep-1);
                % % or
                % AIC_yoked = (N+DF_yoked)/(N-DF_yoked-2) + log(RSS_YokedWeights/N)
                % AIC_indep = (N+DF_indep)/(N-DF_indep-2) + log(RSS_IndependentWeights/N)
                % % Mathworks code library
                % AIC_yoked = (2*DF_yoked)/N + log(RSS_YokedWeights);
                % AIC_indep = (2*DF_indep)/N + log(RSS_IndependentWeights);

                
% % % %                 % Alex's sanity check: manually boost Wves/Wvis to increase
% % % %                 % overweighting and break threshold optimality
% % % %                 x_low{n}(1) = x_low{n}(1) + 0.4*abs(x_low{n}(1));
% % % %                 x_low{n}(2) = x_low{n}(2) - 0.4*abs(x_low{n}(2));              
                
            end
            
            
            %truncate means for Gaussian (resps will also be truncated later)
            if ~strcmp(noisemodel,'Poisson')
                p1_vesminus(p1_vesminus<0) = 0;
                p1_veszero(p1_veszero<0) = 0;
                p1_vesplus(p1_vesplus<0) = 0;
                p1_vislowminus(p1_vislowminus<0) = 0;
                p1_vislowzero(p1_vislowzero<0) = 0;
                p1_vislowplus(p1_vislowplus<0) = 0;
                p1_vishighminus(p1_vishighminus<0) = 0;
                p1_vishighzero(p1_vishighzero<0) = 0;
                p1_vishighplus(p1_vishighplus<0) = 0;

                %and zero sigma could also give weird results
                p2_vesminus(p2_vesminus<=0) = 0.1;
                p2_veszero(p2_veszero<=0) = 0.1;
                p2_vesplus(p2_vesplus<=0) = 0.1;
                p2_vislowminus(p2_vislowminus<=0) = 0.1;
                p2_vislowzero(p2_vislowzero<=0) = 0.1;
                p2_vislowplus(p2_vislowplus<=0) = 0.1;
                p2_vishighminus(p2_vishighminus<=0) = 0.1;
                p2_vishighzero(p2_vishighzero<=0) = 0.1;
                p2_vishighplus(p2_vishighplus<=0) = 0.1;
            end

            
%             fully flexible model:
%             R = x1*Rves + x2*Rvis + x3*Rves^2 + x4*Rvis^2 + x5*Rves*Rvis + x6
            if include_DC
                if mult_nonlin % 4 params: x3 moves to x5, x4 moves to x6, x3 and x4 set to zero
                    x_low{n}(5) = x_low{n}(3); x_low{n}(6) = x_low{n}(4); x_low{n}(3) = 0; x_low{n}(4) = 0;
                    x_high{n}(5) = x_high{n}(3); x_high{n}(6) = x_high{n}(4); x_high{n}(3) = 0; x_high{n}(4) = 0;
                elseif quad_nonlin % 5 params: x5 moves to x6, x5 set to zero
                    x_low{n}(6) = x_low{n}(5); x_low{n}(5) = 0;
                    x_high{n}(6) = x_high{n}(5); x_high{n}(5) = 0;
                else % 3 params: x3 moves to x6, x3:5 set to zero
                    x_low{n}(6) = x_low{n}(3); x_low{n}(3) = 0; x_low{n}(4) = 0; x_low{n}(5) = 0;
                    x_high{n}(6) = x_high{n}(3); x_high{n}(3) = 0; x_high{n}(4) = 0; x_high{n}(5) = 0;
                end  
            else                    
                if mult_nonlin % 3 params: x3 moves to x5; x3,4,6 set to zero
                    x_low{n}(5) = x_low{n}(3); x_low{n}(3) = 0; x_low{n}(4) = 0; x_low{n}(6) = 0;
                    x_high{n}(5) = x_high{n}(3); x_high{n}(3) = 0; x_high{n}(4) = 0; x_high{n}(6) = 0;
                elseif quad_nonlin % 4 params: x5:6 set to zero
                    x_low{n}(5) = 0; x_low{n}(6) = 0;
                    x_high{n}(5) = 0; x_high{n}(6) = 0;
                else % 2 params: x3:6 set to zero
                    x_low{n}(3:6) = 0;
                    x_high{n}(3:6) = 0;                
                end
            end

            % temp (maybe): need Fanos to make fake data the old way, that
            % is random draws from the weighted-sum-of-means instead of
            % weighting and summing the random draws themselves
            for h=1:length(X)
                fanoCombLow(h) = var(resp_trial_all_orig{n}{3,2,1,h}) / mean(resp_trial_all_orig{n}{3,2,1,h});
                fanoCombHigh(h) = var(resp_trial_all_orig{n}{3,2,2,h}) / mean(resp_trial_all_orig{n}{3,2,2,h});
            end
            FFlow = mean(fanoCombLow); FFhigh = mean(fanoCombHigh);
            
            if add_a_constant %piggy-back on unused baseline_rate code
                subtract_baseline = 1;
%                baselines(n) = -abs(randn*22)+5;
                baselines(n) = 20;
                baseline_rate = baselines(n);
            end                
            
            for h = 1:length(xx)
                resp_trial_weightmap{n,v,w}{1,2,1,h} = firing_trial{1,2,1,h}(n,:);
                resp_trial_weightmap{n,v,w}{2,2,1,h} = firing_trial{2,2,1,h}(n,:);
                resp_trial_weightmap{n,v,w}{2,2,2,h} = firing_trial{2,2,2,h}(n,:);
                
                                              % fully flexible model:
                                              % R = x1*Rves + x2*Rvis + x3*Rves^2 + x4*Rvis^2 + x5*Rves*Rvis + x6

                if strcmp(noisemodel,'Poisson')
                    resp_trial_weightmap{n,v,w}{3,1,1,h}(1:num_resp_trials) = x_low{n}(1)*poissrnd(vesminus(h),1,num_resp_trials) + ...
                                                                              x_low{n}(2)*poissrnd(vislowminus(h),1,num_resp_trials) + ...
                                                                              x_low{n}(3)*(poissrnd(vesminus(h),1,num_resp_trials)).^2 + ...
                                                                              x_low{n}(4)*(poissrnd(vislowminus(h),1,num_resp_trials)).^2 + ...
                                                                              x_low{n}(5)*(poissrnd(vesminus(h),1,num_resp_trials).*poissrnd(vislowminus(h),1,num_resp_trials)) + ...
                                                                              x_low{n}(6) + subtract_baseline*baseline_rate;
                    resp_trial_weightmap{n,v,w}{3,2,1,h}(1:num_resp_trials) = x_low{n}(1)*poissrnd(veszero(h),1,num_resp_trials) + ...
                                                                              x_low{n}(2)*poissrnd(vislowzero(h),1,num_resp_trials) + ...
                                                                              x_low{n}(3)*(poissrnd(veszero(h),1,num_resp_trials)).^2 + ...
                                                                              x_low{n}(4)*(poissrnd(vislowzero(h),1,num_resp_trials)).^2 + ...
                                                                              x_low{n}(5)*(poissrnd(veszero(h),1,num_resp_trials).*poissrnd(vislowzero(h),1,num_resp_trials)) + ...
                                                                              x_low{n}(6) + subtract_baseline*baseline_rate;
                    resp_trial_weightmap{n,v,w}{3,3,1,h}(1:num_resp_trials) = x_low{n}(1)*poissrnd(vesplus(h),1,num_resp_trials) + ...
                                                                              x_low{n}(2)*poissrnd(vislowplus(h),1,num_resp_trials) + ...
                                                                              x_low{n}(3)*(poissrnd(vesplus(h),1,num_resp_trials)).^2 + ...
                                                                              x_low{n}(4)*(poissrnd(vislowplus(h),1,num_resp_trials)).^2 + ...
                                                                              x_low{n}(5)*(poissrnd(vesplus(h),1,num_resp_trials).*poissrnd(vislowplus(h),1,num_resp_trials)) + ...
                                                                              x_low{n}(6) + subtract_baseline*baseline_rate;
                                                                          
                    resp_trial_weightmap{n,v,w}{3,1,2,h}(1:num_resp_trials) = x_high{n}(1)*poissrnd(vesminus(h),1,num_resp_trials) + ...
                                                                              x_high{n}(2)*poissrnd(vishighminus(h),1,num_resp_trials) + ...
                                                                              x_high{n}(3)*(poissrnd(vesminus(h),1,num_resp_trials)).^2 + ...
                                                                              x_high{n}(4)*(poissrnd(vishighminus(h),1,num_resp_trials)).^2 + ...
                                                                              x_high{n}(5)*(poissrnd(vesminus(h),1,num_resp_trials).*poissrnd(vishighminus(h),1,num_resp_trials)) + ...
                                                                              x_high{n}(6) + subtract_baseline*baseline_rate;
                    resp_trial_weightmap{n,v,w}{3,2,2,h}(1:num_resp_trials) = x_high{n}(1)*poissrnd(veszero(h),1,num_resp_trials) + ...
                                                                              x_high{n}(2)*poissrnd(vishighzero(h),1,num_resp_trials) + ...
                                                                              x_high{n}(3)*(poissrnd(veszero(h),1,num_resp_trials)).^2 + ...
                                                                              x_high{n}(4)*(poissrnd(vishighzero(h),1,num_resp_trials)).^2 + ...
                                                                              x_high{n}(5)*(poissrnd(veszero(h),1,num_resp_trials).*poissrnd(vishighzero(h),1,num_resp_trials)) + ...
                                                                              x_high{n}(6) + subtract_baseline*baseline_rate;
                    resp_trial_weightmap{n,v,w}{3,3,2,h}(1:num_resp_trials) = x_high{n}(1)*poissrnd(vesplus(h),1,num_resp_trials) + ...
                                                                              x_high{n}(2)*poissrnd(vishighplus(h),1,num_resp_trials) + ...
                                                                              x_high{n}(3)*(poissrnd(vesplus(h),1,num_resp_trials)).^2 + ...
                                                                              x_high{n}(4)*(poissrnd(vishighplus(h),1,num_resp_trials)).^2 + ...
                                                                              x_high{n}(5)*(poissrnd(vesplus(h),1,num_resp_trials).*poissrnd(vishighplus(h),1,num_resp_trials)) + ...
                                                                              x_high{n}(6) + subtract_baseline*baseline_rate;

%                     % return to old method
%                     tempmean = x_low{n}(1)*vesminus(h) + x_low{n}(2)*vislowminus(h) + ...
%                                x_low{n}(3)*vesminus(h)^2 + x_low{n}(4)*vislowminus(h)^2 + ...
%                                x_low{n}(5)*vesminus(h)*vislowminus(h) + x_low{n}(6) + subtract_baseline*baseline_rate;
%                     if tempmean < 1 
%                         tempmean = 1;  % or zero?
%                     end
%                     resp_trial_weightmap{n,v,w}{3,1,1,h}(1:num_resp_trials) = poissrnd(tempmean,1,num_resp_trials);
%                     tempmean = x_low{n}(1)*veszero(h) + x_low{n}(2)*vislowzero(h) + ...
%                                x_low{n}(3)*veszero(h)^2 + x_low{n}(4)*vislowzero(h)^2 + ...
%                                x_low{n}(5)*veszero(h)*vislowzero(h) + x_low{n}(6) + subtract_baseline*baseline_rate;
%                     if tempmean < 1 
%                         tempmean = 1;  % or zero?
%                     end
%                     resp_trial_weightmap{n,v,w}{3,2,1,h}(1:num_resp_trials) = poissrnd(tempmean,1,num_resp_trials);
%                     tempmean = x_low{n}(1)*vesplus(h) + x_low{n}(2)*vislowplus(h) + ...
%                                x_low{n}(3)*vesplus(h)^2 + x_low{n}(4)*vislowplus(h)^2 + ...
%                                x_low{n}(5)*vesplus(h)*vislowplus(h) + x_low{n}(6) + subtract_baseline*baseline_rate;
%                     if tempmean < 1 
%                         tempmean = 1;  % or zero?
%                     end
%                     resp_trial_weightmap{n,v,w}{3,3,1,h}(1:num_resp_trials) = poissrnd(tempmean,1,num_resp_trials);
% 
%                     tempmean = x_high{n}(1)*vesminus(h) + x_high{n}(2)*vishighminus(h) + ...
%                     x_high{n}(3)*vesminus(h)^2 + x_high{n}(4)*vishighminus(h)^2 + ...
%                     x_high{n}(5)*vesminus(h)*vishighminus(h) + x_high{n}(6) + subtract_baseline*baseline_rate;
%                     if tempmean < 1 
%                         tempmean = 1;  % or zero?
%                     end
%                     resp_trial_weightmap{n,v,w}{3,1,2,h}(1:num_resp_trials) = poissrnd(tempmean,1,num_resp_trials);
%                     tempmean = x_high{n}(1)*veszero(h) + x_high{n}(2)*vishighzero(h) + ...
%                                x_high{n}(3)*veszero(h)^2 + x_high{n}(4)*vishighzero(h)^2 + ...
%                                x_high{n}(5)*veszero(h)*vishighzero(h) + x_high{n}(6) + subtract_baseline*baseline_rate;
%                     if tempmean < 1 
%                         tempmean = 1;  % or zero?
%                     end
%                     resp_trial_weightmap{n,v,w}{3,2,2,h}(1:num_resp_trials) = poissrnd(tempmean,1,num_resp_trials);
%                     tempmean = x_high{n}(1)*vesplus(h) + x_high{n}(2)*vishighplus(h) + ...
%                                x_high{n}(3)*vesplus(h)^2 + x_high{n}(4)*vishighplus(h)^2 + ...
%                                x_high{n}(5)*vesplus(h)*vishighplus(h) + x_high{n}(6) + subtract_baseline*baseline_rate;
%                     if tempmean < 1 
%                         tempmean = 1;  % or zero?
%                     end
%                     resp_trial_weightmap{n,v,w}{3,3,2,h}(1:num_resp_trials) = poissrnd(tempmean,1,num_resp_trials);

                                                                          
                elseif strcmp(noisemodel,'Gaussian')

                    resp_trial_weightmap{n,v,w}{3,1,1,h}(1:num_resp_trials) = x_low{n}(1)*normrnd(p1_vesminus(h),p2_vesminus(h),1,num_resp_trials) + ...
                                                                              x_low{n}(2)*normrnd(p1_vislowminus(h),p2_vislowminus(h),1,num_resp_trials) + ...
                                                                              x_low{n}(3)*(normrnd(p1_vesminus(h),p2_vesminus(h),1,num_resp_trials)).^2 + ...
                                                                              x_low{n}(4)*(normrnd(p1_vislowminus(h),p2_vislowminus(h),1,num_resp_trials)).^2 + ...
                                                                              x_low{n}(5)*(normrnd(p1_vesminus(h),p2_vesminus(h),1,num_resp_trials).*normrnd(p1_vislowminus(h),p2_vislowminus(h),1,num_resp_trials)) + ...
                                                                              x_low{n}(6) + subtract_baseline*baseline_rate;
                    resp_trial_weightmap{n,v,w}{3,2,1,h}(1:num_resp_trials) = x_low{n}(1)*normrnd(p1_veszero(h),p2_veszero(h),1,num_resp_trials) + ...
                                                                              x_low{n}(2)*normrnd(p1_vislowzero(h),p2_vislowzero(h),1,num_resp_trials) + ...
                                                                              x_low{n}(3)*(normrnd(p1_veszero(h),p2_veszero(h),1,num_resp_trials)).^2 + ...
                                                                              x_low{n}(4)*(normrnd(p1_vislowzero(h),p2_vislowzero(h),1,num_resp_trials)).^2 + ...
                                                                              x_low{n}(5)*(normrnd(p1_veszero(h),p2_veszero(h),1,num_resp_trials).*normrnd(p1_vislowzero(h),p2_vislowzero(h),1,num_resp_trials)) + ...
                                                                              x_low{n}(6) + subtract_baseline*baseline_rate;
                    resp_trial_weightmap{n,v,w}{3,3,1,h}(1:num_resp_trials) = x_low{n}(1)*normrnd(p1_vesplus(h),p2_vesplus(h),1,num_resp_trials) + ...
                                                                              x_low{n}(2)*normrnd(p1_vislowplus(h),p2_vislowplus(h),1,num_resp_trials) + ...
                                                                              x_low{n}(3)*(normrnd(p1_vesplus(h),p2_vesplus(h),1,num_resp_trials)).^2 + ...
                                                                              x_low{n}(4)*(normrnd(p1_vislowplus(h),p2_vislowplus(h),1,num_resp_trials)).^2 + ...
                                                                              x_low{n}(5)*(normrnd(p1_vesplus(h),p2_vesplus(h),1,num_resp_trials).*normrnd(p1_vislowplus(h),p2_vislowplus(h),1,num_resp_trials)) + ...
                                                                              x_low{n}(6) + subtract_baseline*baseline_rate;
                                                                          
                    resp_trial_weightmap{n,v,w}{3,1,2,h}(1:num_resp_trials) = x_high{n}(1)*normrnd(p1_vesminus(h),p2_vesminus(h),1,num_resp_trials) + ...
                                                                              x_high{n}(2)*normrnd(p1_vishighminus(h),p2_vishighminus(h),1,num_resp_trials) + ...
                                                                              x_high{n}(3)*(normrnd(p1_vesminus(h),p2_vesminus(h),1,num_resp_trials)).^2 + ...
                                                                              x_high{n}(4)*(normrnd(p1_vishighminus(h),p2_vishighminus(h),1,num_resp_trials)).^2 + ...
                                                                              x_high{n}(5)*(normrnd(p1_vesminus(h),p2_vesminus(h),1,num_resp_trials).*normrnd(p1_vishighminus(h),p2_vishighminus(h),1,num_resp_trials)) + ...
                                                                              x_high{n}(6) + subtract_baseline*baseline_rate;
                    resp_trial_weightmap{n,v,w}{3,2,2,h}(1:num_resp_trials) = x_high{n}(1)*normrnd(p1_veszero(h),p2_veszero(h),1,num_resp_trials) + ...
                                                                              x_high{n}(2)*normrnd(p1_vishighzero(h),p2_vishighzero(h),1,num_resp_trials) + ...
                                                                              x_high{n}(3)*(normrnd(p1_veszero(h),p2_veszero(h),1,num_resp_trials)).^2 + ...
                                                                              x_high{n}(4)*(normrnd(p1_vishighzero(h),p2_vishighzero(h),1,num_resp_trials)).^2 + ...
                                                                              x_high{n}(5)*(normrnd(p1_veszero(h),p2_veszero(h),1,num_resp_trials).*normrnd(p1_vishighzero(h),p2_vishighzero(h),1,num_resp_trials)) + ...
                                                                              x_high{n}(6) + subtract_baseline*baseline_rate;
                    resp_trial_weightmap{n,v,w}{3,3,2,h}(1:num_resp_trials) = x_high{n}(1)*normrnd(p1_vesplus(h),p2_vesplus(h),1,num_resp_trials) + ...
                                                                              x_high{n}(2)*normrnd(p1_vishighplus(h),p2_vishighplus(h),1,num_resp_trials) + ...
                                                                              x_high{n}(3)*(normrnd(p1_vesplus(h),p2_vesplus(h),1,num_resp_trials)).^2 + ...
                                                                              x_high{n}(4)*(normrnd(p1_vishighplus(h),p2_vishighplus(h),1,num_resp_trials)).^2 + ...
                                                                              x_high{n}(5)*(normrnd(p1_vesplus(h),p2_vesplus(h),1,num_resp_trials).*normrnd(p1_vishighplus(h),p2_vishighplus(h),1,num_resp_trials)) + ...
                                                                              x_high{n}(6) + subtract_baseline*baseline_rate;

%                     % return to old method
%                     tempmean = x_low{n}(1)*vesminus(h) + x_low{n}(2)*vislowminus(h) + ...
%                                x_low{n}(3)*vesminus(h)^2 + x_low{n}(4)*vislowminus(h)^2 + ...
%                                x_low{n}(5)*vesminus(h)*vislowminus(h) + x_low{n}(6) + subtract_baseline*baseline_rate;
%                     tempvar = tempmean*FFlow;
%                     resp_trial_weightmap{n,v,w}{3,1,1,h}(1:num_resp_trials) = normrnd(tempmean,sqrt(tempvar),1,num_resp_trials);
%                     tempmean = x_low{n}(1)*veszero(h) + x_low{n}(2)*vislowzero(h) + ...
%                                x_low{n}(3)*veszero(h)^2 + x_low{n}(4)*vislowzero(h)^2 + ...
%                                x_low{n}(5)*veszero(h)*vislowzero(h) + x_low{n}(6) + subtract_baseline*baseline_rate;
%                     tempvar = tempmean*FFlow;
%                     resp_trial_weightmap{n,v,w}{3,2,1,h}(1:num_resp_trials) = normrnd(tempmean,sqrt(tempvar),1,num_resp_trials);
%                     tempmean = x_low{n}(1)*vesplus(h) + x_low{n}(2)*vislowplus(h) + ...
%                                x_low{n}(3)*vesplus(h)^2 + x_low{n}(4)*vislowplus(h)^2 + ...
%                                x_low{n}(5)*vesplus(h)*vislowplus(h) + x_low{n}(6) + subtract_baseline*baseline_rate;
%                     tempvar = tempmean*FFlow;
%                     resp_trial_weightmap{n,v,w}{3,3,1,h}(1:num_resp_trials) = normrnd(tempmean,sqrt(tempvar),1,num_resp_trials);
% 
%                     tempmean = x_high{n}(1)*vesminus(h) + x_high{n}(2)*vishighminus(h) + ...
%                     x_high{n}(3)*vesminus(h)^2 + x_high{n}(4)*vishighminus(h)^2 + ...
%                     x_high{n}(5)*vesminus(h)*vishighminus(h) + x_high{n}(6) + subtract_baseline*baseline_rate;
%                     tempvar = tempmean*FFhigh;
%                     resp_trial_weightmap{n,v,w}{3,1,2,h}(1:num_resp_trials) = normrnd(tempmean,sqrt(tempvar),1,num_resp_trials);
%                     tempmean = x_high{n}(1)*veszero(h) + x_high{n}(2)*vishighzero(h) + ...
%                                x_high{n}(3)*veszero(h)^2 + x_high{n}(4)*vishighzero(h)^2 + ...
%                                x_high{n}(5)*veszero(h)*vishighzero(h) + x_high{n}(6) + subtract_baseline*baseline_rate;
%                     tempvar = tempmean*FFhigh;
%                     resp_trial_weightmap{n,v,w}{3,2,2,h}(1:num_resp_trials) = normrnd(tempmean,sqrt(tempvar),1,num_resp_trials);
%                     tempmean = x_high{n}(1)*vesplus(h) + x_high{n}(2)*vishighplus(h) + ...
%                                x_high{n}(3)*vesplus(h)^2 + x_high{n}(4)*vishighplus(h)^2 + ...
%                                x_high{n}(5)*vesplus(h)*vishighplus(h) + x_high{n}(6) + subtract_baseline*baseline_rate;
%                     tempvar = tempmean*FFhigh;
%                     resp_trial_weightmap{n,v,w}{3,3,2,h}(1:num_resp_trials) = normrnd(tempmean,sqrt(tempvar),1,num_resp_trials);

                elseif strcmp(noisemodel,'Gamma');
                    resp_trial_weightmap{n,v,w}{3,1,1,h}(1:num_resp_trials) = x_low{n}(1)*gamrnd(p1_vesminus(h),p2_vesminus(h),1,num_resp_trials) + x_low{n}(2)*gamrnd(p1_vislowminus(h),p2_vislowminus(h),1,num_resp_trials) + x_low{n}(3) + subtract_baseline*baseline_rate;
                    resp_trial_weightmap{n,v,w}{3,2,1,h}(1:num_resp_trials) = x_low{n}(1)*gamrnd(p1_veszero(h),p2_veszero(h),1,num_resp_trials) + x_low{n}(2)*gamrnd(p1_vislowzero(h),p2_vislowzero(h),1,num_resp_trials) + x_low{n}(3) + subtract_baseline*baseline_rate;
                    resp_trial_weightmap{n,v,w}{3,3,1,h}(1:num_resp_trials) = x_low{n}(1)*gamrnd(p1_vesplus(h),p2_vesplus(h),1,num_resp_trials) + x_low{n}(2)*gamrnd(p1_vislowplus(h),p2_vislowplus(h),1,num_resp_trials) + x_low{n}(3) + subtract_baseline*baseline_rate;

                    resp_trial_weightmap{n,v,w}{3,1,2,h}(1:num_resp_trials) = x_high{n}(1)*gamrnd(p1_vesminus(h),p2_vesminus(h),1,num_resp_trials) + x_high{n}(2)*gamrnd(p1_vishighminus(h),p2_vishighminus(h),1,num_resp_trials) + x_high{n}(3) + subtract_baseline*baseline_rate;
                    resp_trial_weightmap{n,v,w}{3,2,2,h}(1:num_resp_trials) = x_high{n}(1)*gamrnd(p1_veszero(h),p2_veszero(h),1,num_resp_trials) + x_high{n}(2)*gamrnd(p1_vishighzero(h),p2_vishighzero(h),1,num_resp_trials) + x_high{n}(3) + subtract_baseline*baseline_rate;
                    resp_trial_weightmap{n,v,w}{3,3,2,h}(1:num_resp_trials) = x_high{n}(1)*gamrnd(p1_vesplus(h),p2_vesplus(h),1,num_resp_trials) + x_high{n}(2)*gamrnd(p1_vishighplus(h),p2_vishighplus(h),1,num_resp_trials) + x_high{n}(3) + subtract_baseline*baseline_rate;                    
                end
                % for nonpoisson dists (really only Gaussian) must truncate negative responses
                if ~strcmp(noisemodel,'Poisson')
                    resp_trial_weightmap{n,v,w}{3,1,1,h}(resp_trial_weightmap{n,v,w}{3,1,1,h} < 0) = 0;
                    resp_trial_weightmap{n,v,w}{3,2,1,h}(resp_trial_weightmap{n,v,w}{3,2,1,h} < 0) = 0;
                    resp_trial_weightmap{n,v,w}{3,3,1,h}(resp_trial_weightmap{n,v,w}{3,3,1,h} < 0) = 0;

                    resp_trial_weightmap{n,v,w}{3,1,2,h}(resp_trial_weightmap{n,v,w}{3,1,2,h} < 0) = 0;
                    resp_trial_weightmap{n,v,w}{3,2,2,h}(resp_trial_weightmap{n,v,w}{3,2,2,h} < 0) = 0;
                    resp_trial_weightmap{n,v,w}{3,3,2,h}(resp_trial_weightmap{n,v,w}{3,3,2,h} < 0) = 0;
                end                
            end

            % store fake data with desired weights (comb), and original singles, as resp_mean_weightmap
            resp_mean_weightmap{n,v,w}{1,2,1} = firing{1,2,1}(n,:);
            resp_mean_weightmap{n,v,w}{2,2,1} = firing{2,2,1}(n,:);
            resp_mean_weightmap{n,v,w}{2,2,2} = firing{2,2,2}(n,:);
            for i = 1:length(xx)
                resp_mean_weightmap{n,v,w}{3,1,1}(i) = mean(resp_trial_weightmap{n,v,w}{3,1,1,i});
                resp_mean_weightmap{n,v,w}{3,2,1}(i) = mean(resp_trial_weightmap{n,v,w}{3,2,1,i});
                resp_mean_weightmap{n,v,w}{3,3,1}(i) = mean(resp_trial_weightmap{n,v,w}{3,3,1,i});

                resp_mean_weightmap{n,v,w}{3,1,2}(i) = mean(resp_trial_weightmap{n,v,w}{3,1,2,i});
                resp_mean_weightmap{n,v,w}{3,2,2}(i) = mean(resp_trial_weightmap{n,v,w}{3,2,2,i});
                resp_mean_weightmap{n,v,w}{3,3,2}(i) = mean(resp_trial_weightmap{n,v,w}{3,3,2,i});
            end

            
%             % new alex request: compute thresholds (ROC?) from Morgan vs. Optimal weights
%             % LOW COH
%             [rr,pp] = corrcoef(xx, resp_mean_weightmap{n,v,w}{3,2,1});
%             line_r = rr(1,2);
%             for i = 1:length(xx)
%                 if line_r > 0
%                     Neuro_correct_low(i) = rocN(resp_trial_weightmap{n,v,w}{3,2,1,i}, resp_trial_weightmap{n,v,w}{3,2,1,26}, 100);
%                 else 
%                     Neuro_correct_low(i) = rocN(resp_trial_weightmap{n,v,w}{3,2,1,26}, resp_trial_weightmap{n,v,w}{3,2,1,i}, 100);
%                 end
%             end
%             fit_data_neuro_low(:,1) = xx';
%             fit_data_neuro_low(:,2) = Neuro_correct_low';
%             fit_data_neuro_low(:,3) = num_resp_trials;
%             [bbb,ttt] = cum_gaussfit_max1(fit_data_neuro_low);
%             Thresh_neu_low(n) = ttt;
%             
%             % HIGH COH
%             [rr,pp] = corrcoef(xx, resp_mean_weightmap{n,v,w}{3,2,2});
%             line_r = rr(1,2);
%             for i = 1:length(xx)
%                 if line_r > 0
%                     Neuro_correct_high(i) = rocN(resp_trial_weightmap{n,v,w}{3,2,2,i}, resp_trial_weightmap{n,v,w}{3,2,2,26}, 100);
%                 else 
%                     Neuro_correct_high(i) = rocN(resp_trial_weightmap{n,v,w}{3,2,2,26}, resp_trial_weightmap{n,v,w}{3,2,2,i}, 100);
%                 end
%             end
%             fit_data_neuro_high(:,1) = xx';
%             fit_data_neuro_high(:,2) = Neuro_correct_high';
%             fit_data_neuro_high(:,3) = num_resp_trials;
%             [bbb,ttt] = cum_gaussfit_max1(fit_data_neuro_high);
%             Thresh_neu_high(n) = ttt;

        end
        
%         centers = [-20 -10 -5 -2.5 -1.25 -.625 -.3125];
%         centers = [centers 0 fliplr(centers)*-1];
%         figure(1919);
%         subplot(2,1,1);
%         bar(hist(Wratio_opt_low(:,100),centers)); set(gca,'XTickLabel',centers)
%         subplot(2,1,2);
%         bar(hist(Wratio_opt_high(:,100),centers)); set(gca,'XTickLabel',centers)
        
        if strcmp(weights_to_apply,'NewMorgan')
            for n = 1:cell_num
                data1(n,:) = [x_low{n}' x_high{n}'];
            end
            data2 = [(R_corr_low').^2 P_corr_low' (R_corr_high').^2 P_corr_high'];
            data3 = [VAF_low' VAF_high'];
            data4 = [R2adj_low' R2adj_high'];
            data_yoked_vs_indep = [p_Ftest AICc_yoked AICc_indep];
%             asdf
        end
        
%         data5 = [Thresh_neu_low' Thresh_neu_high'];
        
%         figure; loglog(tLow(:,2),tLow(:,1),'bo'); xlim([1 100]); ylim([1 100]); axis square;
%         hold on; loglog(tHigh(:,2),tHigh(:,1),'ro',[1 100],[1 100],'k--'); legend('Low Coh','High Coh');
%         xlabel('Neuronal threshold w/ optimal weights'); ylabel('Neuronal threshold w/ actual weights');

%         save(['C:\Documents and Settings\Chris\My Documents\MATLAB\' timefolder(1:end-1) '_Morgans.mat'],'data1','data2','data3','data4')

%         figure; plot(X,mean(residualErrLow));
%         figure; plot(X,mean(residualErrHigh));
%         pause
    end
end

% the weightmap thing is a remnant from the contour plots -- if desired,
% can use indices v and w to cycle through combinations of weights -- now
% they are just 1 and 1 
clear resp_mean_all resp_trial_all
for n = 1:cell_num
    resp_mean_all{n} = resp_mean_weightmap{n,1,1};
    resp_trial_all{n} = resp_trial_weightmap{n,1,1};
end

% enforce goodness of fit criterion:
if strcmp(weights_to_apply,'NewMorgan')
    cells_to_delete = find(P_corr_low>0.05 | P_corr_high>0.05);
    resp_mean_all(cells_to_delete) = [];
    resp_trial_all(cells_to_delete) = [];
    cell_num = length(resp_mean_all)
end

% ------------------------------------------------------------------
else
    if ~keep_real_trials
        % now convert both means and single-trial responses back to my format
        clear resp_trial_all resp_mean_all  % not sure if I have to clear these
        for n = 1:cell_num 
            % convert back from Yong's format to mine
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

    end
    
end % end use_fake_comb

% % % 
% % % % compute weighted sums of *original* spike counts, to test variance etc. (CANNOT BE DONE FOR SIMULATED CONFLICT CONDITIONS!)
% % % xx = X;
% % % for n = 1:cell_num
% % %     for h = 1:length(xx)
% % %         % _ORIG OR NOT?!  DEPENDS ON THE QUESTION: I think YES
% % %         % for 2nd alex request (fano factors)
% % %         maxreps = min([length(resp_trial_all_orig{n}{1,2,1,h}) length(resp_trial_all_orig{n}{2,2,1,h})]);
% % %         resp_trial_comb_weightedSumOfTrials{n}{3,2,1,h} = x_low{n}(1)*resp_trial_all_orig{n}{1,2,1,h}(1:maxreps) + x_low{n}(2)*resp_trial_all_orig{n}{2,2,1,h}(1:maxreps) + subtract_baseline*baseline_rate;
% % %         maxreps = min([length(resp_trial_all_orig{n}{1,2,1,h}) length(resp_trial_all_orig{n}{2,2,2,h})]);
% % %         resp_trial_comb_weightedSumOfTrials{n}{3,2,2,h} = x_high{n}(1)*resp_trial_all_orig{n}{1,2,1,h}(1:maxreps) + x_high{n}(2)*resp_trial_all_orig{n}{2,2,2,h}(1:maxreps) + subtract_baseline*baseline_rate;
% % % 
% % %         % COMPARE VARIANCE OF OPTIMAL-WEIGHTED-SUM TRIALS VS. ACTUAL COMB TRIALS
% % %         varRealLow(n,h) = var(resp_trial_all_orig{n}{3,2,1,h});
% % %             meanRealLow(n,h) = mean(resp_trial_all_orig{n}{3,2,1,h});
% % %         varFakeLow(n,h) = var(resp_trial_comb_weightedSumOfTrials{n}{3,2,1,h});
% % %         varRealHigh(n,h) = var(resp_trial_all_orig{n}{3,2,2,h});
% % %             meanRealHigh(n,h) = mean(resp_trial_all_orig{n}{3,2,2,h});
% % %         varFakeHigh(n,h) = var(resp_trial_comb_weightedSumOfTrials{n}{3,2,2,h});
% % % 
% % %         % Request by Alex: make Poisson draws starting with
% % %         % Morgan-weighted *means*, and compare it to real variance:
% % %         % real variance should be greater if Alex is correct
% % %         % that combined cannot be modeled by adding single-cue
% % %         % means and then Poisson sampling [for this DO NOT
% % %         % replace single-cues with their linear fits, above]
% % %         num_resp_trials = 13;
% % %         resp_trial_comb_drawFromWeightedMeans{n}{3,2,1,h} = poissrnd(x_low{n}(1)*resp_mean_all{n}{1,2,1}(h) + x_low{n}(2)*resp_mean_all{n}{2,2,1}(h),1,num_resp_trials); 
% % %         resp_trial_comb_drawFromWeightedMeans{n}{3,2,2,h} = poissrnd(x_high{n}(1)*resp_mean_all{n}{1,2,1}(h) + x_high{n}(2)*resp_mean_all{n}{2,2,2}(h),1,num_resp_trials); 
% % %         varAlexTestLow(n,h) = var(resp_trial_comb_drawFromWeightedMeans{n}{3,2,1,h});
% % %         varAlexTestHigh(n,h) = var(resp_trial_comb_drawFromWeightedMeans{n}{3,2,2,h});
% % % 
% % %         % then compare Fano factors between single-cues and
% % %         % fake comb (as weighted sums of single-cue *trials*):
% % %         % should be different, given that (a) variance differs
% % %         % between real comb and fake Morgan-mean comb
% % %         % (AlexTest), and (b) variance does not differ between
% % %         % real comb and fake morgan-weighted *trials*
% % %         fanoVes(n,h) = var(resp_trial_all_orig{n}{1,2,1,h}) / mean(resp_trial_all_orig{n}{1,2,1,h});
% % %         fanoVisLow(n,h) = var(resp_trial_all_orig{n}{2,2,1,h}) / mean(resp_trial_all_orig{n}{2,2,1,h});
% % %         fanoVisHigh(n,h) = var(resp_trial_all_orig{n}{2,2,2,h}) / mean(resp_trial_all_orig{n}{2,2,2,h});
% % %         fanoCombLowFake(n,h) = varFakeLow(n,h) / mean(resp_trial_comb_weightedSumOfTrials{n}{3,2,1,h});
% % %         fanoCombHighFake(n,h) = varFakeHigh(n,h) / mean(resp_trial_comb_weightedSumOfTrials{n}{3,2,2,h});                    
% % %         % check real comb again too:
% % %         fanoCombLowReal(n,h) = var(resp_trial_all_orig{n}{3,2,1,h}) / mean(resp_trial_all_orig{n}{3,2,1,h});
% % %         fanoCombHighReal(n,h) = var(resp_trial_all_orig{n}{3,2,2,h}) / mean(resp_trial_all_orig{n}{3,2,2,h});
% % %     end
% % % end
% % % 
% % % % varToPlotVsRealsLow = varFakeLow; varToPlotVsRealsHigh = varFakeHigh;
% % % 
% % % FanoLow = varRealLow./meanRealLow;
% % % FanoHigh = varRealHigh./meanRealHigh;
% % % varToPlotVsRealsLow = varAlexTestLow.*FanoLow; varToPlotVsRealsHigh = varAlexTestHigh.*FanoHigh;
% % % 
% % % figure; subplot(2,1,1); set(gcf,'Position',[500 20 700 900]);
% % % loglog(reshape(varToPlotVsRealsLow,numel(varToPlotVsRealsLow),1), reshape(varRealLow,numel(varRealLow),1), 'bo', [1 1e3], [1 1e3], 'k--');
% % % xlim([1 1e3]); ylim([1 1e3]); axis square;
% % % % xlabel('Variance of optimal weighted sum of single-cue trials');
% % % % xlabel('Variance of Morgan-weighted sum of single-cue trials');
% % % xlabel('Variance of poisson-drawn *means* (Morgan-weighted)');
% % % ylabel('Spike count variance of actual combined responses');
% % % title('Real vs. Simulated Combined Response Variance: Low Coherence');
% % % P = signrank(reshape(varToPlotVsRealsLow,numel(varToPlotVsRealsLow),1), reshape(varRealLow,numel(varRealLow),1));
% % % text(1.25, 600, ['Wilc. matched pairs: p = ' num2str(P)]);
% % % 
% % % subplot(2,1,2);
% % % loglog(reshape(varToPlotVsRealsHigh,numel(varToPlotVsRealsHigh),1), reshape(varRealHigh,numel(varRealHigh),1), 'ro', [1 1e3], [1 1e3], 'k--');
% % % xlim([1 1e3]); ylim([1 1e3]); axis square;
% % % % xlabel('Variance of optimal weighted sum of single-cue trials');
% % % % xlabel('Variance of Morgan-weighted sum of single-cue trials');
% % % xlabel('Variance of poisson-drawn *means* (Morgan-weighted)');
% % % ylabel('Spike count variance of actual combined responses');
% % % title('Real vs. Simulated Combined Response Variance: High Coherence');
% % % P = signrank(reshape(varToPlotVsRealsHigh,numel(varToPlotVsRealsHigh),1), reshape(varRealHigh,numel(varRealHigh),1));
% % % text(1.25, 600, ['Wilc. matched pairs: p = ' num2str(P)]);
% % % 
% % % % fano factors:
% % % figure; subplot(2,1,1); set(gcf,'Position',[500 20 700 900]);
% % % loglog(reshape(fanoVes,numel(fanoVes),1), reshape(fanoCombLowFake,numel(fanoCombLowFake),1), 'ko', [1e-1 1e2], [1e-1 1e2], 'k--');
% % % xlim([1e-1 1e2]); ylim([1e-1 1e2]); axis square; title('Low Coherence');
% % % xlabel('Fano factor - real vestib'); ylabel('Fano factor -- fake combined, weighed sum of spk cts');
% % % P = signrank(reshape(fanoVes,numel(fanoVes),1), reshape(fanoCombLowFake,numel(fanoCombLowFake),1));
% % % text(.15, 70, ['Wilc. matched pairs: p = ' num2str(P)]);
% % % subplot(2,1,2);
% % % loglog(reshape(fanoVisLow,numel(fanoVisLow),1), reshape(fanoCombLowFake,numel(fanoCombLowFake),1), 'mo', [1e-1 1e2], [1e-1 1e2], 'k--');
% % % xlim([1e-1 1e2]); ylim([1e-1 1e2]); axis square; title('Low Coherence');
% % % xlabel('Fano factor - real visual'); ylabel('Fano factor -- fake combined, weighed sum of spk cts');
% % % P = signrank(reshape(fanoVisLow,numel(fanoVisLow),1), reshape(fanoCombLowFake,numel(fanoCombLowFake),1));
% % % text(.15, 70, ['Wilc. matched pairs: p = ' num2str(P)]);
% % % 
% % % figure; subplot(2,1,1); set(gcf,'Position',[550 20 700 900]);
% % % loglog(reshape(fanoVes,numel(fanoVes),1), reshape(fanoCombHighFake,numel(fanoCombHighFake),1), 'ko', [1e-1 1e2], [1e-1 1e2], 'k--');
% % % xlim([1e-1 1e2]); ylim([1e-1 1e2]); axis square; title('High Coherence');
% % % xlabel('Fano factor - real vestib'); ylabel('Fano factor -- fake combined, weighed sum of spk cts');
% % % P = signrank(reshape(fanoVes,numel(fanoVes),1), reshape(fanoCombHighFake,numel(fanoCombHighFake),1));
% % % text(.15, 70, ['Wilc. matched pairs: p = ' num2str(P)]);
% % % subplot(2,1,2);
% % % loglog(reshape(fanoVisHigh,numel(fanoVisHigh),1), reshape(fanoCombHighFake,numel(fanoCombHighFake),1), 'ro', [1e-1 1e2], [1e-1 1e2], 'k--');
% % % xlim([1e-1 1e2]); ylim([1e-1 1e2]); axis square; title('High Coherence');
% % % xlabel('Fano factor - real visual'); ylabel('Fano factor -- fake combined, weighed sum of spk cts');
% % % P = signrank(reshape(fanoVisHigh,numel(fanoVisHigh),1), reshape(fanoCombHighFake,numel(fanoCombHighFake),1));
% % % text(.15, 70, ['Wilc. matched pairs: p = ' num2str(P)]);
% % % 
% % % figure; subplot(2,1,1); set(gcf,'Position',[500 20 700 900]);
% % % loglog(reshape(fanoVes,numel(fanoVes),1), reshape(fanoCombLowReal,numel(fanoCombLowReal),1), 'ko', [1e-1 1e2], [1e-1 1e2], 'k--');
% % % xlim([1e-1 1e2]); ylim([1e-1 1e2]); axis square; title('Fano Comparison - Low Coherence');
% % % xlabel('Fano factor - real vestib'); ylabel('Fano factor -- real combined');
% % % P = signrank(reshape(fanoVes,numel(fanoVes),1), reshape(fanoCombLowReal,numel(fanoCombLowReal),1));
% % % text(.15, 70, ['Wilc. matched pairs: p = ' num2str(P)]);
% % % subplot(2,1,2);
% % % loglog(reshape(fanoVisLow,numel(fanoVisLow),1), reshape(fanoCombLowReal,numel(fanoCombLowReal),1), 'mo', [1e-1 1e2], [1e-1 1e2], 'k--');
% % % xlim([1e-1 1e2]); ylim([1e-1 1e2]); axis square; title('Fano Comparison - Low Coherence');
% % % xlabel('Fano factor - real visual'); ylabel('Fano factor -- real combined');
% % % P = signrank(reshape(fanoVisLow,numel(fanoVisLow),1), reshape(fanoCombLowReal,numel(fanoCombLowReal),1));
% % % text(.15, 70, ['Wilc. matched pairs: p = ' num2str(P)]);
% % % 
% % % figure; subplot(2,1,1); set(gcf,'Position',[550 20 700 900]);
% % % loglog(reshape(fanoVes,numel(fanoVes),1), reshape(fanoCombHighReal,numel(fanoCombHighReal),1), 'ko', [1e-1 1e2], [1e-1 1e2], 'k--');
% % % xlim([1e-1 1e2]); ylim([1e-1 1e2]); axis square; title('Fano Comparison - High Coherence');
% % % xlabel('Fano factor - real vestib'); ylabel('Fano factor -- real combined');
% % % P = signrank(reshape(fanoVes,numel(fanoVes),1), reshape(fanoCombHighReal,numel(fanoCombHighReal),1));
% % % text(.15, 70, ['Wilc. matched pairs: p = ' num2str(P)]);
% % % subplot(2,1,2);
% % % loglog(reshape(fanoVisHigh,numel(fanoVisHigh),1), reshape(fanoCombHighReal,numel(fanoCombHighReal),1), 'ro', [1e-1 1e2], [1e-1 1e2], 'k--');
% % % xlim([1e-1 1e2]); ylim([1e-1 1e2]); axis square; title('Fano Comparison - High Coherence');
% % % xlabel('Fano factor - real visual'); ylabel('Fano factor -- real combined');
% % % P = signrank(reshape(fanoVisHigh,numel(fanoVisHigh),1), reshape(fanoCombHighReal,numel(fanoCombHighReal),1));
% % % text(.15, 70, ['Wilc. matched pairs: p = ' num2str(P)]);



% %------------------------------------------------
% % before running the decoder, just check whether Fisher info is optimal
% % (no information lost; i.e. equal to the sum of single-cue FI)
% for n = 1:cell_num
%     fprime_ves = diff(resp_mean_all{n}{1,2,1}(13:39));
%     fprime_vislow = diff(resp_mean_all{n}{2,2,1}(13:39));
%     fprime_vishigh = diff(resp_mean_all{n}{2,2,2}(13:39));
%     fprime_comblow = diff(resp_mean_all{n}{3,2,1}(13:39));
%     fprime_combhigh = diff(resp_mean_all{n}{3,2,2}(13:39));
% 
%     fisher_ves(n,:) = fprime_ves.^2 ./ resp_mean_all{n}{1,2,1}(13:38);
%     fisher_vislow(n,:) = fprime_vislow.^2 ./ resp_mean_all{n}{2,2,1}(13:38);
%     fisher_vishigh(n,:) = fprime_vishigh.^2 ./ resp_mean_all{n}{2,2,2}(13:38);
% %     fisher_comblow(n,:) = fprime_comblow.^2 ./ resp_mean_all{n}{3,2,1}(13:38);
% %     fisher_combhigh(n,:) = fprime_combhigh.^2 ./ resp_mean_all{n}{3,2,2}(13:38);
%     
%     % can't use simplified FI expression for comb! (nonpoisson)
% %     fisher_comblow(n) = resp_trial_all{n}{3,2,1,26};
% %     fisher_combhigh(n) = fprime_combhigh.^2 / std(resp_trial_all{n}{3,2,2,26})^2;   
%     
%     like = zeros(1,length(heading));
%     for h = 1:length(heading)
%         [counts,binctr] = hist(resp_trial_all{n}{s,2,k,h},30);
%         countProb = counts/max(counts);
%         findBestBin = find(abs(binctr-r(n)) == min(abs(binctr-r(n))));
%         like(h) = countProb(findBestBin(1));
%     end
% %         figure; plot(heading,lik);
%     temp_L(n,:) = interp1(heading,like,hdg,'linear');
%     
%     
%     figure(12); hold on; axis square;
% %     plot(fisher_ves(n,14) + fisher_vislow(n,14), fisher_comblow(n,14),'bo',[0 2e-3], [0 2e-3], 'k--'); 
%     plot(fisher_ves(n,14) + fisher_vishigh(n,14), fisher_combhigh(n,14),'ro',[0 .012], [0 .012], 'k--');
% %     figure(13); hold on; axis square;
% %     plot(mean(fisher_ves(n,:)) + mean(fisher_vislow(n,:)), mean(fisher_comblow(n,:)),'bo',[0 2e-3], [0 2e-3], 'k--'); 
% %     plot(mean(fisher_ves(n,:)) + mean(fisher_vishigh(n,:)), mean(fisher_combhigh(n,:)),'ro');
% 
% %     % sanity check: variance needs to equal mean for singles (it does)
% %     figure(15); hold on; axis square;
% %     plot(sqrt(resp_mean_all{n}{1,2,1}(26)),std(resp_trial_all{n}{1,2,1,26}),'ko');
% %     plot(sqrt(resp_mean_all{n}{2,2,1}(26)),std(resp_trial_all{n}{2,2,1,26}),'ko');
% % 	  plot(sqrt(resp_mean_all{n}{2,2,2}(26)),std(resp_trial_all{n}{2,2,2,26}),'ko');
% %     pause;
%     
% %     figure(14); clf;
% % %     plot(xx,resp_mean_all{n}{1,2,1},'k',xx,resp_mean_all{n}{2,2,1},'m',xx,resp_mean_all{n}{2,2,2},'r',xx,resp_mean_all{n}{3,2,1},'c',xx,resp_mean_all{n}{3,2,2},'b')
% %     errorbar(xx,resp_mean_all{n}{1,2,1},[zeros(1,10) std(resp_trial_all{n}{1,2,1,26}) zeros(1,40)], 'k'); hold on;
% %     errorbar(xx,resp_mean_all{n}{2,2,1},[zeros(1,20) std(resp_trial_all{n}{2,2,1,26}) zeros(1,30)], 'm'); 
% %     errorbar(xx,resp_mean_all{n}{2,2,2},[zeros(1,25) std(resp_trial_all{n}{2,2,2,26}) zeros(1,25)], 'r');
% %     errorbar(xx,resp_mean_all{n}{3,2,1},[zeros(1,30) std(resp_trial_all{n}{3,2,1,26}) zeros(1,20)], 'c');
% %     errorbar(xx,resp_mean_all{n}{3,2,2},[zeros(1,40) std(resp_trial_all{n}{3,2,2,26}) zeros(1,10)], 'b');
% % %     title(['FIsumLow = ' num2str(mean(fisher_ves(n,:)) + mean(fisher_vislow(n,:))) '  FIcombLow = ' num2str(mean(fisher_comblow(n,:)))]);
% % 	title(['FIsumHigh = ' num2str(mean(fisher_ves(n,:)) + mean(fisher_vishigh(n,:))) ' FIcombHigh = ' num2str(mean(fisher_combhigh(n,:)))]); 
% %     xlabel(['Ves/vis weight ratio (low,high) = ' num2str(Wratio_opt{n}(1)) ' ' num2str(Wratio_opt{n}(2))]);
% %     pause;
% 
% end


% n=21 %15 16 21 22 36
% figure; plot(unique_heading,resp_mean_all{n}{1,2,1},'k-',unique_heading,resp_mean_all{n}{2,2,1},'m-',unique_heading,resp_mean_all{n}{2,2,2},'r-',unique_heading,resp_mean_all{n}{3,2,1},'c-',unique_heading,resp_mean_all{n}{3,2,2},'b-');

% % % TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
% % for s = 1:3; % stim type
% %     for j = 1:3 % conflict angle
% %         for k = 1:2 % coherence
% %             param1_interp{s,j,k} = [];
% %             param2_interp{s,j,k} = [];
% %         end
% %     end
% % end
% % % TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP


% regardless of real or fake, safe bet is to re-fit (comb only, really) to get
% proper params to pass in, in case want to use fit for likelihood:
if ~keep_real_trials  %(this is only necessary if interpolating)
if strcmp(likemodel,'Gaussian')
    for n = 1:cell_num  
        n
        for s = 1:3; % stim type
            for j = 1:3 % conflict angle
                for k = 1:2 % coherence
                    if ~isempty(resp_mean_all{n}{s,j,k})
                        for i = 1:length(XI)
                            [muhat,sigmahat] = normfit(resp_trial_all{n}{s,j,k,i});
                            if muhat == 0 && sigmahat == 0 % trials can be all zero, but must have nonzero likelihood
                                muhat = 1; sigmahat = 0.2; % no idea what to put here though
                            end
                            param1_interp{s,j,k}(n,i) = muhat;
                            param2_interp{s,j,k}(n,i) = sigmahat;
                        end
                    end
                end
            end
        end
    end
end
if strcmp(likemodel,'Gamma')
    for n = 1:cell_num
        n
        for s = 1:3; % stim type
            for j = 1:3 % conflict angle
                for k = 1:2 % coherence
                    if ~isempty(resp_mean_all{n}{s,j,k})
                        for i = 1:length(XI)
                            phat = gamfit(resp_trial_all{n}{s,j,k,i});
                            param1_interp{s,j,k}(n,i) = phat(1);
                            param2_interp{s,j,k}(n,i) = phat(2);                  
                        end
                    end
                end
            end
        end
    end
end
end

% % % % Q: why do Morgan-weighted spike counts (trials) not match the real data
% % % % in terms of threshold?  or in other words, how can the real data achieve optimal
% % % % thresh despite sub-optimal neural weights?
% % % % 1) visually inspect fit quality, again
% for n = 1:cell_num
% 
%     fake_low_minus = resp_mean_all{n}{3,1,1};
%     fake_low_zero = resp_mean_all{n}{3,2,1};
%     fake_low_plus = resp_mean_all{n}{3,3,1};
%     fake_high_minus = resp_mean_all{n}{3,1,2};
%     fake_high_zero = resp_mean_all{n}{3,2,2};
%     fake_high_plus = resp_mean_all{n}{3,3,2};
%     
% %     % OR
% %     for i = 1:length(xx)
% %         fake_low_minus(i) = mean(resp_trial_all{n}{3,1,1,i});
% %         fake_low_zero(i) = mean(resp_trial_all{n}{3,2,1,i});
% %         fake_low_plus(i) = mean(resp_trial_all{n}{3,3,1,i});
% %         fake_high_minus(i) = mean(resp_trial_all{n}{3,1,2,i});
% %         fake_high_zero(i) = mean(resp_trial_all{n}{3,2,2,i});
% %         fake_high_plus(i) = mean(resp_trial_all{n}{3,3,2,i});
% %     end
% 
%     figure(13); clf; subplot(2,1,1);
%     plot(xx,interp1(X,resp_mean_all_orig{n}{3,1,1},xx),'r-o'); hold on;
%     plot(xx,interp1(X,resp_mean_all_orig{n}{3,2,1},xx),'b-o');
%     plot(xx,interp1(X,resp_mean_all_orig{n}{3,3,1},xx),'g-o');
%     subplot(2,1,2);
%     plot(xx,fake_low_minus,'r--s',xx,fake_low_zero,'b--s',xx,fake_low_plus,'g--s');
%     R1 = corrcoef(interp1(X,resp_mean_all_orig{n}{3,1,1},xx),fake_low_minus);
%     R2 = corrcoef(interp1(X,resp_mean_all_orig{n}{3,2,1},xx),fake_low_zero);
%     R3 = corrcoef(interp1(X,resp_mean_all_orig{n}{3,3,1},xx),fake_low_plus);
%     [Rall,Pall] = corrcoef([interp1(X,resp_mean_all_orig{n}{3,1,1},xx) interp1(X,resp_mean_all_orig{n}{3,2,1},xx) interp1(X,resp_mean_all_orig{n}{3,3,1},xx)],[fake_low_minus fake_low_zero fake_low_plus]);
%     r1l(n) = R1(1,2); r2l(n) = R2(1,2); r3l(n) = R3(1,2); Rall_low(n) = Rall(1,2); Pall_low(n) = Pall(1,2);
%     title([num2str([R_corr_low(n) P_corr_low(n)]) ' ;; ' num2str([r1l(n) r2l(n) r3l(n)])]);
% %     close;    
%     % sum squared error instead of R^2 for goodness of fit:
%     SSE_low(n) = sum(([interp1(X,resp_mean_all_orig{n}{3,1,1},xx) interp1(X,resp_mean_all_orig{n}{3,2,1},xx) interp1(X,resp_mean_all_orig{n}{3,3,1},xx)]-[fake_low_minus fake_low_zero fake_low_plus]).^2);
%     
%     figure(14); clf; subplot(2,1,1);
%     plot(xx,interp1(X,resp_mean_all_orig{n}{3,1,2},xx),'r-o'); hold on;
%     plot(xx,interp1(X,resp_mean_all_orig{n}{3,2,2},xx),'b-o');
%     plot(xx,interp1(X,resp_mean_all_orig{n}{3,3,2},xx),'g-o');
%     subplot(2,1,2);
%     plot(xx,fake_high_minus,'r--s',xx,fake_high_zero,'b--s',xx,fake_high_plus,'g--s');
%     R1 = corrcoef(interp1(X,resp_mean_all_orig{n}{3,1,2},xx),fake_high_minus);
%     R2 = corrcoef(interp1(X,resp_mean_all_orig{n}{3,2,2},xx),fake_high_zero);
%     R3 = corrcoef(interp1(X,resp_mean_all_orig{n}{3,3,2},xx),fake_high_plus);
%     [Rall,Pall] = corrcoef([interp1(X,resp_mean_all_orig{n}{3,1,2},xx) interp1(X,resp_mean_all_orig{n}{3,2,2},xx) interp1(X,resp_mean_all_orig{n}{3,3,2},xx)],[fake_high_minus fake_high_zero fake_high_plus]); 
%     r1h(n) = R1(1,2); r2h(n) = R2(1,2); r3h(n) = R3(1,2); Rall_high(n) = Rall(1,2); Pall_high(n) = Pall(1,2);
%     title([num2str([R_corr_high(n) P_corr_high(n)]) ' ;; ' num2str([r1h(n) r2h(n) r3h(n)])]);
% %     close;
%     % sum squared error instead of R^2 for goodness of fit:
%     SSE_high(n) = sum(([interp1(X,resp_mean_all_orig{n}{3,1,2},xx) interp1(X,resp_mean_all_orig{n}{3,2,2},xx) interp1(X,resp_mean_all_orig{n}{3,3,2},xx)]-[fake_high_minus fake_high_zero fake_high_plus]).^2);
%     
%     pause
% end

% % underestimating slope at all conflicts, or just delta=0? (A: all conflicts -- why doesn't DC fix that?!)
% for n = 1:cell_num
% %     figure(15); clf; subplot(2,1,1); j=1;
% %     plot(xx,interp1(X,resp_mean_all_orig{n}{3,j,1},xx),'r-o', xx, resp_mean_all{n}{3,j,1},'m-s');
% %     legend('data','fit');
% %     subplot(2,1,2);
% %     plot(xx,interp1(X,resp_mean_all_orig{n}{3,j,2},xx),'r-o', xx, resp_mean_all{n}{3,j,2},'m-s');
% %     
% %     figure(16); clf; subplot(2,1,1); j=2;
% %     plot(xx,interp1(X,resp_mean_all_orig{n}{3,j,1},xx),'b-o', xx, resp_mean_all{n}{3,j,1},'c-s');
% %     legend('data','fit');
% %     subplot(2,1,2);
% %     plot(xx,interp1(X,resp_mean_all_orig{n}{3,j,2},xx),'b-o', xx, resp_mean_all{n}{3,j,2},'c-s');
% % 
% %     figure(17); clf; subplot(2,1,1); j=3;
% %     plot(xx,interp1(X,resp_mean_all_orig{n}{3,j,1},xx),'g-o', xx, resp_mean_all{n}{3,j,1},'k-s');
% %     legend('data','fit');
% %     subplot(2,1,2);
% %     plot(xx,interp1(X,resp_mean_all_orig{n}{3,j,2},xx),'g-o', xx, resp_mean_all{n}{3,j,2},'k-s');
% % 
% %     % let's see the fits as they were intended:
% %     resp_mean_all_temp{n}{3,2,1} = x_low{n}(1)*resp_mean_all_orig{n}{1,2,1} + x_low{n}(2)*resp_mean_all_orig{n}{2,2,1} + x_low{n}(6);
% %     resp_mean_all_temp{n}{3,2,2} = x_high{n}(1)*resp_mean_all_orig{n}{1,2,1} + x_high{n}(2)*resp_mean_all_orig{n}{2,2,2} + x_high{n}(6);
% %     figure(18); clf; subplot(2,1,1);
% %     plot(X,resp_mean_all_orig{n}{3,2,1},'b-o', X, resp_mean_all_temp{n}{3,2,1},'c-s');
% %     legend('data','fit');
% %     subplot(2,1,2);
% %     plot(X,resp_mean_all_orig{n}{3,2,2},'b-o', X, resp_mean_all_temp{n}{3,2,2},'c-s');
% % 
% %     pause;
% 
%         
%     temp = interp1(X,resp_mean_all_orig{n}{3,2,1},xx);
%     slope_data_low(n) = abs(mean(diff(temp(2:end-1))));
%     slope_fit_low(n) = abs(mean(diff(resp_mean_all{n}{3,2,1}(2:end-1))));
%     temp = interp1(X,resp_mean_all_orig{n}{3,2,2},xx);
%     slope_data_high(n) = abs(mean(diff(temp(2:end-1))));
%     slope_fit_high(n) = abs(mean(diff(resp_mean_all{n}{3,2,2}(2:end-1))));
% end
% 
% figure; subplot(2,1,1); plot(slope_data_low, slope_fit_low,'bo',[0 0.4],[0 0.4],'k--'); % axis square;
% xlabel('slope of data');ylabel('slope of fits')
% subplot(2,1,2); plot(slope_data_high, slope_fit_high,'ro',[0 0.8],[0 0.8],'k--');  % axis square;
% xlabel('slope of data');ylabel('slope of fits')
% P=signrank(slope_data_low,slope_fit_low)
% P=signrank(slope_data_high,slope_fit_high)
% 
% mean(SSElow)  % 102.7276
% mean(SSEhigh) % 133.7985, for de-interp after fitting, versus 82.36 and 104.36 for de-interp before fitting (less of a slope problem there)
% asdf

% now it's 149.1062 and 160.6167 for de-interp after fitting,
% vs 134.8437 and 135.7033 for de-interp before fitting
% vs 82.3622 and 104.3659 for latter w/ DC (and slope diff is nonsig)

% figure;
% subplot(3,2,1); hist(r1l); title(num2str(median(r1l)));
% subplot(3,2,3); hist(r2l); title(num2str(median(r2l)));
% subplot(3,2,5); hist(r3l); title(num2str(median(r3l)));
% subplot(3,2,2); hist(r1h); title(num2str(median(r1h)));
% subplot(3,2,4); hist(r2h); title(num2str(median(r2h)));
% subplot(3,2,6); hist(r3h); title(num2str(median(r3h)));



% data = [(Rall_low').^2 Pall_low' (Rall_high').^2 Pall_high' SSE_low' SSE_high'];

% % LOAD DATA (old morgans) AND DATA_NEW (new morgans) FROM SPREADSHEET
% figure;
% subplot(2,2,1); plot(data(:,1),data_new(:,1),'bo',[0 1],[0 1],'k--');
% xlabel('R-squared, Old Morgans'); ylabel('R-squared, New Morgans'); title('Low Coherence');
% subplot(2,2,2); plot(data(:,3),data_new(:,3),'ro',[0 1],[0 1],'k--');
% xlabel('R-squared, Old Morgans'); ylabel('R-squared, New Morgans'); title('High Coherence');
% subplot(2,2,3); loglog(data(:,5),data_new(:,5),'bs',[10 10000],[10 10000],'k--');
% xlabel('SSE, Old Morgans'); ylabel('SSE, New Morgans'); title('Low Coherence'); xlim([10 10000]); ylim([10 10000]);
% subplot(2,2,4); loglog(data(:,6),data_new(:,6),'rs',[10 10000],[10 10000],'k--'); xlim([10 10000]); ylim([10 10000]);
% xlabel('SSE, Old Morgans'); ylabel('SSE, New Morgans'); title('High Coherence');


% for n = 1:cell_num
%     figure(906); title(['#' num2str(n) ': ' files{n}]); hold on;
%     plot(unique_heading, resp_mean_all{n}{1,2,1}, 'k', 'LineWidth', 3);
%     plot(unique_heading, resp_mean_all{n}{2,2,1}, 'm', 'LineWidth', 3);
%     plot(unique_heading, resp_mean_all{n}{2,2,2}, 'r', 'LineWidth', 3);
%     plot(unique_heading, resp_mean_all{n}{3,2,1}, 'c', 'LineWidth', 3);
%     plot(unique_heading, resp_mean_all{n}{3,2,2}, 'b', 'LineWidth', 3);
% %     xlabel(['slopeRatio-low = ' num2str(slopeRatio_low(n)) ' meanRatio-low = ' num2str(meanRatio_low(n))]);
% %     title(['Cell #' num2str(n) ' WoptLow = ' num2str(x_low(n,1)) ' ; WoptHigh = ' num2str(x_high(n,1))]);
%     title(['Cell #' num2str(n)]);
%     pause; clf;
% end



%------------------------------------------------
for n = 1:num_stim_reps
    tic
    n
    for t = 1:length(condition_list)
        trial = (n-1)*length(condition_list)+t;
        stim(trial) = stimtype(t);
        coh(trial) = coherence(t);
        hdg(trial) = heading(t);
        delt(trial) = delta(t);
        
        if unique_stim_type == 1
            deltnum = 2;
        else
            deltnum = find(unique(delta)==delt(trial));
        end
        cohnum = find(unique(coherence)==coh(trial));
        hdgnum = find(unique(heading)==hdg(trial));

% %         since debug mode is not working:
        s=stim(trial);j=deltnum;k=cohnum;i=hdgnum;

        if s<3
            likemodel_to_pass_in = noisemodel;
        else
            likemodel_to_pass_in = likemodel;
        end



%         for asdf = 1:4

%         %%%%%%%
%         for abcd = 1:20
%         s=1;j=2;k=1;i=26;likemodel_to_pass_in = noisemodel;
%         %%%%%%%

        [r,L] = compute_likelihood_backup(unique_heading',resp_trial_all,resp_mean_all,s,j,k,i,use_logL,param1_interp{s,2,k},param2_interp{s,2,k},likemodel_to_pass_in);
        L = L/sum(L); % normalize
        
%         Post = L; % no prior for now        
%         figure(1); clf; plot(x,Post); 
%         title([num2str( sum(Post(1:find(x==0)-1))) '   ' num2str(sum(Post(find(x==0)+1:end)))]);
%         tempchoice = (sum(Post(find(x==0)+1:end)) > sum(Post(1:find(x==0)-1))) * 2 - 1;
%         xlabel(num2str(tempchoice));
%         pause
%         %%%%%%%
%         end % abcd
%         %%%%%%%               

%         AAAA(:,asdf)=L';
%         end
%         figure;plot(xi,AAAA); legend('1','2','3','4');

        Post = L; % no prior for now

        if sum(Post)==0 || isnan(sum(Post))
            disp('Likelihood is all zero or NaN -- WTF?');
%             figure(80); clf; plot(L); pause;
            choice(trial) = sign(randn);
        else
            if Jan % alternate decision rule: compare area under posterior <0 vs. >0
                choice(trial) = (sum(Post(find(x==0)+1:end)) > sum(Post(1:find(x==0)-1))) * 2 - 1;
            else % simple MAP rule
                findmap = find(Post==max(Post));
                MAP = x(findmap(1));
                choice(trial) = sign(MAP);
            end
            if choice(trial)==0
                choice(trial) = sign(randn);
            end
        end
        
    end
    toc

end % end num_stim_reps

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
                elseif sum(trials_select) ~= num_stim_reps
                    disp('ERROR: num_stim_reps');
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
        

%------------------------------------------------------------------------
%Resample and re-fit data N times to build distributions of bias and threshold.
%These will be drawn from to compute confidence intervals for the weights.
if run_resample
N_boots = 1000; prog = 0; endprog = 1 + length(unique_coherence) + length(unique_delta)*length(unique_coherence);
Thresh_psy_resamp = []; Bias_psy_resamp = [];
for s = 1:length(unique_stim_type)
    for j = 1:length(unique_delta)
        for k = 1:length(unique_coherence)
            if fit_data_psycho_cum{s,j,k}(1,3) == 0 % identifies and skips invalid condition combinations (e.g., vestibular only with a nonzero conflict angle)
                Thresh_psy_resamp{s,j,k} = NaN;
                Bias_psy_resamp{s,j,k} = NaN;
            else
                prog = prog+1; tic;
                clear fit_data_resamp; fit_data_resamp(:,1) = fit_data_psycho_cum{s,j,k}(:,1);
                for i = 1:length(fit_data_psycho_cum{s,j,k}(:,1))
                    data_to_resample{i} = zeros(1,fit_data_psycho_cum{s,j,k}(i,3));
                    num_rightward_choices(i) = round(fit_data_psycho_cum{s,j,k}(i,2).*fit_data_psycho_cum{s,j,k}(i,3));
                    data_to_resample{i}(1:num_rightward_choices(i)) = 1;
                    data_to_resample{i} = data_to_resample{i}(randperm(length(data_to_resample{i})));
                end
                disp(['Progress = ' num2str(prog) ' out of ' num2str(endprog)]);
                for n = 1:N_boots
%                     n
                    for i = 1:length(fit_data_psycho_cum{s,j,k}(:,1))
                        clear randtemp randsamp
                        for h = 1:length(data_to_resample{i});
                            randtemp = randperm(length(data_to_resample{i}));
                            randsamp(h) = randtemp(1);
                        end
                        fit_data_resamp(i,2) = sum(data_to_resample{i}(randsamp))/length(data_to_resample{i});
                        fit_data_resamp(i,3) = length(data_to_resample{i});
                    end
                    if wichmann
                        wichman_resamp = pfit(fit_data_resamp,'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'sens',0,'compute_stats','false','verbose','false');  
                        Thresh_psy_resamp{s,j,k}(n) = wichman_resamp.params.est(2);
                        Bias_psy_resamp{s,j,k}(n) = wichman_resamp.params.est(1);
% %                   % to save time, use a different fit method:      %NO, this doesn't work here -- threshold is overestimated% %                     [bb,tt] = cum_gaussfit_max1(fit_data_resamp);
% %                     Bias_psy_resamp{s,j,k}(n) = bb;
% %                     Thresh_psy_resamp{s,j,k}(n) = tt;
                    else
                        [bb,tt] = cum_gaussfit_max1(fit_data_resamp);
                        Thresh_psy_resamp{s,j,k}(n) = tt;
                        Bias_psy_resamp{s,j,k}(n) = bb;
                    end
                end
                toc;
            end
        end
    end
end

% find confidence intervals on the weights/thresh by drawing
%(with replacement) from bootstrapped bias/thresh distributions
for k = 1:length(unique_coherence)
    for m = 1:N_boots
%                 m
        randtemp = randperm(N_boots);
        randdraw = randtemp(1);        % NOTE: must add "1 - " to weights for Iolaus_2I
        Wves_actual_minus_boot{k}(m) = ((Bias_psy_resamp{3,1,k}(randdraw) - Bias_psy_resamp{3,2,k}(randdraw)) - (-unique_delta(1)/2)) / unique_delta(1);
        Wves_actual_plus_boot{k}(m) = ((Bias_psy_resamp{3,3,k}(randdraw) - Bias_psy_resamp{3,2,k}(randdraw)) - (-unique_delta(3)/2)) / unique_delta(3);
        Wves_actual_boot{k}(m) = (Wves_actual_minus_boot{k}(m) + Wves_actual_plus_boot{k}(m)) / 2;
        Wves_predicted_boot{k}(m) = Thresh_psy_resamp{2,2,k}(randdraw)^2/(Thresh_psy_resamp{1,2,1}(randdraw)^2+Thresh_psy_resamp{2,2,k}(randdraw)^2);
        thresh_predicted_boot{k}(m) = sqrt((Thresh_psy_resamp{2,2,k}(randdraw)^2*Thresh_psy_resamp{1,2,1}(randdraw)^2)/(Thresh_psy_resamp{2,2,k}(randdraw)^2+Thresh_psy_resamp{1,2,1}(randdraw)^2));
        thresh_actual_d0_boot{k}(m) = Thresh_psy_resamp{3,2,k}(randdraw);
        thresh_actual_all_boot{k}(m) = (Thresh_psy_resamp{3,1,k}(randdraw)+Thresh_psy_resamp{3,2,k}(randdraw)+Thresh_psy_resamp{3,3,k}(randdraw)) / 3;
        thresh_actual_pm_boot{k}(m) = (Thresh_psy_resamp{3,1,k}(randdraw)+Thresh_psy_resamp{3,3,k}(randdraw)) / 2;
        bias_delta0_boot{k}(m) = Bias_psy_resamp{3,2,k}(randdraw);
        thresh_ves_boot{k}(m) = Thresh_psy_resamp{1,2,1}(randdraw);
        thresh_vis_boot{k}(m) = Thresh_psy_resamp{2,2,k}(randdraw);
    end
    
%     % OPTIONAL: clip actual weights at 0 and 1
%     Wves_actual_minus_boot{k}(Wves_actual_minus_boot{k}<0) = 0;
%     Wves_actual_minus_boot{k}(Wves_actual_minus_boot{k}>1) = 1;
%     Wves_actual_plus_boot{k}(Wves_actual_plus_boot{k}<0) = 0;
%     Wves_actual_plus_boot{k}(Wves_actual_plus_boot{k}>1) = 1;
%     Wves_actual_boot{k}(Wves_actual_boot{k}<0) = 0;
%     Wves_actual_boot{k}(Wves_actual_boot{k}>1) = 1;

    boot_sorted = sort(Wves_actual_minus_boot{k});
    Wves_actual_minus_lower(k) = boot_sorted(round(.025*N_boots)+1);
    Wves_actual_minus_upper(k) = boot_sorted(round(.975*N_boots)-1);

    boot_sorted = sort(Wves_actual_plus_boot{k});
    Wves_actual_plus_lower(k) = boot_sorted(round(.025*N_boots)+1);
    Wves_actual_plus_upper(k) = boot_sorted(round(.975*N_boots)-1);

    boot_sorted = sort(Wves_actual_boot{k});
    Wves_actual_lower(k) = boot_sorted(round(.025*N_boots)+1);
    Wves_actual_upper(k) = boot_sorted(round(.975*N_boots)-1);

    boot_sorted = sort(Wves_predicted_boot{k});
    Wves_predicted_lower(k) = boot_sorted(round(.025*N_boots)+1);
    Wves_predicted_upper(k) = boot_sorted(round(.975*N_boots)-1);

    boot_sorted = sort(thresh_predicted_boot{k});
    thresh_predicted_lower(k) = boot_sorted(round(.025*N_boots)+1);
    thresh_predicted_upper(k) = boot_sorted(round(.975*N_boots)-1);

    boot_sorted = sort(thresh_actual_d0_boot{k});
    thresh_actual_d0_lower(k) = boot_sorted(round(.025*N_boots)+1);
    thresh_actual_d0_upper(k) = boot_sorted(round(.975*N_boots)-1);

    boot_sorted = sort(thresh_actual_all_boot{k});
    thresh_actual_all_lower(k) = boot_sorted(round(.025*N_boots)+1);
    thresh_actual_all_upper(k) = boot_sorted(round(.975*N_boots)-1);

    boot_sorted = sort(thresh_actual_pm_boot{k});
    thresh_actual_pm_lower(k) = boot_sorted(round(.025*N_boots)+1);
    thresh_actual_pm_upper(k) = boot_sorted(round(.975*N_boots)-1);

    boot_sorted = sort(bias_delta0_boot{k});
    bias_delta0_lower(k) = boot_sorted(round(.025*N_boots)+1);
    bias_delta0_upper(k) = boot_sorted(round(.975*N_boots)-1);
    
    boot_sorted = sort(thresh_ves_boot{k});
    thresh_ves_lower(k) = boot_sorted(round(.025*N_boots)+1);
    thresh_ves_upper(k) = boot_sorted(round(.975*N_boots)-1);
    
    boot_sorted = sort(thresh_vis_boot{k});
    thresh_vis_lower(k) = boot_sorted(round(.025*N_boots)+1);
    thresh_vis_upper(k) = boot_sorted(round(.975*N_boots)-1);
end

end



if plotflag
%--------------------------------------------------------------------------
% plot psychometric function here
clear F
H{1,1} = 'bo'; H{2,1} = 'b^'; H{3,1} = 'bs'; H{4,1} = 'b*'; F{1} = 'b-';
H{1,2} = 'go'; H{2,2} = 'g^'; H{3,2} = 'gs'; H{4,2} = 'g*'; F{2} = 'g-';
H{1,3} = 'ro'; H{2,3} = 'r^'; H{3,3} = 'rs'; H{4,3} = 'r*'; F{3} = 'r-';
H{1,4} = 'co'; H{2,4} = 'c^'; H{3,4} = 'cs'; H{4,4} = 'c*'; F{4} = 'c-';
H{1,5} = 'mo'; H{2,5} = 'm^'; H{3,5} = 'ms'; H{4,5} = 'm*'; F{5} = 'm-';
H{1,6} = 'yo'; H{2,6} = 'y^'; H{3,6} = 'ys'; H{4,6} = 'y*'; F{6} = 'y-';
H{1,7} = 'ko'; H{2,7} = 'k^'; H{3,7} = 'ks'; H{4,7} = 'k*'; F{7} = 'k-';
t=1;
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
text(45,10, ['repeats = ' num2str(num_stim_reps)]);
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
    text(45,10, ['repeats = ' num2str(num_stim_reps)]);
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

Thresh_ves(2) = Thresh_ves(1);
Thresh_vis = Thresh_vis';
Thresh_ves = Thresh_ves';

if run_resample
    figure; hold on;
    errorbar([1 2], W_ves_predicted, W_ves_predicted-Wves_predicted_lower', Wves_predicted_upper'-W_ves_predicted, 'bo-')
    errorbar([1 2], W_ves_actual, W_ves_actual-Wves_actual_lower', Wves_actual_upper'-W_ves_actual, 'ro-')
    xlim([0.5 2.5]); set(gca,'xtick', [1 2], 'xticklabel',[16 60]);
    title(['Weights, wich=' num2str(wichmann)]); legend('Predicted','Actual');
    
    figure; hold on;
    errorbar([1 2], Thresh_predicted, Thresh_predicted-thresh_predicted_lower', thresh_predicted_upper'-Thresh_predicted, 'co-')
    errorbar([1 2], Thresh_actual, Thresh_actual-thresh_actual_d0_lower', thresh_actual_d0_upper'-Thresh_actual, 'bo-');
    errorbar([1 2], Thresh_vis, Thresh_vis-thresh_vis_lower', thresh_vis_upper'-Thresh_vis, 'rx-');
    errorbar([1 2], Thresh_ves, Thresh_ves-thresh_ves_lower', thresh_ves_upper'-Thresh_ves, 'kx-');
    xlim([0.5 2.5]); set(gca,'xtick', [1 2], 'xticklabel',[16 60]);
    title(['Thresholds, wich=' num2str(wichmann)]); legend('Predicted','Actual','Visual','Vestib');
    AAorigin1 = [W_ves_predicted W_ves_predicted-Wves_predicted_lower' Wves_predicted_upper'-W_ves_predicted W_ves_actual W_ves_actual-Wves_actual_lower' Wves_actual_upper'-W_ves_actual];
    AAorigin2 = [Thresh_predicted Thresh_predicted-thresh_predicted_lower' thresh_predicted_upper'-Thresh_predicted Thresh_actual Thresh_actual-thresh_actual_d0_lower' thresh_actual_d0_upper'-Thresh_actual Thresh_ves Thresh_ves-thresh_ves_lower' thresh_ves_upper'-Thresh_ves Thresh_vis Thresh_vis-thresh_vis_lower' thresh_vis_upper'-Thresh_vis];
    AAorigin = [AAorigin1 AAorigin2];
else
    figure; hold on;
    plot([1 2], W_ves_predicted, 'bo-', [1 2], W_ves_actual, 'ro-');
    xlim([0.5 2.5]); set(gca,'xtick', [1 2], 'xticklabel',[16 60]);
    title(['Weights, wich=' num2str(wichmann)]); legend('Predicted','Actual');
    
    figure; hold on;
	plot([1 2], Thresh_predicted, 'co-', [1 2], Thresh_actual, 'bo-');
    plot([1 2], Thresh_vis, 'rx-', [1 2], Thresh_ves, 'kx-');
    xlim([0.5 2.5]); set(gca,'xtick', [1 2], 'xticklabel',[16 60]);
    title(['Thresholds, wich=' num2str(wichmann)]); legend('Predicted','Actual','Visual','Vestib');    
    AAorigin1 = [W_ves_predicted [0;0] W_ves_actual [0;0]];
    AAorigin2 = [Thresh_predicted [0;0] Thresh_actual [0;0] Thresh_ves [0;0] Thresh_vis [0;0]];
    AAorigin = [AAorigin1 AAorigin2];
end

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
save(['likelihood_results_' batch_name '_' noisemodel '_' datestamp '_sponsub' num2str(subtract_baseline) '.mat'])
% save(['C:\Documents and Settings\Chris\My Documents\MATLAB\' timefolder(1:end-1) '_all_keepreals_' noisemodel '.mat'])
end

clock
