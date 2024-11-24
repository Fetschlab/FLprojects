% population_likelihood.m 
% Trial-by-trial simulation of the behavioral experiment, computing
% likelihood functions by drawing spike counts from real neuronal data.

% CF rebooted in 2021 based on old code from 2011 paper; Fun!

% CF 07/2023: added previous trial conditionalizaton

clear
close all

cd /Users/chris/Documents/MATLAB/Projects/LikelihoodDecodingNN
% ^ change to your folder

load('Fetsch_et_al_NatNeuro_2011.mat')
load('condition_list_neuro.mat')

unique_stim_type = unique(data.modality);
unique_coherence = unique(data.coherence);
unique_heading = unique(data.heading);
unique_delta = unique(data.delta);

num_stim_reps = 50; % number of simulated repetitions of each condition

stimtype = condition_list(:,1); % what we used to call modality
coherence = condition_list(:,2);
heading = condition_list(:,3);
    heading(heading==1.225) = 1.23; heading(heading==-1.225) = -1.23; % fix a rounding error
delta = condition_list(:,4);

total_tr = length(stimtype) * num_stim_reps;
stim = nan(total_tr,1);
coh = nan(total_tr,1);
hdg = nan(total_tr,1);
delt = nan(total_tr,1);
choice = nan(total_tr,1);

noisemodel = 'Poisson';    % refers to generation of interpolated single-cue trials
% noisemodel = 'Gaussian'
% noisemodel = 'Gamma'

likemodel = 'Poisson';  % refers to 'noise' distribution assumed by likelihood calculation
% likemodel = 'Gaussian'
% likemodel = 'Gamma'
% likemodel = 'histogram'  % no need to use histogram unless generating fake combs as weighted sums of single trials

% these are needed for the interpolation in the case of Gaussian or Gamma,
% but don't remember how they were set (just stick w Poisson for simplicity)
param1_interp = [];
param2_interp = [];

use_logL = 1;
Jan = 1; % Drugowitsch must have told me to do area-under-posterior! hah
 
plotflag = 1;
repeat_sim = 1; % not sure what this was used for


%% create resp mats in old format

ufile = unique(data.filename);
resp_mean_all = cell(length(ufile),1);
resp_trial_all = cell(length(ufile),1);
resp_trial_subset = cell(length(ufile),1);
nreps = nan(length(ufile),length(unique_stim_type),length(unique_delta),length(unique_coherence),length(unique_heading));


% To condition on previous trial type using logical vectors.
% A key question is do we only need to do this for comb + nonzero delta
% trials, since that's the only time the weights are relevant
% (weights may change on zero delta trials too, but this will have no
% effect on the firing rates--or will it???)
cndn = false(size(data.filename)); % cndn = "conditionalization"
for d = 3:length(data.filename)
    % here are some previous trial conditions ("PTC"):
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% going back two trials %%% 
    % previous two trials were vestib (forget it, it's like 1%)
%     PTC = data.modality(d-1)==1 && data.modality(d-2)==1;
    
    % previous two trials were high coh (and not vestib) [works! still noisy]
%     PTC = (data.coherence(d-1)==60 && data.coherence(d-2)==60) && (data.modality(d-1)~=1 && data.modality(d-2)~=1);

    % previous two trials were low coh (and not vestib) [works! still noisy]
    PTC = (data.coherence(d-1)==16 && data.coherence(d-2)==16) && (data.modality(d-1)~=1 && data.modality(d-2)~=1);

    % for each of the above, make sure last two trials were in the same block
    PTC = PTC && strcmp(data.filename{d},data.filename{d-1}) && strcmp(data.filename{d},data.filename{d-2}); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% going back one trial %%% 
    % previous trial was vestib (about 11%, actually works but SE of weights will be too large to be sure of anything)
%     PTC = data.modality(d-1)==1;

    % previous trial was high coh (and not vestib):
%     PTC = data.coherence(d-1)==60 && data.modality(d-1)~=1 

    % previous trial was low coh (and not vestib):
%     PTC = data.coherence(d-1)==16 && data.modality(d-1)~=1 && strcmp(data.filename{d},data.filename{d-1});

    % for each of the above, make sure last trial was in the same block
    % PTC = PTC && strcmp(data.filename{d},data.filename{d-1}); % for each, make sure last trial was in the same block
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % the logic is:
    % among all comb nonzero delta trials, if previous N trials are the same block AND meet the PTC, then keep.
    % also keep all non-comb trials and zero-delta comb trials, both of which are covered by delta==0
    if (data.modality(d)==3 && data.delta(d)~=0 && PTC) || data.delta(d)==0
        cndn(d) = true;
    end
end
sum(cndn & data.modality==3 & data.delta~=0)/sum(data.modality==3 & data.delta~=0) % what pct of comb-conflict trials are included?


%%

tic
for n = 1:length(ufile)
    resp_mean_all{n} = cell(length(unique_stim_type),length(unique_delta),length(unique_coherence));
    resp_trial_all{n} = cell(length(unique_stim_type),length(unique_delta),length(unique_coherence),length(unique_heading));
    resp_trial_subset{n} = cell(length(unique_stim_type),length(unique_delta),length(unique_coherence),length(unique_heading));
    
    for s = 1:length(unique_stim_type)
        for j = 1:length(unique_delta)
            for k = 1:length(unique_coherence)
                resp_mean_all{n}{s,j,k} = nan(1,length(unique_heading));
                for i = 1:length(unique_heading)
                    
                    I = strcmp(data.filename,ufile(n)) & data.modality==unique_stim_type(s) & data.delta==unique_delta(j) & data.coherence==unique_coherence(k) & data.heading==unique_heading(i);
                    if sum(I)==0
                        resp_mean_all{n}{s,j,k}(i) = NaN;
                        resp_trial_all{n}{s,j,k,i} = NaN;
                        resp_trial_subset{n}{s,j,k,i} = NaN;
                    else
                        resp_mean_all{n}{s,j,k}(i) = mean(data.spikeRates(I));
                        resp_trial_all{n}{s,j,k,i} = data.spikeRates(I);
                        J = I & cndn; % here's where we bring in the PTC
                        resp_trial_subset{n}{s,j,k,i} = data.spikeRates(J);                        
                    end
                    
                end
            end
        end
    end
end
toc




save temp.mat % because the above is slow



%%

clear
load temp.mat % because the above is slow



% %% plot (some) tuning curves - optional
% 
% for n = 1:10 % length(ufile)
%     figure(1); clf; set(gcf,'Position', [200 200 560 780]); clf;
%     
%     subplot(2,1,1); title([num2str(n) ' - low coh']); % title doesn't show up for some reason
%     plot(unique_heading,resp_mean_all{n}{1,2,1},'ko-'); hold on;
%     plot(unique_heading,resp_mean_all{n}{2,2,1},'ro-'); 
%     plot(unique_heading,resp_mean_all{n}{3,2,1},'go-'); 
% 
%     subplot(2,1,2); title([num2str(n) ' - high coh']);
%     plot(unique_heading,resp_mean_all{n}{1,2,1},'ks-'); hold on;
%     plot(unique_heading,resp_mean_all{n}{2,2,2},'rs-'); 
%     plot(unique_heading,resp_mean_all{n}{3,2,2},'gs-');
%     
%     % an excercise: make y-axis ranges the same, and fix titles, add error
%     % bars if wanted
% 
%     pause;
% end

%% Likelihood-based decoding

% for all trials:
% resp_trial_subset = resp_trial_all;
% otherwise subset will be used.

x = -10:0.1:10;

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
     
        if stim(trial)<3 % ??? I guess this was related to interpolation? weird
            likemodel_to_pass_in = noisemodel;
        else
            likemodel_to_pass_in = likemodel;
        end

        [r,L] = compute_likelihood_backup(unique_heading',resp_trial_subset,resp_trial_all,resp_mean_all,stim(trial),deltnum,cohnum,hdgnum,use_logL,[],[],likemodel_to_pass_in);
        Post = L/sum(L); % no prior for now; posterior is normalized likelihood

        if sum(Post)==0 || isnan(sum(Post))
            disp('Likelihood is all zero or NaN -- WTF?');
            figure(80); clf; plot(L); pause;
            pause
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

end % end likelihood decoding


%% parse data for plots & weight/threshold calculation

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
                elseif sum(trials_select) ~= num_stim_reps
                    disp('ERROR: num_stim_reps mismatch');
                    return;
                else
                    right_trials = (trials_select & (choice == RIGHT));
                    right_pct{s,j,k}(i) = sum(right_trials) / sum(trials_select);
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
                [bbb,ttt] = cum_gaussfit_max1(fit_data_psycho_cum{s,j,k,repeat_sim});
                Bias_psy{s,j,k} = bbb;
                Thresh_psy{s,j,k} = ttt;
                psy_perf{s,j,k} = [bbb,ttt];
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
        


%% plot

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
text(15,10, ['stimtype = 1+2']);
text(45,10, ['repeats = ' num2str(num_stim_reps)]);
text(0,8.3, 'stimtype');
text(10,8.3, 'coherence');
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
    text(15,10, ['stimtype = ' num2str(unique_stim_type(s))]);
    text(45,10, ['repeats = ' num2str(num_stim_reps)]);
    text(0,8.3, 'delta');
    text(8,8.3, 'coherence');
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



