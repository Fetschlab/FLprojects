% Trial-by-trial simulation of the behavioral experiment, computing
% likelihood functions by drawing spike counts from real neuronal data.
% The specific dataset is defined within poisson_likelihood_realdata.m.

% close all
clear all

% global Tvislow Tvishigh Tves Tcomblow Tcombhigh repeat_sim
global W_ves_actual W_ves_predicted Thresh_actual Thresh_predicted repeat_sim
repeat_sim = 1;

while repeat_sim < 2

clear
tic
load('C:\Program Files\MATLAB\R2007a\work\condition_list_neuro.mat');

wichmann = 1;
plotflag = 1; % for psychometrics
prior = 0;
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
num_reps = 20;

%------------------------------------------------
% prior #3: choose a continuous distribution of which our heading angles 
% could be considered an approximation when drawn from using log-spaced
% bins (i.e., a distribution whose area-under-curve increases linearly with
% each doubling (or 1/0.35 -ing) on the X-axis).  But what dist. does this?

% lambda = 0.16;
% % PrDist = 1/log(1/0.35) .* 1./(x(x>0));
% PrDist = lambda*exp(-lambda*x(x>0));
% Pr = [fliplr(PrDist) PrDist(1) PrDist];
% % Pr = decimate(conv(Pr,normpdf(x,0,sig)/sum(normpdf(x,0,sig))), 2);
% Pr = Pr / sum(Pr);
% Pr(Pr<0)=0;

% or,
sig = 5;
Pr = normpdf(x,0,5);
Pr = Pr / sum(Pr);

%------------------------------------------------
total_tr = length(stimtype) * num_reps;
stim = zeros(1,total_tr); coh = zeros(1,total_tr);
hdg = zeros(1,total_tr); delt= zeros(1,total_tr);
choice = zeros(1,total_tr);
for n = 1:num_reps
    n
    cond = randperm(length(stimtype));
    for t = 1:length(cond)
        trial = (n-1)*length(cond)+t;
        stim(trial) = stimtype(cond(t));
        coh(trial) = coherence(cond(t));
        hdg(trial) = heading(cond(t));
        delt(trial) = delta(cond(t));
        deltnum = find(unique(delta)==delt(trial));
        cohnum = find(unique(coherence)==coh(trial));
        hdgnum = find(unique(heading)==hdg(trial));

        L = poisson_likelihood_realdata(stim(trial),deltnum,cohnum,hdgnum);
%         sum(L)
        L = L/sum(L); % normalize
        if prior
            Post = L.*Pr; % multiply by prior
        else
            Post = L;
        end
        if sum(Post)>0
            Post = Post/sum(Post); % normalize
            findmap = find(Post==max(Post));
            MAP = x(findmap(1));
            choice(trial) = sign(MAP);
        else
            disp('Likelihood is all zero -- WTF?');
            choice(trial) = sign(randn);
        end
    end
end

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
                fit_data_psycho_cum{s,j,k}(i,1) = unique_heading(i);
                fit_data_psycho_cum{s,j,k}(i,2) = right_pct{s,j,k}(i);
                fit_data_psycho_cum{s,j,k}(i,3) = sum(trials_select);
            end
        end
	end
end


%%%%%% use either pfit or cum_gaussfit_max1 to estimate threshold and bias
for s = 1:length(unique_stim_type)
    for j = 1:length(unique_delta)
        for k = 1:length(unique_coherence)
            if fit_data_psycho_cum{s,j,k}(1,3) == 0 % identifies and skips invalid condition combinations (e.g., vestibular only with a nonzero conflict angle)
                Thresh_psy{s,j,k} = NaN;
                Bias_psy{s,j,k} = NaN;
                psy_perf{s,j,k} = [NaN , NaN];
            else
                if wichmann
                    wichman_psy = pfit(fit_data_psycho_cum{s,j,k},'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'sens',0,'compute_stats','false','verbose','false');
                    Bias_psy{s,j,k} = wichman_psy.params.est(1);
                    Thresh_psy{s,j,k} = wichman_psy.params.est(2);
                    psy_perf{s,j,k} = [wichman_psy.params.est(1),wichman_psy.params.est(2)];
                else
                    [bbb,ttt] = cum_gaussfit_max1(fit_data_psycho_cum{s,j,k});
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
if length(unique_delta) > 1
    if length(unique_stim_type) > 1
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
            thresh_actual_d0(k) = Thresh_psy{1,1,k};
            thresh_actual_all(k) = NaN;
            thresh_actual_pm(k) = NaN;
            bias_delta0(k) = Bias_psy{1,1,k};
        end
    end
else
    if length(unique_stim_type) > 1
        for k = 1:length(unique_coherence)
            Wves_actual_minus(k) = NaN;
            Wves_actual_plus(k) = NaN;
            Wves_actual(k) = NaN;
            Wves_predicted(k) = NaN;
            thresh_predicted(k) = sqrt((Thresh_psy{2,1,k}^2*Thresh_psy{1,1,1}^2)/(Thresh_psy{2,1,k}^2+Thresh_psy{1,1,1}^2));
            thresh_actual_d0(k) = Thresh_psy{3,1,k};
            thresh_actual_all(k) = NaN;
            thresh_actual_pm(k) = NaN;
            bias_delta0(k) = Bias_psy{3,1,k};
        end
    else
        Wves_actual_minus(k) = NaN;
        Wves_actual_plus(k) = NaN;
        Wves_actual(k) = NaN;
        Wves_predicted(k) = NaN;
        thresh_predicted(k) = NaN;
        thresh_actual_d0(k) = Thresh_psy{1,1,k};
        thresh_actual_all(k) = NaN;
        thresh_actual_pm(k) = NaN;
        bias_delta0(k) = Bias_psy{1,1,k};
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
            plot(unique_heading, fit_data_psycho_cum{s,j,k}(:,2), H{k,s}, x, cum_gaussfit(psy_perf{s,j,k}, x), F{s});
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
        plot(unique_heading, fit_data_psycho_cum{s,j,k}(:,2), H{k,j}, x, cum_gaussfit(psy_perf{s,j,k}, x),  F{j} );
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
        text(37,8-j,num2str(num2str(mean(mean([1-fit_data_psycho_cum{s,j,k}(1:3,2) fit_data_psycho_cum{s,j,k}(5:7,2)])))));
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
    
end % plotflag

% global Tvislow Tvishigh Tves Tcomblow Tcombhigh repeat_sim
% Tvislow(repeat_sim) = Thresh_psy{2,2,1};
% Tvishigh(repeat_sim) = Thresh_psy{2,2,2};
% Tves(repeat_sim) = Thresh_psy{1,2,1};
% Tcomblow(repeat_sim) = Thresh_psy{3,2,1};
% Tcombhigh(repeat_sim) = Thresh_psy{3,2,2};

global W_ves_actual W_ves_predicted Thresh_actual Thresh_predicted repeat_sim Tvis Tves e_cell_psy_raw e_cell_psy_fit
W_ves_actual(:,repeat_sim) = Wves_actual;
W_ves_predicted(:,repeat_sim) = Wves_predicted;
Thresh_actual(:,repeat_sim) = thresh_actual_all;
Thresh_predicted(:,repeat_sim) = thresh_predicted;
Tvis(repeat_sim, :) = horzcat(Thresh_psy{2,2,:});
Tves(repeat_sim, :) = horzcat(Thresh_psy{1,2,:});

e_cell_psy_raw{repeat_sim} = [unique_heading right_pct{2,2,1}' right_pct{2,2,2}' right_pct{1,2,1}' right_pct{3,1,1}' right_pct{3,2,1}' right_pct{3,3,1}' right_pct{3,1,2}' right_pct{3,2,2}' right_pct{3,3,2}'];
xi = -10:0.1:10;
e_cell_psy_fit{repeat_sim} = [xi' cum_gaussfit(psy_perf{2,2,1}, xi)' cum_gaussfit(psy_perf{2,2,2}, xi)' cum_gaussfit(psy_perf{1,2,1}, xi)' cum_gaussfit(psy_perf{3,1,1}, xi)' cum_gaussfit(psy_perf{3,2,1}, xi)' cum_gaussfit(psy_perf{3,3,1}, xi)' cum_gaussfit(psy_perf{3,1,2}, xi)' cum_gaussfit(psy_perf{3,2,2}, xi)' cum_gaussfit(psy_perf{3,3,2}, xi)'];

repeat_sim = repeat_sim + 1
toc
end % repeat_sim

if repeat_sim > 2
    figure; plot(W_ves_predicted(1,:),W_ves_actual(1,:),'x',W_ves_predicted(2,:),W_ves_actual(2,:),'o',[0 1],[0 1],'k--');
    axis square; title('Weights');
    figure; plot(Thresh_predicted(1,:),Thresh_actual(1,:),'x',Thresh_predicted(2,:),Thresh_actual(2,:),'o',[0 2],[0 2],'k--');
    axis square; title('Thresholds');

    figure; errorbar([1 2], mean(W_ves_predicted'), std(W_ves_predicted')/sqrt(repeat_sim-1), 'bo-')
    hold on; errorbar([1 2], mean(W_ves_actual'), std(W_ves_actual')/sqrt(repeat_sim-1), 'ro-');
    xlim([0.5 2.5]); set(gca,'xtick', [1 2], 'xticklabel',[16 60]);
    title('Weights'); legend('Predicted','Actual');
    
    figure; errorbar([1 2], mean(Thresh_predicted'), std(Thresh_predicted')/sqrt(repeat_sim-1), 'bo-')
    hold on; errorbar([1 2], mean(Thresh_actual'), std(Thresh_actual')/sqrt(repeat_sim-1), 'ro-');
    xlim([0.5 2.5]); set(gca,'xtick', [1 2], 'xticklabel',[16 60]);
    title('Thresholds'); legend('Predicted','Actual');
end

AAorigin1 = [mean(W_ves_predicted')' std(W_ves_predicted')'/sqrt(repeat_sim-1) mean(W_ves_actual')' std(W_ves_actual')'/sqrt(repeat_sim-1)];
% AAorigin1 = [mean(W_ves_actual')' std(W_ves_actual')'/sqrt(repeat_sim-1)];
AAorigin2 = [mean(Thresh_predicted')' std(Thresh_predicted')'/sqrt(repeat_sim-1) mean(Thresh_actual')' std(Thresh_actual')'/sqrt(repeat_sim-1)];
% AAorigin2 = [mean(Thresh_actual')' std(Thresh_actual')'/sqrt(repeat_sim-1)];

repeat = repeat_sim-1;
e_cell_sum_fit = zeros(size(e_cell_psy_fit{repeat}));
e_cell_sum_raw = zeros(size(e_cell_psy_raw{repeat}));
for h = 1:repeat
	e_cell_sum_fit = e_cell_sum_fit + e_cell_psy_fit{h};
    e_cell_sum_raw = e_cell_sum_raw + e_cell_psy_raw{h};
end
AAAe_cell_mean_fit = e_cell_sum_fit / repeat;
AAAe_cell_mean_raw = e_cell_sum_raw / repeat;

% for scatter plot
% p = 1;
% for q = 1:2:(repeat_sim-1)*2
%     Aorigin(q,1) = W_ves_predicted(1,p);
%     Aorigin(q,2) = W_ves_actual(1,p);
%     Aorigin(q,3) = unique_coherence(1);
%     Aorigin(q,4) = Thresh_predicted(1,p);
%     Aorigin(q,5) = Thresh_actual(1,p);
%     Aorigin(q,6) = unique_coherence(1);
%     
%     Aorigin(q+1,1) = W_ves_predicted(2,p);
%     Aorigin(q+1,2) = W_ves_actual(2,p);
%     Aorigin(q+1,3) = unique_coherence(2);
%     Aorigin(q+1,4) = Thresh_predicted(2,p);
%     Aorigin(q+1,5) = Thresh_actual(2,p);
%     Aorigin(q+1,6) = unique_coherence(2);
%     p = p+1
% end

% % for line plot
% AAAAAA1 = [mean(W_ves_predicted,2) std(W_ves_predicted')'/sqrt(repeat_sim-1) mean(W_ves_actual,2) std(W_ves_actual')'/sqrt(repeat_sim-1)];
% % convert real data from origin (scatter) to here, then back again (to line) with AAAAAA
% i=1;
% for ii = 1:2:length(datas)
%     W_ves_predicted(1,i) = datas(ii,1);
%     W_ves_predicted(2,i) = datas(ii+1,1);
%     W_ves_actual(1,i) = datas(ii,2);
%     W_ves_actual(2,i) = datas(ii+1,2);
%     i = i+1;
% end
% 
% AAAAAA2 = [mean(Thresh_predicted,2) std(Thresh_predicted')'/sqrt(repeat_sim-1) mean(Thresh_actual,2) std(Thresh_actual')'/sqrt(repeat_sim-1)];
% % convert real data from origin (scatter) to here, then back again (to line) with AAAAAA
% i=1;
% for ii = 1:2:length(datas)
%     Thresh_predicted(1,i) = datas(ii,1);
%     Thresh_predicted(2,i) = datas(ii+1,1);
%     Thresh_actual(1,i) = datas(ii,2);
%     Thresh_actual(2,i) = datas(ii+1,2);
%     i = i+1;
% end

save(['likelihood_results_all_cells_2-17-10_40reps20repeats.mat']);
% save(['likelihood_results_all_cells_WichOFF_2-16-10.mat']);
% save(['likelihood_results_good_cells_only.mat']);
% save(['likelihood_results_good_cells_only_WichOFF.mat']);
% save(['likelihood_results_all_cells_wPrior_WichOff.mat']);
