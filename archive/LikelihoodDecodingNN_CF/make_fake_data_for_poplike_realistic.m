function [resp_mean_all, makefake_baselines] = make_fake_data_for_poplike_realistic

N = 48;
hdg = -10:0.1:10;
g = 1;
slope_start = -0.48;
slope_end = 0.48;
plotflag = 0;

% these values are from fitting mean tuning curves
mult_low = 0.62; % 0.62 normally
mult_high = 1.96; % 1.96 normally

m_ves = slope_start: (slope_end-slope_start)/(N-1) : slope_end;
m_vislow = (slope_start: (slope_end-slope_start)/(N-1) : slope_end) * mult_low;  
m_vishigh = (slope_start: (slope_end-slope_start)/(N-1) : slope_end) * mult_high;

N = length(m_ves);
heading = [-10 -3.5 -1.3 0 1.3 3.5 10];

for n = 1:N
    % to create a more realistic population, (1) randomize baselines...
%    spons = gamrnd([3 3 3],[4 3 4.5]);
%	 spons = gamrnd([3 3 3],[6 3 9]); % a more drastic difference in DC across conditions
    spons = [25 10 60]; % even more drastic
    spons(spons<5)=5;
    spon_ves = spons(1);
    spon_vislow = spons(2);
    spon_vishigh = spons(3);
    makefake_baselines(n) = mean(spons); % this is a kluge -- the conditions act as though they have separate baseline firing rates, but obviously a cell has only one, so I'll just take the mean
    
    f_ves(n,:) = m_ves(n)*hdg + spon_ves;
    
    % ... (2) randomize slopes a little bit ...
    f_vislow(n,:) = m_vislow(n)*hdg*(rand+0.5) + spon_vislow;
    f_vishigh(n,:) = m_vishigh(n)*hdg*(rand+0.5) + spon_vishigh;
    
    % ... and (3) sprinkle in a few higher-slope cells for good measure
    if round(n/8) == n/8
        f_ves(n,:) = m_ves(n)*hdg*(rand+1) + spon_vislow+8;
    end
    if round(n/6) == n/6
        f_vislow(n,:) = m_vislow(n)*hdg*(rand+1) + spon_vislow+6;
        f_vishigh(n,:) = m_vishigh(n)*hdg*(rand+1.75) + spon_vishigh+10;
    end
    % shift up negative-going cells
    if min(f_vislow(n,:)) < 0
        f_vislow(n,:) = f_vislow(n,:) + min(f_vislow(n,:)) + 1;
    end
    if min(f_vishigh(n,:)) < 0
        f_vishigh(n,:) = f_vishigh(n,:) - min(f_vishigh(n,:)) + 1;
    end
    
    if plotflag
        figure(10);plot(hdg,f_ves(n,:),'k');hold on;
        figure(11);plot(hdg,f_vislow(n,:),'m');hold on;
        figure(12);plot(hdg,f_vishigh(n,:),'r');hold on;
    end
    
    for i = 1:length(heading)
        i_index = find(round(hdg*10) == round(heading(i)*10)); % fix matlab rounding problem
        i_index = i_index(1);
        resp_mean_all{n}{1,2,1}(i) = f_ves(n,i_index);
        resp_mean_all{n}{2,2,1}(i) = f_vislow(n,i_index);
        resp_mean_all{n}{2,2,2}(i) = f_vishigh(n,i_index);
    end
end
% figure; plot(hdg,sum(f_vislow));

% % temp: run plot_tunings code to see how these differ from real data
% for n = 1:N
%     [rr,pp] = corrcoef(heading, resp_mean_all{n}{2,2,2});
%     if rr(1,2) < 0
%         tuning_vestib(n,:) = fliplr(resp_mean_all{n}{1,2,1});
%         tuning_vislow(n,:) = fliplr(resp_mean_all{n}{2,2,1});
%         tuning_vishigh(n,:) = fliplr(resp_mean_all{n}{2,2,2});
% %         tuning_comblow(n,:) = fliplr(resp_mean_all{n}{3,2,1});
% %         tuning_combhigh(n,:) = fliplr(resp_mean_all{n}{3,2,2});
%     else
%         tuning_vestib(n,:) = resp_mean_all{n}{1,2,1};
%         tuning_vislow(n,:) = resp_mean_all{n}{2,2,1};
%         tuning_vishigh(n,:) = resp_mean_all{n}{2,2,2};
% %         tuning_comblow(n,:) = resp_mean_all{n}{3,2,1};
% %         tuning_combhigh(n,:) = resp_mean_all{n}{3,2,2};
%     end
% end
% 
% figure; hold on;
% plot(heading, mean(tuning_vestib), 'k', 'LineWidth', 3);
% plot(heading, mean(tuning_vislow), 'm', 'LineWidth',  3);
% plot(heading, mean(tuning_vishigh), 'r', 'LineWidth', 3);
% legend('Vestib','Vis-Low','Vis-High','Location','NorthWest');
% xlabel('Heading Angle (deg)'); ylabel('Firing Rate (sp/s)'); xlim([-11 11]);
% 
% % fit mean tuning curve to find parameters for make_fake_data
% P = polyfit(heading,mean(tuning_vestib),1);
% slope_vestib = P(1)
% P = polyfit(heading,mean(tuning_vislow),1);
% slope_vislow = P(1)
% P = polyfit(heading,mean(tuning_vishigh),1);
% slope_vishigh = P(1)
% 
% slope_ratio_low = slope_vislow/slope_vestib
% slope_ratio_high = slope_vishigh/slope_vestib

return




