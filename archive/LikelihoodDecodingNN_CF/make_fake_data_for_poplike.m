function [resp_mean_all, makefake_baselines] = make_fake_data_for_poplike

N = 48;
hdg = -10:0.1:10;
g = 1;
slope_start = -0.48;
slope_end = 0.48;
plotflag = 0;

m_ves = slope_start: (slope_end-slope_start)/(N-1) : slope_end;
m_vislow = (slope_start: (slope_end-slope_start)/(N-1) : slope_end) * 0.2;  % * 0.62;  % these values are from 
m_vishigh = (slope_start: (slope_end-slope_start)/(N-1) : slope_end) * 5;   % * 1.96; % fitting mean tuning curves

N = length(m_ves);
spon_ves = 50;
spon_vislow = 50;
spon_vishigh = 50;
makefake_baselines(1:N) = 50;

heading = [-10 -3.5 -1.3 0 1.3 3.5 10];

for n = 1:N
    f_ves(n,:) = m_ves(n)*hdg + spon_ves;
    
    % EFFECT OF COHERENCE:
    % (1) Slope change:
    f_vislow(n,:) = m_vislow(n)*hdg + spon_vislow;
    f_vishigh(n,:) = m_vishigh(n)*hdg + spon_vishigh;
%     % (2) Ma style (gain)
%     f_vislow(n,:) = f_ves(n,:) * 0.5;
%     f_vishigh(n,:) = f_ves(n,:) * 2;

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

return
