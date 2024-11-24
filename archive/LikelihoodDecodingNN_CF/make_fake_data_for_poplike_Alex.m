function [resp_mean_all, makefake_baselines, Wratio_opt] = make_fake_data_for_poplike_Alex

% clear
plotflag = 0;

% N = 41;
N = 32;
% N = 512;
hdg = -10:0.1:10;
heading = [-10 -3.5 -1.3 0 1.3 3.5 10];
coh= [.16 0.6];
nu_coh = length(coh);

% mean_vis = 15;
% slope_vis = 4;
% offset_vis = 15;
% offset_coh = 15;
% mean_ves = 30;
% slope_ves = 1.5;

% mean_vis = 60;
% slope_vis = 6;
% offset_vis = 0;
% offset_coh = 20; 
% mean_ves = 30;
% slope_ves = 2.5;

% % this matches the mean cong tunings very well:
% mean_vis = 15;
% slope_vis = 1.6;
% offset_vis = 19;
% offset_coh = -8;
% mean_ves = 14.8;
% slope_ves = .41;

% NO! that was incorrect! re-check plot_tunings! this is better:
mean_vis = 16;
slope_vis = 1.7;
offset_vis = 11;
offset_coh = 0;
mean_ves = 21;
slope_ves = .41;


f_1 = mean_vis*coh + offset_vis + offset_coh*(1-coh);
f_1_prime = slope_vis*coh;
f_2 = ones(1,nu_coh)*mean_ves;
f_2_prime = ones(1,nu_coh)*slope_ves;

% figure; plot(hdg,f_2_prime(1)*hdg + f_2(1),'k',hdg,f_1_prime(1)*hdg + f_1(1),'m',hdg,f_1_prime(2)*hdg + f_2(2),'r');

% weight_2 = (f_1.*f_2_prime)./(f_1_prime.*f_2);
    
for n = 1:N

	makefake_baselines(n) = mean([mean_ves mean_vis]);

    % randomize overall slopes, but relative slope across conditions is
    % fixed (all cells are congruent and should have same optimal weight)
%     randomize = rand-0.5;
%     randomize = rand/2+0.2 * sign(randn);
%     randomize = (rand+1) * sign(randn);
    
    randomize = sign(randn)*(randn+1);
    
% 	randomize = (randn+1);
    % .08, 08, .2, .4 = standard set
	% occasionally make one cue nearly flat, to simulate blowing up of Wopt
    if n < 5
        f_ves(n,:) = f_2_prime(1)*hdg*randomize*.08 + f_2(1) + round((rand-0.5)*15);
    else
        f_ves(n,:) = f_2_prime(1)*hdg*randomize + f_2(1) + round((rand-0.5)*15);
    end
    if n > 4 && n < 9
        f_vislow(n,:) = f_1_prime(1)*hdg*randomize*.08 + f_1(1) + round((rand-0.5)*15);
        f_vishigh(n,:) = f_1_prime(2)*hdg*randomize*.2 + f_1(2) + round((rand-0.5)*15);
    elseif n > 8 && n < 15
        f_vislow(n,:) = f_1_prime(1)*hdg*randomize*.2 + f_1(1) + round((rand-0.5)*15);
        f_vishigh(n,:) = f_1_prime(2)*hdg*randomize*.4 + f_1(2) + round((rand-0.5)*15);
    else
        f_vislow(n,:) = f_1_prime(1)*hdg*randomize + f_1(1) + round((rand-0.5)*15);
        f_vishigh(n,:) = f_1_prime(2)*hdg*randomize + f_1(2) + round((rand-0.5)*15);
    end    
    
    % shift up negative-going cells
    if min(f_ves(n,:)) < 0
        f_ves(n,:) = f_ves(n,:) - min(f_ves(n,:)) + 1;
    end
    if min(f_vislow(n,:)) < 0
        f_vislow(n,:) = f_vislow(n,:) - min(f_vislow(n,:)) + 1;
    end
    if min(f_vishigh(n,:)) < 0
        f_vishigh(n,:) = f_vishigh(n,:) - min(f_vishigh(n,:)) + 1;
    end

    % optimal weights from Alex's equation: f1*f2prime / f1prime*f2 ; REMEMBER CUE 2 IS VESTIBULAR IN ALEX'S STUFF!
    weight_2_low(n,:) = (f_vislow(n,2:end).*diff(f_ves(n,:))) ./ (diff(f_vislow(n,:)).*f_ves(n,2:end));
    weight_2_high(n,:) = (f_vishigh(n,2:end).*diff(f_ves(n,:))) ./ (diff(f_vishigh(n,:)).*f_ves(n,2:end));
    Wratio_opt{n} = [mean(weight_2_low(n,100)) mean(weight_2_high(n,100))]; % 65:135 is -3.5 to 3.5 hdg, 100 is zero
    
    a = 1; % this is somewhat arbitrary -- could also draw a random value
    f_comblow(n,:) = a*(median(weight_2_low(n,:))*f_ves(n,:) + f_vislow(n,:)); 
    f_combhigh(n,:) = a*(median(weight_2_high(n,:))*f_ves(n,:) + f_vishigh(n,:));
    
%     if min(f_comblow(n,:)) < 0
%         f_comblow(n,:) = f_comblow(n,:) - min(f_comblow(n,:)) + 1;
%     end
%     if min(f_vishigh(n,:)) < 0
%         f_combhigh(n,:) = f_combhigh(n,:) - min(f_combhigh(n,:)) + 1;
%     end
    
    if plotflag
%         figure(10);plot(hdg,f_ves(n,:),'k');hold on;
%         figure(11);plot(hdg,f_vislow(n,:),'m');hold on;
%         figure(12);plot(hdg,f_vishigh(n,:),'r');hold on;
%         and/or:
%         if n==1
            figure(102); clf; plot(hdg,f_ves(n,:),'k'); hold on;
            plot(hdg,f_vislow(n,:),'m');
            plot(hdg,f_vishigh(n,:),'r');
%             plot(hdg,f_comblow(n,:),'c');
%             plot(hdg,f_combhigh(n,:),'b');
            title(num2str([median(weight_2_low(n,:)) median(weight_2_high(n,:))]));
            pause;
%         end
    end
    
    for i = 1:length(heading)
        i_index = find(round(hdg*10) == round(heading(i)*10)); % fix matlab rounding problem
        i_index = i_index(1);
        resp_mean_all{n}{1,2,1}(i) = f_ves(n,i_index);
        resp_mean_all{n}{2,2,1}(i) = f_vislow(n,i_index);
        resp_mean_all{n}{2,2,2}(i) = f_vishigh(n,i_index);
%         resp_mean_all{n}{3,2,1}(i) = f_comblow(n,i_index);
%         resp_mean_all{n}{3,2,2}(i) = f_combhigh(n,i_index);
    end
end
% figure; plot(hdg,mean(f_vislow),'m',hdg,mean(f_ves),'k',hdg,mean(f_vishigh),'r');
% ylim([10 35]);

% n=7;
% n=16;
% figure; plot(hdg,f_vislow(n,:),'m',hdg,f_ves(n,:),'k',hdg,f_vishigh(n,:),'r');

centers = [-20 -10 -5 -2.5 -1.25 -.625 -.3125];
centers = [centers 0 fliplr(centers)*-1];
figure;
subplot(2,1,1);
bar(hist(weight_2_low(:,100),centers)); set(gca,'XTickLabel',centers)
subplot(2,1,2);
bar(hist(weight_2_high(:,100),centers)); set(gca,'XTickLabel',centers)

% centers = [-20 -10 -5 -2.5 -1.25 -.625 -.3125];
% centers = [centers 0 fliplr(centers)*-1];
% figure;
% subplot(2,1,1);
% bar(hist(asdf(1:2:63),centers)); set(gca,'XTickLabel',centers)
% subplot(2,1,2);
% bar(hist(asdf(2:2:64),centers)); set(gca,'XTickLabel',centers)

return




