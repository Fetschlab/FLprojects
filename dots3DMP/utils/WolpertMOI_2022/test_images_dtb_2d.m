% P =  images_dtb_2d(R)
% Method of images solution to bounded drift diffusion of two races with
% correlated (-1/sqrt(0.5)) noise  i.e.  2D drift to bound.
% The winner of the race is the one that reaches the upper bound first
% There a no lower bound and the bounds are flat
%
% ~~~~~~~~~~~~~~~
% Inputs is a structure R (all in SI units)
% drift     vector of drift rates [length ndrift]
% t         vector time series to simulate in seconds [length nt]
% Bup       bound height [scalar]
% k_urg    time dependent urgencey signal added to DV
% grid     grid for x and y
%
clear all; close all;
R.notabs_flag = 0;
clf
for k=1 %:2
    kappa=12;
    coh=[0 0.032 0.064 0.128 0.256 0.512];
    R.drift=kappa*coh;
    R.t=linspace(0,2,500)';
    R.Bup=.7;
    R.k_urg=0;
    R.lose_flag=0;  %do we need the pdf of the losing race - usually not
    R.grid=linspace(-7,0,500);  %dv values can change lower

    if k==1
        R.low_th=-Inf;  % this is the lower threshold
    elseif k==2
        R.low_th=-R.Bup-R.Bup/4;  % this is the lower threshold
    else
        R.low_th=-R.Bup-R.Bup/4;  % this is the lower threshold
    end

    tnd=0.0;

%     P =  images_dtb_2d_new (R);
    P =  images_dtb_2d(R);

    if 1
        figure(k);

        subplot(2,2,1)
        plot(1:6,P.up.p,'ko-');
        hold on
        plot(1:6,P.lo.p,'ro-');
        plot(1:6,P.up.p+P.lo.p,'go-');
        ylabel('P.up (blk), P.lo (red)')
        xlabel('Drift')

        subplot(2,2,2)
        plot(1:6,P.up.mean_t+tnd,'ko-');
        hold on
        plot(1:6,P.lo.mean_t+tnd,'ro-');
        xlabel('Drift')
        ylabel('RT (s)')

        subplot(2,2,3)
        plot(P.t,P.up.pdf_t,'k');
        hold on
        plot(P.t,P.lo.pdf_t,'r');
        aa=axis;
        axis([0 P.t(end) aa(3:4)])
        ylabel('pdf')
        xlabel('t (s)')

        subplot(2,2,4)
        plot(P.t,P.up.cdf_t,'k');
        hold on
        plot(P.t,P.lo.cdf_t,'r');
        aa=axis;
        axis([0 P.t(end) aa(3:4)])
        ylabel('cdf')
        xlabel('t (s)')

        hold on
    end
    
end




% 
% 
%     %%
%     figure(5)
%     
%     if k==1
%         subplot(3,2,1)
%         plot(P.t,P.up.notabs.mean','k-');
%         hold on
%         plot(P.t,P.lo.notabs.mean','r-');
%         ylabel('Mean')
%         xlabel('t')
%         axis([0 0.2 -1.8 0])
% 
%     else
% 
%         subplot(3,2,2)
%         plot(P.t,P.up.notabs.mean','k--');
%         hold on
%         plot(P.t,P.lo.notabs.mean','r--');
%         ylabel('Mean')
%         xlabel('t')
%         axis([0 0.2 -1.8 0])
% 
%     end
% 
%     subplot(3,2,3)
%     plot(P.t,P.up.notabs.var','k-');
%     hold on
%     plot(P.t,P.lo.notabs.var','r-');
%     ylabel('Var')
%     xlabel('t')
% 
% 
% 
%     figure(6)
%     plot(P.y,squeeze(P.lo.distr_loser(end,100,:))','k');
% 
% 
%     shg
