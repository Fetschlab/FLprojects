function P =  images_dtb_2d(R)
% P =  images_dtb_2d(R)
% Method of images solution to bounded drift diffusion of two races with
% correlated (-1/sqrt(0.5)) noise  i.e.  2D drift to bound.
% The winner of the race is the one that reaches its upper bound first
% There is  no lower bound and the bounds are flat
% the grid size (ngrid) for simulation is set at 500 in DV space but can be changed
%
% ~~~~~~~~~~~~~~~
% Inputs is a structure R (all in SI units)
% drift:      vector of drift rates [length ndrift]
% t:          vector time series to simulate [length nt]
% Bup:        bound height [scalar]
% lose_flag:  whether we want the losing densities as well (slower)

% ~~~~~~~~~~~~~~~
% Outputs:  a structure P (which includes the inputs as well)
% P.up and P.lo (for two bounds) contains
%       pdf_t:  pdf of the winner [ndrift x nt]
%       cdf_t:  cdf of the winner [ndrift x nt]
%       p:  total prob up/lo  [ndrift]
%       mean_t:  mean time of up/lo [ndrift]
%       distr_loser: pdf of the losing race [ndrift x nt x ngrid ]
% y:    dv values of grid (bounds are at zero and intial dv is -Bup) [ngrid]
% dy:   grid spacing

if ~isfield(R,'lose_flag')
    R.lose_flag=0;
end
if ~isfield(R,'plotflag')
    R.plotflag=0;
end
P=R; % copy inputs to outputs

ndrift=length(R.drift); %number of drifts
nt=length(R.t); %number of time points to simulate
dt=R.t(2)-R.t(1); %sampling interval

ngrid=500; %grid size can change
g=linspace(-4*R.Bup,0,ngrid); %plausible range of dv values
dg=g(2)-g(1); %grid spacing

%  pdf of the losing races if needed
if R.lose_flag
    P.up.distr_loser = zeros(ndrift,nt,ngrid);
    P.lo.distr_loser = zeros(ndrift,nt,ngrid);
end

for id=1:ndrift %loop over drifts
    for k=2:nt  %loop over time starting at t(2)
        if R.lose_flag
            [P.up.pdf_t(id,k) P.lo.pdf_t(id,k) P.up.distr_loser(id,k,:) P.lo.distr_loser(id,k,:)]...
                = flux_img7(R.Bup,R.drift(id),R.t(k),g,g);
        else
            [P.up.pdf_t(id,k) P.lo.pdf_t(id,k) ]=flux_img7(R.Bup,R.drift(id),R.t(k),g,g); %#ok<*NCOMMA>
        end
    end
end
P.up.pdf_t= P.up.pdf_t*dt;
P.lo.pdf_t= P.lo.pdf_t*dt;

if R.lose_flag %put any missing density into the most negative bin
    P.up.distr_loser = P.up.distr_loser*dg*dt;
    P.lo.distr_loser = P.lo.distr_loser*dg*dt;
    
    P.up.distr_loser(:,:,1)= P.up.pdf_t- nansum(P.up.distr_loser(:,:,2:end),3);
    P.lo.distr_loser(:,:,1)= P.lo.pdf_t- nansum(P.lo.distr_loser(:,:,2:end),3);
end

P.up.cdf_t = cumsum(P.up.pdf_t,2);
P.lo.cdf_t = cumsum(P.lo.pdf_t,2);

P.up.p = P.up.cdf_t(:,end);
P.lo.p = P.lo.cdf_t(:,end);

try
    P.up.mean_t = P.up.pdf_t * R.t  ./P.up.p;
    P.lo.mean_t = P.lo.pdf_t * R.t  ./P.lo.p;
catch
    P.up.mean_t = P.up.pdf_t * R.t'  ./P.up.p;
    P.lo.mean_t = P.lo.pdf_t * R.t'  ./P.lo.p;
end

P.y=g;
P.dy=dg;

