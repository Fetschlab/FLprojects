function [pUpAbs, rtUp, rtLo, upDist, loDist,pLoAbs,dvNotAbs,pGT0,xtDist] = runFPonSeries(uvect,t,Bup,Blo,y0,sigma,urgencyMax)
dbstop if error

% Runs the Fokker-Planck to bounds through a series of stimulus difficulties.
% Just a wrapper to run through a coherence series to generate
% probability surface in time under boundary conditions.

% Original author: Michael Shadlen (mns)

% Input args
% ~~~~~~~~~~
% uvect vector of the mean drift rates
% t     time axis. Typically [0:dt:tmax] where dt is 0.05 and tmax is
%       whatever you choose
% Bup and Blo are the upper and lower bounds (same dimension as t)
% y0    starting value of the accumulation; zero by default.
%
% dependencies: FP4.mex (from Roozbeh Kiani) or a slightly modified version
% FP4_BoundCrossZero.mex. I think either will work now because roozbeh made
% some adjustments fto FP4.
%
% return args
% ~~~~~~~~~~~
% pUpAbs  for each x member of uvect, prob of absorption at upper bound
% rtUp    for each x member of uvect, the mean time of absorption at upper bound
% rtLo    for each x member of uvect, the mean time of absorption at lower bound
% upDist  for each x member of uvect, the probability distribution of bound
%         crossings at upper bound (cell array of vectors size of t)
% loDist  for each x member of uvect, the probability distribution of bound
%         crossings at upper bound (cell array of vectors size of t)
% pLoAbs  see pUpAbs (see below)
% dvNotAbs  for each x member of uvect, the distribution of the
%           non-absorbed decision variable at t=tmax. This is a column
%           vector that corresponds to the values of xmesh, which is the
%           grid formed within the FP routine. It's not available, but I
%           haven't needed it.
% pGT0     for each x member of uvect, a vector that is the total
%           probability above 0 at each time step. I think this includes
%           absorbed probability. It's useful for calculations in non-RT
%           experiments when a decision is based on the sign of the dv
%           (cell array of vectors size of t)
% xtDist   for each x member of uvect, a matrix containing the
%          probabability of the dv at each time. I can't remember what's
%          stored at the bounds. This is a huge matrix because its
%          dimensions are size of t by the mesh that separates the upper
%          and lower bounds. (cell array of matrices)


% I added the pLoAbs return argument later mainly to help cross check that
% the sum of the total choice probability is 1. 
% 1 - rtUp(i) + rtLo(i) should be equal to the total probabability of
% nonAbsorbed decision variable at tmax, which is sum(dvNotAbs{i}) 

% 5/30/08 mns wrote it
% 12/13/08 mns changed pcor to pUpAbs; added pLoAbs and dvNotAbs
% 6/23/09 mns added pGT0 argout
% 8/26/09 updated documentation

m = length(uvect); %MVL: Length of the Mean Drift Rate for each Coherence (So length of #of Coh)

% if any(Bup <=0) || any(Blo>=0)
%     warning('bounds cross zero')
% end

% Don't need next bit if using RK's method to generate xmesh
% length of mesh, ensure it's odd
% nmesh = round((2*maxB + diff(b_margin))/dx);
% if rem(nmesh,2)==0
%     nmesh = nmesh+1;
% end
% xmesh = linspace(minBlo-b_margin(1), maxBup+b_margin(2), nmesh)';

% length(xmesh)

% Initialization of the probability surface
if nargin<5
    y0 = 0;
end

% Roozbeh wants to see bounds passed as a single scalar from which the mesh
% is constructed. The real upper bound is the 2nd column of b_change added
% to B. The real lower bound is the 1st column of b_change added to -B.
% b_change = [Blo(:) - minBlo, Bup(:)-maxBup]; 
% figure(10), plot(t,b_change(:,1),t,b_change(:,2))


%allocate vectors and cells 
%MVL: I think 'nans' is suppose to be 'NaN'
pUpAbs = NaN(m,1); %MVL: Probability of hitting upper bond
pLoAbs = NaN(m,1); %MVL: Probability of hitting lower bond
rtUp = NaN(m,1); %MVL: Time in which upper bound was hit
rtLo = NaN(m,1); %MVL: Time in which lower bound was hit
upDist = cell(m,1); %MVL: Probability distribution of Upper bound crossing for coherence -> notice cell not variable
loDist = cell(m,1); %MVL: Probability distribution of Lower bound crossing for coherence
dvNotAbs = cell(m,1);
if nargout>7
    pGT0 = cell(m,1);
end
if nargout > 8
    xtDist = cell(m,1); %MVL: The most important matrix (Holds all info)
end

for i = 1:m %MVL: Iterates through Mean Drift Coefficient (# of Coherences)
    % All arguments to FP4 are corrupted by FP4. Reset them w/in the loop,
    % before each call.
    
    % time
    dt  = t(2)-t(1); %MVL: We already stored dt, is it done like this because of the corruption
%     max_dur = max(t);
    
    % ygrid, what Roozbeh calls the mesh
%     sigma = 1;
    maxBup = max(Bup(:)); %MVL: Max upper bound value
    minBlo = min(Blo(:)); %MVL: Min lower bound value
    maxB = max(abs([maxBup minBlo])); %MVL: Max bound value overall
    % B = maxB; % temporary to use RK's nomenclature
    dx = min(0.1,maxB/100); %MVL: highest value for dx bound, no higer than .1
    b_margin = [4*sigma; 4*sigma]; %MVL: Variance of Drift Rate * 4, for some reason. Why in this vector form? Most likely because of matrix multiplication purposes.
    xmesh = (minBlo-b_margin(1)+dx : dx : maxBup+b_margin(2)-dx)'; % RK's method %MVL: Constructing the xmesh, or big matrix (I think its actually a vector) that will contian probability values
    % adjust the mesh so that it contains a zero.
    if ~any(xmesh==0) %MVL: Check if xmesh has any nonzeros 
        % adjust dx so that the mesh has zero and is therefore
        % symmetric
        [~, I] = min(abs(xmesh)); %MVL: Find smallest value in xmesh
        I = min(I(:)); % in case there's a tie
        % xmesh(I) is closest to 0.
        b_margin = b_margin+xmesh(I); %MVL: Grab the smallest value in the xmesh and add by the b_margin values and make that the new b_margin values
        xmesh = (minBlo-b_margin(1)+dx : dx : maxBup+b_margin(2)-dx)'; % recompute the mesh
        xmesh(abs(xmesh)<1e-12) = 0; % now and then you don't get exactly zero
        % check to be sure that 0 is in the mesh
        if sum(xmesh==0) ~= 1
            save crashfile.mat
            error('My logic is wrong in the xmesh adjustment. Saving crashfile.mat')
        end
    end

    %MVL: 11/28/23 I might be wrong but this might be for dealing with collapsing bounds
    b_change = [Blo(:) - minBlo, Bup(:)-maxBup];  % Roozbeh wants this argument

    % Start position
    % MVL: Matrix that contains a value of 1 for where the starting point
    % is
    uinit = zeros(size(xmesh));
    [~, I] = min(abs(xmesh-y0));
    uinit(I) = 1;   % this should be a 1 where the mesh is closest to y0
    if sum(uinit)~=1
        error('defective initialization or mesh')
    end
        
    % [ufinal_,Pt_,Ptb_] = FP4(xmesh,uinit,uvect(i),sigma,b_change,b_margin,dt);
    %MVL: Inputs: xmesh which is the vector containing values between Up
    %Boundary and Low Boundary with steps dx. Keep in mind its the Y
    %values. uinit vector with starting point as 1, unbiased should be
    %right in the middle. uvect - drift mean values, iterate through this.
    %sigma contains variance of drift rate. b_change difference between min
    %Bound and Low Bound vector, relevant when Bounds change through
    %trials. b_margin variance of drift rate times 4 (why?). dt - time step
    %in x axis.
    if nargout<8
        [ufinal,~,Ptb] = FP4(xmesh,uinit,uvect(i),sigma,b_change,b_margin,dt);
    elseif nargout<9
        [ufinal,~,Ptb,Pg0] = FP4(xmesh,uinit,uvect(i),sigma,b_change,b_margin,dt);
    else
        try 
            [ufinal,~,Ptb,Pg0,Pxt] = FP4(xmesh,uinit,uvect(i),sigma,b_change,b_margin,dt); 
        catch
            % So the problem arises because b_margin appears to be to
            % similar to a bound, or maybe too small so when this happens
            % make b_margin = [4; 4];
            try 
                b_margin = [4*sigma; 4*sigma];
                [ufinal,~,Ptb,Pg0,Pxt] = FP4(xmesh,uinit,uvect(i),sigma,b_change,b_margin,dt); 
            catch % This next part seems to work
                try % This should only work when the grid is on the edge
                    % temp = xmesh(1) > xmesh(end);
                    temp = abs(xmesh(1)) > xmesh(end);
                    if ~temp
                        b_margin = abs([xmesh(3); xmesh(3)]);
                    else
                        b_margin = [xmesh(end-2); xmesh(end-2)];
                    end
                    [ufinal,~,Ptb,Pg0,Pxt] = FP4(xmesh,uinit,uvect(i),sigma,b_change,b_margin,dt); 
                catch %if another error occurs, god help me!
                    keyboard
                end
            end
        end
    end

    if nargin > 6 && urgencyMax > 0
        % First outline how the linear bound works
        newBup = maxBup - linspace(0, urgencyMax, size(t,2));
        newBlo = minBlo + linspace(0, urgencyMax, size(t,2));
        % Prevent from it crossing zero
        newBup(newBup <= 0) = newBup(find(newBup<=0, 1)-1); %.001;
        newBlo(newBlo >= 0) = newBlo(find(newBlo>=0, 1)-1); %-.001;
        % Now calculate the crossing of bound from these new bounds
        newPtb = zeros(length(newBup), 2);
        for ti = 2:length(newBup) % Loop through time
            % Now find the components of the bound that match with the mesh
            newPtb(ti,2) = sum(Pxt(find(xmesh >= newBup(ti), 1):end, ti));
            newPtb(ti,1) = sum(Pxt(1:find(xmesh <= newBlo(ti), 1, 'last'), ti));
        end
        % Now add it with the rest of things
        % Original code from Kiani
        plo = Ptb(2:end,1) + newPtb(:,1);        %probability of crossing the lower bound at each moment
        pup = Ptb(2:end,2) + newPtb(:,2);        %probability of crossing the upper bound
        % everything else looks the same
        pUpAbs(i) = sum(pup);
        pLoAbs(i) = sum(plo);
        % These means only make sense if there is no mass in ufinal. It's the
        % mean rt only for absorbed.
        rtUp(i) = sum(pup(1:end) .* t(:)) ./ pUpAbs(i);
        rtLo(i) = sum(plo(1:end) .*t(:)) / pLoAbs(i);
        upDist{i} = pup(1:end); % in RK's code, the 1st element of Ptb is absorption at t=0, that is 1 time step before our 1st
        loDist{i} = plo(1:end); % Because I have accounted this above start at 1 instead of the usual 2
        dvNotAbs{i} = ufinal;
        if exist('Pg0','var')
            pGT0{i} = Pg0;
        end
        if exist('Pxt','var')
            xtDist{i} = Pxt;
        end
    else % The above only works for the sequential model, everything else... No
        % Original code from Kiani
        plo = Ptb(:,1);         %probability of crossing the lower bound at each moment
        pup = Ptb(:,2);        %probability of crossing the upper bound
        % survival_ = Pt_;            %survivor function
        % Pt_ = -diff(Pt_);           %the total probability of crossing either of the bounds at each moment
        pUpAbs(i) = sum(pup);
        pLoAbs(i) = sum(plo);
        % These means only make sense if there is no mass in ufinal. It's the
        % mean rt only for absorbed.
        rtUp(i) = sum(pup(2:end) .* t(:)) ./ pUpAbs(i);
        rtLo(i) = sum(plo(2:end) .*t(:)) / pLoAbs(i);
        upDist{i} = pup(2:end); % in RK's code, the 1st element of Ptb is absorption at t=0, that is 1 time step before our 1st
        loDist{i} = plo(2:end);
        dvNotAbs{i} = ufinal;
        if exist('Pg0','var')
            pGT0{i} = Pg0;
        end
        if exist('Pxt','var')
            xtDist{i} = Pxt;
        end
    end
        
end

