function [indexHighWager] = calculating_MoreFlexibleTheta(logOddsMap, theta0, timePoint, DistOfOtherThres)

if timePoint > 3000 %ms %optimal
    indexHighWager = logOddsMap > theta0;
    if nargin == 4
        indexHighWager(logOddsMap > theta0+DistOfOtherThres) = false; % A-x-Time matrix
        % Clean up leftover
        for i = 1:size(indexHighWager, 1)
            p = find(indexHighWager(i,:) == 1);
            if ~isempty(p)
                diffP = diff(p);
                deleteP = find(diffP > 1, 1, 'last');
                deleteP = p(deleteP);
                indexHighWager(i, 1:deleteP) = false;
            end
        end
    end
elseif timePoint == 0 %flat
    indexHighWager = logOddsMap > theta0; % A-x-Time matrix
    [a, ~] = find(indexHighWager == 1);
    highPoint = min(a); % find the highest point
    indexHighWager(highPoint:end, :) = true;
else %hybrid 
    if 0 % Optimal at first but after a certain time just choose low!
        indexHighWager = logOddsMap > theta0; % A-x-Time matrix 
        % select where in time you want to flatten the curve
        indexHighWager(:, timePoint:end) = false; %vector from mat
    elseif 1 % Linear changing theta
        theta = theta0 + [1:size(logOddsMap, 2)].* timePoint; %-.001; %range -.005 to .005?
        indexHighWager = zeros(size(logOddsMap));  % A-x-Time matrix 
        for t = 1:size(logOddsMap,2)
            indexHighWager(:,t) = logOddsMap(:,t) > theta(t);
        end
        indexHighWager = logical(indexHighWager);
    elseif  false %(Optimal at first but after t=timePoint flat)
        indexHighWager = logOddsMap > theta0; % A-x-Time matrix 
        % select where in time you want to flatten the curve
        vectorIHW = indexHighWager(:, timePoint); %vector from mat
        indexOnes = find(vectorIHW == 1, 1); %find highest section of that vector with the value of 1
        for i = timePoint: size(indexHighWager,2) % only after the time
            indexHighWager(indexOnes:end, i) = true;
        end
    else % Flat at first, then optimal (How to model this?)
        indexHighWager = logOddsMap > theta0; % A-x-Time matrix 
        % select where in time you want to flatten the curve
        vectorIHW = indexHighWager(:, timePoint); %vector from mat
        indexOnes = find(vectorIHW == 1, 1); %find highest section of that vector with the value of 1
        for i = 1:timePoint % only before the time
            indexHighWager(1:(indexOnes-1), i) = false;
        end
    end
    
end
end