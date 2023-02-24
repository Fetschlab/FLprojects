% %%

% options.fitMethod = 'fms';
% options.fitMethod = 'global';
% options.fitMethod = 'multi';
% options.fitMethod = 'pattern';
% options.fitMethod = 'bads';


% %% check sample sizes for each trial type
% 
% % reminder:
% % n = nan(length(mods),length(cohs),length(deltas)+1,length(hdgs));
% %                                % add extra column^ for pooling all trials irrespective of delta
% 
% % could reshape the 4D array to a column vector using sort(n(:)), but need
% % to know the dim ordering.
% 
% % conceptually easier (if inelegant) to create arrays of same size as N
% % which jointly identify the unique trial type. Then we can reshape them
% % the same way and pass in the index array from sort(n(:)) to get trial type.
% for m = 1:length(mods)
% for c = 1:length(cohs)
% for d = 1:length(deltas)+1 % add extra column for all trials irrespective of delta
% for h = 1:length(hdgs)
%     Mod(m,c,d,h) = m;
%     Coh(m,c,d,h) = c;
%     Delta(m,c,d,h) = d;
%     Hdg(m,c,d,h) = h;    
% end
% end
% end
% end
% Mod = Mod(:); Coh = Coh(:); Delta = Delta(:); Hdg = Hdg(:);
% 
% [ntrSorted,ind] = sort(n(:));
% 
% % first verify that indices with zero trials are only the invalid
% % combinations:
% zeroInds = ind(ntrSorted==0);
% conds = [Mod(zeroInds) Coh(zeroInds) Delta(zeroInds) Hdg(zeroInds)];
% condsWithZeroTr = sortrows(conds,[1 3]) % this should show that all zeroInd conds
%                                         % are mod 1 or 2 and delta 1 or 3
%                               
% % now look at the conds with fewest trials to see if something's amiss
% lowInds = ind(ntrSorted>0 & ntrSorted<30);
% conds = [Mod(lowInds) Coh(lowInds) Delta(lowInds) Hdg(lowInds)];
% condsWithLowN = sortrows(conds,[1 2 3 4]) % a random sprinkling of low-coh tr...
%                                           % double check pldaps code
% 
% 
% 
% ks = 17;
% sigma = 0.03;
% B = 1.5;
% 
% guess = [ks sigma B];
% 
% % plot error trajectory (prob doesn't work with parallel fit methods)
% options.ploterr = 0;
% options.fh=500;


%% fit DDM

% options.fitMethod = 'fms';
% options.fitMethod = 'global';
% options.fitMethod = 'multi';
% options.fitMethod = 'pattern';

% params: 

% % % 
% % % 
% % % % Drugowitsch model has 12 params, plus 8 for biases and lapse rates
% % % % we'll skip the latter, and we can drop the 3 Tnd terms until we are fitting RT
% % % % so we have 9 params:
% % % 
% % %     %    aVis gammaVis bVis thetaVis kVes thetaVes gammaCom bCom thetaCom
% % % fixed = [0    0        0    0        0    0        0        0    0       ];
% % % 
% % % % per Drugowitsch, variance of momentary evidence scales with coherence as:
% % % % var(c) ~ 1 + bVis*coh^gammaVis;
% % % 
% % % % similarly,
% % % % kVis ~ aVis*cVis^gammaVis
% % % 
% % % % the bound gets a free parameter (theta) for each modality, but this is
% % % % only because the variance is what really changes across conditions, and
% % % % this is absorbed into the definition of (normalized) bounds. I'm still
% % % % not sure about this...
% % %  
% % % % initial guess (or hand-tuned params)
% % % aVis = 0.5; % sensitivity parameter, multiplies coh
% % % gammaVis = 1; % determines scaling of sensitivity (and variance) by coh
% % % bVis = 30; % bound height
% % % theta = 1.6; % criterion (in log odds correct) for betting high
% % % alpha = 0.1; % base rate of low-bet choices
% % % 
% % % guess = [aVis gammaVis bVis thetaVis kVes thetaVes gammaCom bCom thetaCom];
% % % 
% % % 
% % % 
% % % 

options.RTtask = 1;
options.conftask = 1;

options.fitMethod = 'fms';
% options.fitMethod = 'global';
% options.fitMethod = 'multi';
% options.fitMethod = 'pattern';
% options.fitMethod = 'bads';

% params that work well for Dots
% % % k = 7;
% % % B = 1.5;
% % % theta = 1;
% % % alpha = 0.2; % base rate of low-bet choices
% % % Tnd = 0.3; % non-decision time
% % % sigma = 0.1; (not a free param)

% old params that don't work here
% % % kves = 1.2;
% % % kvisMult = 4; % will be multiplied by coh to get kvis (this simplifies parameterization)
% % % B = 70;
% % % muTnd = 0.3; % s


% initial guess (or hand-tuned params)
% ks = 28;
% sigma = 0.05;
% B = 3;
% muTnd = 0.5;
%     % this isn't bad but converges to what must be a local min

ks = 30;
sigma = 0.05;
B = 3;
muTnd = 0.3;

    %    ks sigma B  muTnd
fixed = [0  0     0  0   ]; % can be used to fix some params and not others
guess = [ks sigma B muTnd];

% ************************************
% set all fixed to 1 for hand-tuning:
fixed(:)=1;
% ************************************

% plot error trajectory (prob doesn't work with parallel fit methods)
options.ploterr = 0;

close all


[X, err_final, fit, fitInterp] = dots3DMP_fitDDM(data,options,guess,fixed);

% plot it!
fitgauss = 1; % 'smooth' choppy model predictions by fitting (c)gaussians

dots3DMP_plots_fit(data,fitInterp,conftask,RTtask,fitgauss)



