% dots3DMP test likelihoods for each parameter

% for each parameter, hold all others fixed at generative parameters, and see how LL for each of the
% three outcomes (choice, RT, conf)


%% sim data

cd /Users/stevenjerjian/Desktop/FetschLab/Analysis/data/dots3DMP_DDM
load tempsim_sepConfMaps_1000reps.mat

%% set some vars

% set allNonHB == 0 to ignore unabsorbed probability, ie no need for weighted sum with max_dur
% in RT calculation, e.g. when sim excluded trials that failed to hit bound
if ~exist('allowNonHB','var'); allowNonHB=0; end

mods   = unique(data.modality);
cohs   = unique(data.coherence);
deltas = unique(data.delta);

options.RTtask   = 1; 
options.conftask = 2; % 1=continuous/rating, 2=PDW


%% select model + options

% ==== method ====
modelID=1; % unused for now
options.errfcn    = @dots3DMP_errfcn_DDM_2D_wConf_noMC; % 2D DDM aka anticorrelated race, for RT+conf [Kiani 14 / van den Berg 16 (uses Wolpert's images_dtb_2d (method of images, from Moreno-Bote 2010))]
options.fitMethod = 'global'; %'fms','global','multi','pattern','bads'
options.whichFit  = {'choice','conf','RT'}; % choice, conf, RT, multinom (choice+conf)

% ==== implementation ====
options.dummyRun = 0; % dummyRun=1
options.ignoreUnabs = ~allowNonHB;
options.useVelAcc = 0; % use stimulus physical profiles in drift rates
% options.confModel = 'evidence+time'; % unused for now

% ==== diagnostics and output ====
options.plot      = 0;      % plot confidence maps
options.feedback  = 1;      % 0 - none, 1 - text feedback on LLs/err
options.runInterpFit = 0;   % model predictions for interpolated headings? for nice plotting

paramNames = {'kmult','B',sprintf('%s_Ves',char(952)),sprintf('%s_Vis',char(952)),sprintf('%s_Comb',char(952)),...
    'alpha','TndVes','TndVis','TndComb'};

guess = [origParams.kmult, origParams.B, origParams.theta, origParams.alpha, origParams.TndMean/1000];
numParams = length(guess);

fixed = zeros(1,numParams);
fixed(:) = 1;




%% evaluate the model at these parameters

nTestVals = 10;
testvals(:,1) = linspace(10,100,nTestVals);
testvals(:,2) = linspace(0.5,1.5,nTestVals);
testvals(:,3) = linspace(0.5,1.5,nTestVals);
testvals(:,4) = linspace(0.5,1.5,nTestVals);
testvals(:,5) = linspace(0.5,1.5,nTestVals);
testvals(:,6) = linspace(0.0,0.1,nTestVals);
testvals(:,7) = linspace(0.0,0.9,nTestVals);
testvals(:,8) = linspace(0.0,0.9,nTestVals);
testvals(:,9) = linspace(0.0,0.9,nTestVals);


all_LL = nan(nTestVals,numParams,4);
for p = 1:numParams
    % p is the param we are going to vary
    
    X = guess;
    
    for n = 1:nTestVals
        X(p) = testvals(n,p);
        fprintf('varying %s, current = %.1f\n',paramNames{p},X(p))
        [err_final, ~, ~, LLs] = feval(options.errfcn,X,X,fixed,data,options);
        all_LL(n,p,1) = LLs.choice;
        all_LL(n,p,2) = LLs.conf;
        all_LL(n,p,3) = LLs.RT;
        all_LL(n,p,4) = err_final;
%         fprintf('err_final = %.2f\n\n',err_final)

    end

end


%%
titles = {'Choice','Conf','RT','ALL'};

figure('position',[100 100 1200 1200])
for p=1:9
for c=1:4
    row = p; col = c;
    subplot(9,4,c+4*(p-1)); hold on
    if c<4
        y = -all_LL(:,p,c);
    else
        y = all_LL(:,p,c);
    end
    y = y-min(y);
    plot(testvals(:,p),y,'bo','linestyle','-');
    plot([guess(p) guess(p)],get(gca,'ylim'),'r--')
    if p==1
        title(titles{c});
    end
    if c==1
        ylabel(paramNames{p});
    end
end
end

%% what about varying 2 params at the same time...let's just do k and B



