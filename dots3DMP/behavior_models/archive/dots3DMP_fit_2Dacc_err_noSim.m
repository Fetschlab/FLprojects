function [err,fit] = dots3DMP_fit_2Dacc_err_noSim(param, guess, fixed, data, options)

% SJ 10-11-2021 no Monte Carlo simulation for fitting, it's redundant! just use model
% predictions directly

global call_num %paramVals errVals

if mod(call_num,100)==0
    plot_flag=1;
else
    plot_flag=0;
end

% set parameters for this run based on guess and fixed params flag
param = getParam(param, guess, fixed);

% data.heading = abs(data.heading);

mods   = unique(data.modality)';
cohs   = unique(data.coherence)'; 
hdgs   = unique(data.heading)';
deltas = unique(data.delta)';

% assume symmetrical
ushdgs = hdgs(hdgs>=0);

% don't 'fit' deltas unless dummyRun is set (for predicted curves)
if options.dummyRun==0
%     mods = [1 2]; % don't fit combined condition at all!
    deltas = 0; 
end

duration = 2000; % stimulus duration (ms)

kves  = param(1)     .*100;
kvis  = param([2 3]) .*100;

if ~options.sepbounds
    % single bound
    B     = abs(param(4)); % don't accept negative bound heights
    Tnds  = param(5:7);
    ttc      = param(8);
    
    if options.conftask==1
        cL = param(9:11);
    elseif options.conftask==2
        data.conf = data.PDW;
        theta  = param(9);
        cL = param(10:12);
    else
        cL = [0 0 0];
    end
else
    % separate bounds
    B     = abs(param(4:6)); % don't accept negative bound heights
    Tnds  = param(7:9);
    ttc      = param(10);
    
    if options.conftask==1
        cL = param(9:11);
    elseif options.conftask==2
        data.conf = data.PDW;
        theta  = param(9);
        cL = param(10:12);
    else
        cL = [0 0 0];
    end
end

paramNames = options.paramNames;

duration = duration + ttc;
% sdTnd = 60; % fixed SD

% assume the mapping for confidence is based on an equal amount of experience with the 
% *three* levels of reliability (ves, vis-low, vis-high) hence k is the
% mean 
% need separate logOddsMaps if bound heights are different

k = mean([kves kvis]);

clear logOddsMap_All
for b=1:length(B)
    if b>1,plot_flag=0; end
    R.t = 0.001:0.001:duration/1000;
    R.Bup = B(b);
    R.drift = k * sind(ushdgs); % takes only unsigned drift rates
    R.lose_flag = 1;
    R.plotflag = plot_flag; % 1 = plot, 2 = plot and export_fig
    P =  images_dtb_2d(R);
    logOddsMap_All(:,:,b) = P.logOddsCorrMap;
end
if length(B)==1
    logOddsMap_All = repmat(logOddsMap_All,1,1,3);
    B = [B B B];
end

% now compute the P we need for model choices and RTs (modality-specific)
% use signed headings here!

RVes.t = 0.001:0.001:duration/1000;
RVes.Bup = B(1);
RVes.drift = kves * sind(hdgs); 
RVes.lose_flag = 1;
RVes.plotflag = plot_flag; % 1 = plot, 2 = plot and export_fig
PVes =  images_dtb_2d(RVes);

clear PVis PComb
for c = 1:length(cohs)
    RVis.t = 0.001:0.001:duration/1000;
    RVis.Bup = B(2);
    RVis.drift = kvis(c) * sind(hdgs); 
    RVis.lose_flag = 1;
    RVis.plotflag = 0; % 1 = plot, 2 = plot and export_fig
    PVis{c} =  images_dtb_2d(RVis);
    
    wVes = sqrt( kves^2 / (kves^2 + kvis(c)^2) );
    wVis = sqrt( kvis(c)^2 / (kves^2 + kvis(c)^2) );
    
    for d = 1:length(deltas)
        % positive delta defined as ves to the left, vis to the right
        muVes = kves    * sind(hdgs-deltas(d)/2);
        muVis = kvis(c) * sind(hdgs+deltas(d)/2);
%         kcomb(c) = sqrt(kves^2 + kvis(c)^2);
%         R.drift = kcomb(c) * sind(ushdgs); 
        R.drift = wVes.*muVes + wVis.*muVis;

        R.t = 0.001:0.001:duration/1000;
        R.Bup = B(3);
        R.lose_flag = 1;
        R.plotflag = 0; % 1 = plot, 2 = plot and export_fig
        PComb{c}{d} =  images_dtb_2d(R);
    end
end

%%
RTdata = data.RT;
confdata = data.conf;

if options.conftask==1
    confdata = confdata*100;
end

usetrs_data  = data.correct | data.heading==0; % use only correct (and 0 hdg) trials for RT/SEP fits
choiceD  = logical(data.choice-1);
% choiceD  = data.correct;

pRight_model = nan(length(data.heading),1);
pHigh_model  = nan(length(data.heading),1);

n = nan(length(mods),length(cohs),length(deltas),length(hdgs));

pRight_fit = n;     pHigh_fit = n;
meanRT_fit = n;     meanRT_data = n;    sigmaRT = n;
meanConf_fit = n;   meanConf_data = n;  sigmaConf = n;
nCor = n;

% if options.dummyRun, keyboard, end

for m = 1:length(mods)
    
    if options.conftask>0
        if length(cL)>1, confLapse = cL(m);
        else confLapse = cL;
        end
    end
    
for c = 1:length(cohs)
for d = 1:length(deltas)

    for h = 1:length(hdgs)
%         uh = abs(hdgs(h))==ushdgs;
        
        Jdata  = data.modality==mods(m) & data.coherence==cohs(c) & data.heading==hdgs(h) & data.delta==deltas(d);
        
        if sum(Jdata)==0, continue, end
        
        n(m,c,d,h) = sum(Jdata);

        % pRight_model stores the predictions made by the model for the
        % probability of a rightward choice on each trial condition
        % and then assigned to the actual data trials of each condition in pRight_model for
        % log likelihood calculation
        
        if m==1,     Ptemp = PVes; muTnd = Tnds(m); logOddsMap = logOddsMap_All(:,:,m);
        elseif m==2, Ptemp = PVis{c}; muTnd = Tnds(m); logOddsMap = logOddsMap_All(:,:,m);
        elseif m==3, Ptemp = PComb{c}{d}; muTnd = Tnds(m); logOddsMap = logOddsMap_All(:,:,m);
        end
        
        % SJ 10-11-2021 replaced redundant Monte Carlo simulation with computations
        % from method of images 
        % for choices, P.up.p/P.lo.p contains probabilities of bound hit
        % for signed headings, up is rightward, lo is leftward

        pRight_model(Jdata) = Ptemp.up.p(h);
        pRight_fit(m,c,d,h) = Ptemp.up.p(h);
        
        % Palmer 2005 method
%         pRight_model(Jdata) = 1 ./ (1 + exp(-2*k*B*sind(hdgs(h))));
           
        % density of losing DV for this condition 
        Pxt = squeeze(Ptemp.up.distr_loser(h,:,:))' .* Ptemp.up.p(h);
        Pxt2= squeeze(Ptemp.lo.distr_loser(h,:,:))' .* Ptemp.lo.p(h);
%         Pxt = squeeze(Ptemp.up.distr_loser(h,:,:))';
%         Pxt2= squeeze(Ptemp.lo.distr_loser(h,:,:))';
        
        nCor(m,c,d,h) = sum(Jdata & usetrs_data);
        
        % CONF
        if options.conftask == 1 % sacc endpoint
%             meanConf_fit(m,c,d,h) = nansum(nansum(Pxt.*(2*logistic(logOddsMap)-1)));
            meanConf_fit(m,c,d,h) = (nansum(nansum(Pxt .* (2*logistic(logOddsMap)-1)) + ...
                nansum(Pxt2 .* (2*logistic(logOddsMap)-1))));
            
            % confidence lapse
            meanConf_fit(m,c,d,h) = (1-confLapse)*meanConf_fit(m,c,d,h) + confLapse*1/2; % random report
%             meanConf_fit(m,c,d,h) = meanConf_fit(m,c,d,h)+confLapse*(1-meanConf_fit(m,c,d,h)); % 'lapse' high conf
            
            meanConf_data(m,c,d,h) = mean(confdata(Jdata & usetrs_data));
            sigmaConf(m,c,d,h) = std(confdata(Jdata & usetrs_data)) / sqrt(nCor(m,c,d,h));
            
        elseif options.conftask == 2 % PDW
%             pHigh = sum(Pxt.*(logOddsMap>theta)); % sum of density above theta [but still is a function of time!]   
            pHigh = sum(Pxt.*(logOddsMap>theta)) + sum(Pxt2.*(logOddsMap>theta));
            pHigh = sum(pHigh);  % marginalize over time
            pHigh = pHigh+confLapse*(1-pHigh); % 'lapse' high bet
%             
            pHigh_model(Jdata) = pHigh; 
            pHigh_fit(m,c,d,h) = pHigh; % store by condition as well for ease of plotting later
        end
        
        % RT
        if options.RTtask           
%             cRT = Ptemp.up.mean_t(h);
            cRT = Ptemp.up.mean_t(h)*Ptemp.up.p(h) + Ptemp.lo.mean_t(h)*Ptemp.lo.p(h); % is this correct? 
            meanRT_fit(m,c,d,h) = cRT+muTnd;
            
            meanRT_data(m,c,d,h) = mean(RTdata(Jdata & usetrs_data));
            sigmaRT(m,c,d,h) = std(RTdata(Jdata & usetrs_data)) / sqrt(nCor(m,c,d,h));
            
%             meanRT_data(m,c,d,h) =  mean(RTdata(Jdata));
%             sigmaRT(m,c,d,h) = std(RTdata(Jdata)) / sqrt(n(m,c,d,h));


        end
        
    end
    
end
end
end



%% calculate log likelihoods

meanConf_fit = meanConf_fit * 100;
meanRT_fit = meanRT_fit * 1000;
meanRT_data = meanRT_data * 1000;


if options.dummyRun
   err = NaN;
   
else
      
% % likelihood of rightward choice on each trial, under binomial assumptions
% based on Palmer et al. 2005 coh version, but this ignores mod/coh in 3DMPtask
% Pr_model = 1 ./ (1 + exp(-2*k*B*sind(hdg)));
% kluge to avoid log(0) issues
pRight_model(pRight_model==0) = eps; 
pRight_model(pRight_model==1) = 1-eps;

% log likelihood of choice is summed log likelihood across all trials (log
% probability of observing choice given model)
LL_choice = nansum(log(pRight_model(choiceD))) + nansum(log(1-pRight_model(~choiceD)));

% log likelihood of mean RTs, under Gaussian approximation
LL_RT = 0;
if options.RTtask
    L_RT = 1./(sigmaRT.*sqrt(2*pi)) .* exp(-(meanRT_fit - meanRT_data).^2 ./ (2*sigmaRT.^2));
    L_RT(L_RT==0) = eps; %min(L_RT(L_RT~=0));
    LL_RT = nansum(log(L_RT(:))); % sum over all conditions
end

% likelihood of mean confidence (for SEP), or log probabilities for PDW
LL_conf = 0; % initialize, if no conf task will stay as 0
if options.conftask == 1 % SEP
    L_conf = 1./(sigmaConf*sqrt(2*pi)) .* exp(-(meanConf_fit - meanConf_data).^2 ./ (2*sigmaConf.^2));
    L_conf(L_conf==0) = eps; %min(L_conf(L_conf~=0));
    LL_conf = nansum(log(L_conf(:))); % sum over all conditions
    
elseif options.conftask == 2 % PDW
    pHigh_model(pHigh_model==0) = eps; 
    pHigh_model(pHigh_model==1) = 1-eps;
    PDW = logical(data.conf);
%     LL_conf = nansum(log(pHigh_model(PDW))) + nansum(log(1-pHigh_model(~PDW)));   
    LL_conf = nansum(log(pHigh_model(PDW & usetrs_data))) + nansum(log(1-pHigh_model(~PDW & usetrs_data)));

end

% LLs above are negative, and closer to 0 = higher log likelihood
% so sum the negatives here, then flip sign to get error in positive, and
% minimizing error
if options.whichFit == 3 || ~isfield(options,'whichFit')
%     err = -(LL_choice + LL_conf + LL_RT); % negate for minimization

    % LL_RT is order of magnitude smaller, so standardize to help joint minimization
    if options.conftask==2
        err = -(LL_choice/10 + LL_conf/10 + LL_RT); % negate for minimization
    else
        err = -(LL_choice + LL_conf + LL_RT); % negate for minimization
    end

elseif options.whichFit == 2
    err = -(LL_choice + LL_conf/10);
elseif options.whichFit == 1
     err = -(LL_choice + LL_RT);
elseif options.whichFit == 0
    err = -LL_choice;
end

end

%%

% assign same conditions for fit data, will be useful for plotting later!
fit.heading = data.heading;
fit.modality = data.modality;
fit.coherence = data.coherence;
fit.delta = data.delta;

fit.pRight = pRight_fit;
if options.conftask==1
    fit.confMean = meanConf_fit/100;
elseif options.conftask==2
    fit.confMean = pHigh_fit;
end
fit.RTmean = meanRT_fit/1000;


% copy vestib-only data to all coherences, to aid plotting
for c=1:length(cohs)
    fit.pRight(1,c,:,:)   = fit.pRight(1,1,:,:);
    fit.confMean(1,c,:,:) = fit.confMean(1,1,:,:);
    fit.RTmean(1,c,:,:)   = fit.RTmean(1,1,:,:);
end

%% print progress report!

if ~options.dummyRun
fprintf('\n\n\n****************************************\n');
fprintf('run %d\n', call_num);
fprintf('param values\n')
for i=1:length(param)
    fprintf('%s:\t%g\n',paramNames{i},param(i))
end
fprintf('\n');

fprintf('err:\t%.4f\n', err);
fprintf('choice:\t%.2f\nconf:\t%.2f\nRT:\t%.2f\n',LL_choice,LL_conf,LL_RT)


if options.ploterr && contains(options.fitMethod,'fms')
    figure(options.fh);
    set(gcf,'color','white');
    if call_num==1,clf; end
    %subplot(311); 
    hold on;
    xlabel('run'); ylabel('total error'); title('Log Likelihood Error');
    xlim([0 max(call_num,5)]);
    plot(call_num, err, '.','MarkerSize',16,'color','k');
    try tidyaxes; catch; end
    %{
    if call_num>1
        paramVals(:,call_num) = param;
        normX = norm(paramVals(:,call_num)-paramVals(:,call_num-1),Inf);
        tolX  = 0.1 * (1+norm(paramVals(:,call_num-1),Inf));
        
        errVals(call_num) = err;
        normF = norm(errVals(call_num)-errVals(call_num-1),Inf);
        tolF  = 0.1 * (1+errVals(call_num-1));

        subplot(312); hold on; xlabel('run'); ylabel('Xdiff'); 
        plot(call_num, normX, '.','MarkerSize',16,'color','k');
        plot(call_num, tolX, 'x','MarkerSize',16,'color','r');
        xlim([0 max(call_num,5)]);

        subplot(313); hold on; xlabel('run'); ylabel('Fdiff'); 
        plot(call_num, normF, '.','MarkerSize',16,'color','k');
        plot(call_num, tolF, 'x','MarkerSize',16,'color','r');

        xlim([0 max(call_num,5)]);
    end
    %}
    
    drawnow;
end
call_num = call_num + 1;


end

end


% retrieve the full parameter set given the adjustable and fixed parameters 
function param2 = getParam ( param1 , guess , fixed )
  
  if all(fixed), param2 = param1; return; end

  param2(fixed==0) = param1;            %get adjustable parameters from param1
  param2(fixed==1) = guess(fixed==1);   %get fixed parameters from guess (initial point)
end
