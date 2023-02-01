function [X, err_final, fit, parsedFit, fitInterp] = dots3DMP_fitDDM(data,options,guess,fixed)

% parameter bounds for fitting

try
    LB = options.lowerbound;
    UB = options.upperbound;
catch
    LB = guess/5;
    UB = guess*5;
end

% also set the "plausible" lower/upper bounds used by BADS
PLB = guess/2;
PUB = guess*2;

global call_num; call_num=1;

global paramVals; paramVals = guess';
global errVals; errVals = NaN;

options.dummyRun = 0;

if all(fixed)
    X = guess;
else
    
    switch options.fitMethod
        case 'fms'
            fitOptions = optimset('Display', 'iter', 'MaxFunEvals', 100*sum(fixed==0), 'MaxIter', ... 
                100*sum(fixed==0), 'TolX',1e-1,'TolFun',1e-1,'UseParallel','Always');
            [X, fval, ~] = fminsearch(@(x) feval(options.errfcn,x,guess,fixed,data,options), guess(fixed==0), fitOptions);
    %         [X, fval, exitflag] = fminunc(@(x) feval(options.errfcn,x,guess,fixed,data,options), guess(fixed==0), fitOptions);
    %         fprintf('fval: %f\n', fval);

        case 'fmsbnd'
            fitOptions = optimset('Display', 'iter', 'MaxFunEvals', 100*sum(fixed==0), 'MaxIter', ... 
                100*sum(fixed==0), 'TolX',1e-3,'TolFun',1e-3,'UseParallel','Always');
            [X, fval, ~] = fminsearchbnd(@(x) feval(options.errfcn,x,guess,fixed,data,options), guess(fixed==0), LB(fixed==0), UB(fixed==0), fitOptions);

        case 'global'
            % GlobalSearch from Global Optimization Toolbox
            fitOptions = optimoptions(@fmincon,'Display','iter',...
                'Algorithm','interior-point',...
                'FinDiffType','central',...
                'FinDiffRelStep',1e-4,...
                'UseParallel','always');
            problem = createOptimProblem('fmincon','x0',guess(fixed==0),'objective',...
                @(x) feval(options.errfcn,x,guess,fixed,data,options),'lb',LB(fixed==0),...
                'ub',UB(fixed==0),'options',fitOptions);
            gs = GlobalSearch;
            [X,~,~,~,~] = run(gs,problem);      

        case 'multi'
            % MultiStart from Global Optimization Toolbox
            fitOptions = optimoptions(@fmincon,'Display','iter',...
                'Algorithm','interior-point',...
                'FinDiffType','central',...
                'FinDiffRelStep',1e-1,...
                'UseParallel','always');
            problem = createOptimProblem('fmincon','x0',guess(fixed==0),'objective',...
                @(x) feval(options.errfcn,x,guess,fixed,data,options),'lb',LB(fixed==0),...
                'ub',UB(fixed==0),'options',fitOptions);
            ms = MultiStart('FunctionTolerance',1e3,'XTolerance',1,...
                'StartPointsToRun','bounds-ineqs','UseParallel','always');
            [X,~,~,~,~] = run(ms,problem,200);    

        case 'pattern'
            % Pattern Search from Global Optimization Toolbox
            fitOptions = psoptimset(@patternsearch);
            fitOptions = psoptimset(fitOptions,'CompletePoll','on','Vectorized','off',...
                'UseParallel','always');
    %         [X,~,~,~] = patternsearch(@(x) dots3DMP_fitDDM_err(x,data,options),...
    %             guess(fixed==0),[],[],[],[],LB(fixed==0),UB(fixed==0),[],fitOptions);

            [X,~,~,~] = patternsearch(@(x) feval(options.errfcn,x,guess,fixed,data,options),...
                guess(fixed==0),[],[],[],[],LB(fixed==0),UB(fixed==0),[],fitOptions);

        case 'bads'

            % We define the negative log-likelihood function via a call to IBSLIKE
            % (IBSLIKE provides a vectorized implementation of IBS for MATLAB; check
            % out IBS_BASIC for a bare bones implementation)

            % Here's how this works, with reference to Luigi's ibs_example.m: 
            % R = psycho_gen(theta,S) is his simulator, taking S (a vector of
            % orientations) and generating responses (-1 or 1)

            % Then he writes the function handle to pass into bads:
            % nllfun_ibs = @(theta) ibslike(@psycho_gen,theta,R,S,options_ibs);
            % where ibslike takes arguments: 
            % (fun,params,respMat,designMat,options,varargin)

            % In the example, designMat is the same as S, because the stimuli
            % were generated randomly from a continuous distribution and not
            % repeated (recall that IBS only requires responses to be discrete,
            % not stimuli) 

            % This made me think that designMat should not be the entire matrix
            % of trials, including repeats, but only the unique trial types...
            % Now I don't think that's right; it really should include all trials.

            % Another confusing thing is that ibslike says designMat is
            % optional, and if omitted, it works with simulated data... by 
            % simply passing trial indices (1,2,3...) into fun?? That's weird!
            % How would fun be expected to generate reasonable responses with
            % trial indices as input, unless you wrote in a way for it to do so?
            % (his psycho_gen doesn't do this)

            % Anyway, what we need is to turn data into two matrices, R and S:
            S = [data.modality data.heading data.coherence data.delta];

            % also need to discretize RT and conf
            nbins = 10;
            RTedges = [0 quantile(data.RT,nbins-1) inf];
            RTbinned = discretize(data.RT,RTedges);

            confEdges = [-inf quantile(data.conf,nbins-1) inf]; % no need to assume 0..1
            confBinned = discretize(data.conf,confEdges);

    %         % temp: debugging
    %         [R,~] = dots3DMP_simDDM_forIBS(guess,S,'RTedges',RTedges,'confEdges',confEdges);

            R = [data.choice RTbinned confBinned];

            options_ibs.Nreps = 10;

            % temp: debug with ibs_basic
    %         ibs_basic(@dots3DMP_simDDM_forIBS,guess,R,S,'RTedges',RTedges,'confEdges',confEdges);

            % temp: running ibslike even once is prohibitively slow!
            ibslike(@dots3DMP_simDDM_forIBS,guess,R,S,options_ibs,'RTedges',RTedges,'confEdges',confEdges)
            % (ibslike takes varargin to pass other variables to fun beyond
            % just param (theta) and S, which we need to do for the bin edges)
            keyboard


            % if it did work, we'd proceed with:

            nllfun_ibs = @(x) ibslike(@dots3DMP_simDDM_forIBS,x,R,S,options_ibs,'RTedges',RTedges,'confEdges',confEdges);

            % As a starting point for the optimization, we draw a sample inside
            % the plausible box (in practice you should use multiple restarts!)
            guess = rand(size(LB)).*(PUB-PLB) + PLB;

            X = bads(nllfun_ibs,guess,LB,UB,PLB,PUB); % 

    end

    % re-insert initial fixed params in the appropriate places for final error
    % calcs (since we have only fit the non-fixed ones).
    try
        temp = X;
        X = zeros(size(fixed));
        X(fixed==0) = temp;
        X(fixed==1) = guess(fixed==1);
    catch
        keyboard
    end

end
    
%%
% run err func again at the fitted/fixed params to generate a final
% error value and model-generated data points (trial outcomes)

options.ploterr = 0;

options.whichFit  = {'multinom','RT'}; % choice, conf, RT, multinom (choice+conf)
[err_final, fit, parsedFit] = feval(options.errfcn,X,X,true(size(X)),data,options);


% THIS SHOULD ALL BECOME OBSOLETE, NO MORE MC!
if options.runInterpFit

    % *** now run it one more time with interpolated headings
    % to generate smooth plots ***
        % (this ordinarily works well to show the fits vs. data, but in the
        % current version the model fits/predictions are generated by monte
        % carlo simulation and are thus too noisy to be worth repeating at
        % finely spaced headings (you can try uncommenting the next two lines
        % to see what this looks like). So for now hdgs will be the actual hdgs\
        % from data and this section will seem wholly redundant, but will be
        % useful later)

nsteps = 33; % should be odd so there's a zero
hdgs = linspace(min(data.heading),max(data.heading),nsteps);
% hdgs(hdgs==0) = [];
% hdgs = unique(data.heading)'; % TEMP, see above comment
%                                   ^ no longer an issue, we just fit the fits with (c)gauss


% Dfit is a dummy dataset with the same proportions of all trial types,
% just repeated for the new interpolated headings (and some arbitrary nreps)
mods = unique(data.modality);
cohs = unique(data.coherence);
deltas = unique(data.delta);

% nreps = 200;
nreps = 1;
[Dfit.heading, Dfit.modality, Dfit.coherence, Dfit.delta, ~] = dots3DMP_create_trial_list(hdgs,mods,cohs,deltas,nreps,0);

    % the observables are just placeholders because the err func still needs to
    % calculate err even though in this case it's meaningless. What will be
    % plotted is fitInterp, not Dfit

    % in other words, we are using the already calculated fit parameters to
    % plot curves for simulated trial conditions

Dfit.choice  = ones(size(Dfit.heading));
Dfit.RT      = ones(size(Dfit.heading));
Dfit.conf    = ones(size(Dfit.heading));
Dfit.correct = ones(size(Dfit.heading));

if options.conftask==2
    Dfit.PDW = ones(size(Dfit.heading));
end

% [~,fitInterp] = dots3DMP_fitDDM_err(X,Dfit);
fixed = true(size(X)); % fix all params
options.dummyRun = 1;
[~, fitInterp] = feval(options.errfcn,X,X,fixed,Dfit,options);

else
    fitInterp = fit;
end



