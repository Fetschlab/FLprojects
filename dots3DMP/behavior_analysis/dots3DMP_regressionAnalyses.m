% function dots3DMP_regressionAnalyses(data,mods,cohs,deltas,conftask,RTtask)

% SJ 06/2020 started
% SJ 10/2021 modified significantly
% SJ 04/2022 modified again

% Choice/Correct data is going to be fitted by logistic regression (MLE method
% under binomial assumptions i.e. trials are a list of Bernoulli coin tosses)
% Same with PDW
% SEP will require linear reg

% Regression models

% 1.  effect of combined condition on accuracy
% 1b. effect of heading strength on accuracy
% 2. effect of heading strength and RT on high bet probability/SEP
% (logistic for PDW, multiple linear reg for SEP)
% 3. effect of heading strength and RT on probability correct (logistic)
% 4. improvement in accuracy when high bet is chosen


subjs = unique(data.subj);

% general tstatistic function for logistic regression to calculate p-value
tstatfun = @(beta,se) 1-tcdf(beta./se,length(beta)-1);

mods = unique(data.modality);
cohs = unique(data.coherence);
deltas = unique(data.delta);

%% 1. improvement in accuracy with combined condition (non-conflict only)

% can use choice and signed heading or correct and unsigned heading?
J = data.delta==0;
I = data.modality(J)==3;
hdg = data.heading(J);

D = [ones(size(hdg)), I, hdg.*I, hdg, data.choice(J)-1];
[beta_AccMod,llik,pred,se] = logistfit_se(D);
p_AccMod = tstatfun(beta_AccMod,se);

outputStats.AccMod.beta = beta_AccMod;
outputStats.AccMod.p    = p_AccMod;

%% 1b. effect of heading strength on P(correct)

% isn't this just a subset of the above?
% J = data.delta==0;
% D = [ones(size(data.heading(J))), data.heading(J), data.choice(J)-1];
% [beta_AccHdg,llik,pred,se] = logistfit_se(D);
% p_AccHdg = tstatfun(beta_AccHdg,se);
% 
% outputStats.AccHdg.beta = beta_AccHdg;
% outputStats.AccHdg.p    = p_AccHdg;

%% relationship of RT and choice certainty
% see Figure 2, equation 4 Kiani 2014

clear lm*
for s = 1:length(subjs)+1
    temp = data;
    
    if s<length(subjs)+1
        fnames = fieldnames(data);
        for f=1:length(fnames)
            temp.(fnames{f})(~strcmp(data.subj,subjs{s})) = [];
        end
    end
    
    if conftask==1
        tbl = table(abs(temp.heading),temp.RT,temp.conf,temp.modality,'VariableNames',{'Hdg','RT','Conf','Mod'});
        lm1 = fitlm(tbl,'interactions','ResponseVar','Conf','PredictorVars',{'Hdg','RT','Mod'},'CategoricalVar','Mod');
        lmHdgAll(s) = lm1.Coefficients.pValue(3);

        for m=1:length(mods)
            for c=1:length(cohs)
                if m==1 && c>1, continue,end
                J = temp.modality==mods(m) & temp.coherence==cohs(c) & temp.delta==0;
                tbl = table(abs(temp.heading(J)),temp.RT(J),temp.conf(J),'VariableNames',{'Hdg','RT','Conf'});
                lm2 = fitlm(tbl,'interactions','ResponseVar','Conf','PredictorVars',{'Hdg','RT'});
                lmHdgModCoh(m,c,s) = lm2.Coefficients.pValue(3);
            end
        end
    elseif conftask==2
        % P(High Bet) = [1 + exp(-b0 + b1*hdg + b2*RT) ] ^ -1
        
        D = [ones(size(data.heading(J))), abs(data.heading(J)), data.PDW(J)];
        [beta_ConfHdg,llik,pred,se] = logistfit_se(D);
        
        pConf = tstatfun(beta_ConfHdg,se);
        
    end
end

    

%% 3b. effect of confidence on P(correct)

if conftask==1
elseif conftask==2
    D = [ones(size(data.heading)), data.PDW, abs(data.heading), data.correct];
    [beta_AccConf,llik,pred,se] = logistfit_se(D);
    p_AccConf = tstatfun(beta_AccConf,se);
end

%% 5.are confidence judgements associated with different reaction times

% separately for each mod and coh!


p = nan(length(mods),length(cohs),3);
for m=1:length(mods)
    for c=1:length(cohs)      
        J = data.modality==mods(m) & data.coherence == cohs(c) & data.delta == 0;
        
        if sum(J)==0, continue, end
        [p(m,c,:),t,stats,terms] = anovan(data.RT(J),[abs(data.heading(J)) data.PDW(J)],'model','interaction','display','off');
    end
end


%%

%{



%%
% 1. quantify change in sensitivity associated with combined condition

I = data.modality==3;

% Pright = 1 + exp(-(b0 + b1*I + b2*I*hdg + b3*hdg))
D = [ones(size(data.heading)), I, I.*data.heading, data.heading data.choice-1];
%                                                               ^ subtract 1 to make choices 0 or 1                                                                

[beta,llik,pred,se] = logistfit_se(D);

% null hypothesis is that combined condition has no effect, i.e. b2 = 0
% tstat = (m-mu)/se  (mu=0 so m-mu is just beta)
p = tstatfun(beta,se);

% reject null if p(3) < 0.05 --> combined condition interacts with hdg to
% affect choice, and beta coefficient is positive, so combined condition significantly increases sensitivity  
% compare results when indicator variable refers to ves or vis only!

%% RT and conf (c.f. Kiani 2009/2014 analyses)

% 2. fit relationship between heading, RT and conf, separately for correct and
% error trials

goodRT = true(size(data.RT));

% goodRT = data.RT~=mode(data.RT); % kluge, in hand-tuned sim there is a risk that a
% lot of RTs will max out at stimulus time, so that too many weak stim (and
% error) trials have the same fixed RT - exclude these for now

% normalize, esp if there are different ranges across sessions/subjects
%RT = data.RT/max(data.RT); 
RT = data.RT;

% excluding zero heading trials, if any
data.correct = (data.choice==2 & data.heading>0) | (data.choice==1 & data.heading<0);
cor = data.correct & goodRT; 
err = ~data.correct & goodRT;

% correct trials
S = data.conf(cor);
% D = [ones(size(data.heading(cor))), abs(data.heading(cor)), RT(cor)];
% [Bcor,Bint_cor,~,~,statscor] = regress(S,D);

% use fitlm instead of regress (Bcor should match Estimate from fitlm)
tbl = table(abs(data.heading(cor)),RT(cor),S,'VariableNames',{'Hdg','RT','Conf'});
lmCor = fitlm(tbl,'Conf~Hdg+RT');

% repeat for error trials
S = data.conf(err);
% D = [ones(size(data.heading(err))), abs(data.heading(err)), RT(err)];
% [Berr,Bint_err] = regress(S,D);

tbl = table(abs(data.heading(err)),RT(err),S,'VariableNames',{'Hdg','RT','Conf'});
lmErr = fitlm(tbl,'Conf~Hdg+RT');

% subjects are more confident for faster decisions, both correct and
% incorrect


%% 3. does slope of regression change for correct and error trials?

% null hypothesis, b5 = 0
%I = abs(data.heading)<=2.5; % use only low hdg trials
I = abs(data.heading)>0.01; % use all trials
I = I & goodRT;
S = data.conf(I); 

% D = [ones(size(data.heading(I))), abs(data.heading(I)), RT(I),...
%         err(I), err(I).*abs(data.heading(I)), err(I).*RT(I)];
% [Bcorerr,Bint_corerr] = regress(S,D);

% or use fitlm
tbl = table(abs(data.heading(I)),RT(I),err(I),S,'VariableNames',{'Hdg','RT','Corr','Conf'});
lmCorErr = fitlm(tbl,'Conf~Hdg+RT+Corr+Corr*Hdg+Corr*RT');
 
% no - RT:Corr interaction coefficient is not different from 0

%% 4. relationship between hdg and certainty for errors
% I = abs(data.heading)<=6; % use only low hdg trials (not enough after RT exclusion!)
S = data.conf(err & I); 

% D = [ones(size(data.heading(err & I))), abs(data.heading(err & I))];
% [Bhdgerr,Bint_hdgerr] = regress(S,D);

tbl = table(abs(data.heading(err&I)),S,'VariableNames',{'Hdg','Conf'});
lmHdgErr = fitlm(tbl,'Conf~Hdg');

% significant +ve coefficient - on error trials, subjects are more certain
% for larger heading directions


%% 5. effect of RT and hdg on P(correct)

% Pcor= 1 / (1 + exp(-(b0 + b1*RT + b2*hdg)))
D = [ones(size(data.heading(goodRT))), RT(goodRT), abs(data.heading(goodRT)), data.correct(goodRT)];
[beta,llik,pred,se] = logistfit_se(D);

% null hypothesis is that RT has no effect on Pcor , i.e. b1 = 0
p = 1-tcdf(beta./se,length(beta)-1);


%}