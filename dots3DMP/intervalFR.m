function [meanFR,semFR,trialFR,tInt,groups,groupIC] = intervalFR(spiketimes,ev1,ev2,groupVar,resort)
%
% calculate firing rate within specified intervals for given set of
% spiketimes, over indep variable and grouping variable (indepVar will be
% the first column in groupVar)
%
% groupVar should be an array with each column containing trial condition to
% group over (one column per indepVar) 
% the final column should be the variable to plot on x-axis! e.g. heading
%
% ev1 should contain the starts of the intervals (one column per interval)
% ev2 should contain the ends, 1 entry per trial
%
% e.g. ev1 = [fixation stimOn stimOff];
%      ev2 = [targsOn  stimOff reward];
%
% columns can also be specified relative to an event e.g. stimOn+2
%
% TO DO clearer documentation on inputs + outputs and usage

if nargin<5,resort=0;end

[ntrs,nInt] = size(ev1);
if nargin<4 || isempty(groupVar)
    groupVar = ones(1,ntrs);
end
[groups,groupIA,groupIC] = unique(groupVar,'rows');

ntrsPerCond = hist(groupIC,unique(groupIC));

[trialFR,tInt] = deal(nan(ntrs,nInt));
[meanFR,semFR] = deal(nan(size(groups,1),nInt));

% calculate firing rate for each trial and interval (trialFR), then calculate mean
% and sem for each condition group
for ii=1:nInt
   
    tInt(:,ii) = ev2(:,ii) - ev1(:,ii);
        
    xgt = bsxfun(@gt,spiketimes,ev1(:,ii)');
    xlt = bsxfun(@lt,spiketimes,ev2(:,ii)');
    
    spkCount = sum(xgt & xlt,1);
 
    trialFR(:,ii) = spkCount' ./ tInt(:,ii);
    
    meanFR(:,ii) = accumarray(groupIC,trialFR(:,ii),size(groupIA),@nanmean);
    semFR(:,ii)  = accumarray(groupIC,trialFR(:,ii),size(groupIA),@nanstd) ./ sqrt(ntrsPerCond');
end

% reshape so that 'heading' or main independent variable is separate axis
% to other conditions, makes plotting easier in some cases
if resort
uhdgs=unique(groupVar(:,end));
meanFR = reshape(meanFR,length(uhdgs),[],size(meanFR,2));
semFR  = reshape(semFR,length(uhdgs),[],size(semFR,2));
groups = reshape(groups,length(uhdgs),[],size(groups,2));
end
