% Miguel Vivar-Lazo
% 08/16/2021
% Fetsch Lab: Average Values Conditioned On Something
function [avgDepVar, N, stderr, uniqIndVar] = behavioralAverages(depVariable, indVariable)
uniqIndVar = unique(indVariable);
avgDepVar = nan(1,length(uniqIndVar));
N = nan(1,length(uniqIndVar));
stderr = nan(1,length(uniqIndVar));
for i = 1:length(uniqIndVar)
avgDepVar(i) = mean(depVariable(uniqIndVar(i) == indVariable), 'omitnan');
N(i) = nansum(uniqIndVar(i) == indVariable);
stderr(i) = std(depVariable(uniqIndVar(i) == indVariable), 'omitnan')/sqrt(N(i));
end

