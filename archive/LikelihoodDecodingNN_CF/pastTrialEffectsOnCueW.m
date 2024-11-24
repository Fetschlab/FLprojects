function [vesW vesW_avg] = pastTrialEffectsOnCueW(tr,data,LI,trPerDelta,boostraps) % xAxisPos,clr --> seperate plotting function?
% xAxis and clr should be same length as size(LI,2)
%   if not, either throw error AT START of function, or have default settings for both vars

data.delta(data.delta==4) = 3;
data.delta(data.delta==-4) = -3;
mod = data.modality(tr);
coh = data.coherence(tr);
delta = data.delta(tr);
% conf = data.PDW(tr);

mods = unique(data.modality);
hdngs = unique(data.heading);
cohs = unique(data.coherence);
deltas = unique(data.delta);

% trTypes = [];  % trTypes(:,I) = [trTypes(:,I) randsample(LI(I) & deltas(d),trPerDelta)];

trTypes = zeros(trPerDelta*length(deltas),size(LI,2));
for b = 1:boostraps
    for I = 1:size(LI,2)
        for d = 1:length(deltas)
            trTypes( (d-1)*trPerDelta+1:(d-1)*trPerDelta+trPerDelta, I) = randsample(find(delta(LI(:,I)) == deltas(d)),trPerDelta,true); % Create indices to trials of interest
        end
    end

    fn = fieldnames(data);
    for da = 1:size(trTypes,2) % for each data struct,
        for f = 1:length(fn)    % fill in each field
             datas(da).(fn{f}) = data.(fn{f})(trTypes(:,da)); % Build one struct at a time, 'da', and one field at a time,fn{f}, results in 'datas' as an array of parallel structs. Is there a way to build all structs simultaneously? Such that one data.field(LI), operates column by column and stores appropriately each column in diff struct ...
        end
    end

    for da = 1:size(trTypes,2)
        gfit_datas(da) = dots3DMP_fit_cgauss_NN(datas(da),mods,cohs,deltas);
        vesW{da}(b,:) = dots3DMP_wgts_thres_NN(gfit_datas(da).muPMF,gfit_datas(da).sigmaPMF,cohs,deltas); % Each cell is one trial type (one column of LI), each row in each cell is one bootstrap itiration 
    end
end

vesW_avg = cellfun(@(x) mean(x, 1), vesW, 'UniformOutput', false); % Takes average across rows of each cell in vesWkkkkkkkkkkkkkkkkkkkk and stores each in a cell

end
    
% gfit_datas = arrayfun(@(da) dots3DMP_fit_cgauss_NN(datas(da), mods,cohs, deltas), 1:size(trTypes, 2));
% vesW = arrayfun(@(da) dots3DMP_wgts_thres_NN(gfit_datas(da).muPMF, gfit_datas(da).sigmaPMF, cohs, deltas), 1:size(trTypes, 2), 'UniformOutput', false);


% Testing cellFun
% xx = {[2,3], [4,5], [6,7]};
% averageCell = cellfun(@(x) mean(x, 2), xx, 'UniformOutput', false)

    % ChatGPT assisting in removing loops and linearizing operations...
    % w.i.p.
% %     gfit_datas = arrayfun(@(da) dots3DMP_fit_cgauss_NN(datas(da), mods,cohs, deltas), 1:size(trTypes, 2)); % 'datas' is an array of structs. Not sure if dots3DMP function can handle having all those structs accessed simultaneously
% %     vesW = arrayfun(@(da) dots3DMP_wgts_thres_NN(gfit_datas(da).muPMF, gfit_datas(da).sigmaPMF, cohs, deltas), 1:size(trTypes, 2), 'UniformOutput', false);
% 2 inputs into arrayfun, first is a function, second is the elements to
% which the specified function should be applied to. Here we are using @(a)
% to specify a function that indicates operations should occur in the
% indices specified by 'a'. Using arrayfun, the second input now gets input
% into 'a'. Thus, the indices 1:size(trTypes,2) replace 'a' in the function dots3DMP!
% gfits_datas output will be of dimensionality defined by 1:size(trTypes,2)
% Looks good, only worry is whether 3DMP functions can deal with having
% multiple structs indexed simultaneously, will it know to take it one by
% one, applying elementwise?
    




% 
% 
% 
% % Vectorized alternative
% % Simply make tr the full vector and apply all at once, no need to loop
% % through individually... who knew... ty chatGPT
%     % Past trial vesO + Hi Conf vs. Lo Conf
% tr = 3:length(coh);
% vesOnlyX1_HighConf = find(mod(tr-1)==mods(2) & conf(tr-1)==1) +2; % +2 to align index to original vectors
% vesOnlyX1_LoConf = find(mod(tr-1)==mods(1) & conf(tr-1)==0) +2;
% precX2_cntrlMod1 = find(mod(tr-1)~=mod(tr-2) & coh(tr-1)~=coh(tr-2)) +2;
% 
% 
% 
% 
% 
% 
% 
% 
% % Another alternative is staggering vectors so position 1 aligns with pos 2
% % and compare values stored by consecutive positions in a vector.
% % 1) mod on past trial is mod1 & past trial conf is Hi
% mod1nHiConf = mod(1:end-1)==1 & conf(1:end-1)==1 & mod(2:end)==3;
% 
% 
% %2) Past 2 trials are mod1
% % mod(1:end-2) == mod(2:end-1); --> outputs binary variable of length(mod)
% % -1. Use this binary vector as index into mod(3:end) to indicate trials
% % preceded by mod1 x2
% 
% 
% 
% 
% 
% 
% 
% 
% Inputs: 
% data_file, file_location, LIs, trPerDelta, bootstrapNum, clr{}
% 
% %plots L --> R along LI columns
% 
% 
% if LI
%     cond1 = true;
% else false
% end
% if  LI2
%     cond2 = true;
% else
%     cond2=false;
% end
% 
% trPerDelta
% 
% 
% cond = nan(size(LI);
% 
% LI(tr)