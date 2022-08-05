% scriptified the post-processing after simDDM (1D and 2D), to standardize


% adjust wager probability for base rate of low bets, as seen in data
% ('compresses' the curve, not an offset, because P(high) varies with coh)
pdw_preAlpha = pdw; % save original PDW for plotting choice/RT splits (relationship between conf and choice/RT is, by construction, independent of alpha)
pdw(pdw==1 & rand(length(pdw),1)<alpha) = 0;

% add non-decision time (truncated normal dist)
Tnd = zeros(ntrials,1);
for n = 1:ntrials
    while Tnd(n)<=TndMin || Tnd(n)>=TndMax % simple trick for truncating, requires integers (ms)
        Tnd(n) = round(normrnd(TndMean,TndSD));
    end
end
DT = RT; % rename this 'decision time'
RT = DT+Tnd;

% quick sanity check that params are reasonable
pCorrect_total = sum(sign(choice)==sign(coh)) / ntrials

% should be >0.95 or else maxDur isn't long enough (or need urgency!)
hitBoundPct = sum(hitBound)/length(hitBound)

%% format data as in experimental data files and generate output structs

coh(coh==0) = sign(randn)*eps; % should have no actual zeros, but if so, sign them randomly;
                               % this is just to assign a direction and correct/error
data.correct = choice==sign(coh);
data.direction = nan(ntrials,1);
data.direction(coh>0) = 0;
data.direction(coh<0) = 180;
% coh(abs(coh)<1e-6) = 0; % now go back to one 'zero' [OR NOT!]
data.coherence = abs(coh);
data.scoh = coh;

data.choice = choice;
data.choice(data.choice==-1) = 0; % code elsewhere assumes 0s and 1s
data.RT = RT/1000; % convert to seconds
data.PDW = pdw;
data.PDW_preAlpha = pdw_preAlpha;
data.conf = conf;
data.dur = ones(size(data.coherence))*max_dur;

conftask = 2; % pdw (2) for now, even though conf rating can be generated here
RTtask = 1; RTCorrOnly = 0;
parsedData = Dots_parseData(data,conftask,RTtask,RTCorrOnly);

% plot
cohs = unique(coh); wFit = 0; forTalk = 0;
Dots_plot(parsedData,cohs,conftask,RTtask,wFit,forTalk)


%% temp: save data, e.g. for param recovery

save tempsim.mat data origParams

% if needed for debugging fitting code:
% Psim = P; Rsim = R;
% save tempdata.mat Psim Rsim
