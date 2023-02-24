% 

% Not all of the days in this dataset have the cue-conflict condition, but
% keeping them in because we can still look at choice biases and thresholds
% over days

clear; clc

cd /Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/dataStructs/
load lucio_20211101-20220809_partClean.mat

% partClean because I have not re-assigned cohs and headings etc, so that I
% can go through each day individually and make notes on performance


alldata = data; % re-cast as alldata, so that we can store individual blocks in data

conftask = 2; % 1=colorbars, 2=PDW
RTtask = 1;

%%

mods   = unique(data.modality); 
cohs   = unique(data.coherence); 
deltas = unique(data.delta);
hdgs   = unique(data.heading);

days = unique(data.date); %date
nDays = length(days);

files = unique(data.filename); %  blocks
nFiles = length(files);

%%

fnames = fieldnames(alldata);

% loop over days

clear nCorrect* nTrials*
clear parsedData gfit wves deltas_day
for d = 1:nDays  %1:nFiles
    
    clear data;
    trials_day = alldata.date == days(d);
%     trials_day = strcmp(alldata.filename,files{d});

    % cutoff for the min number of trials to get reliable estimates
    
%     if sum(trials_day) > 600
        
        for F = 1:length(fnames)
            data.(fnames{F}) = alldata.(fnames{F})(trials_day);
        end
        nTrials(d,1) = sum(trials_day);
        nCorrect(d,1) = sum(data.correct);
        
        nTrialsConf(d,1) = sum(data.PDW==1 & ~data.oneTargConf);
        nTrialsConf(d,2) = sum(data.PDW==0 & ~data.oneTargConf);
        nTrialsConf(d,3)  = sum(data.oneTargConf);
        
        nCorrectConf(d,1) = sum(data.correct & data.PDW==1 & ~data.oneTargConf);
        nCorrectConf(d,2) = sum(data.correct & data.PDW==0 & ~data.oneTargConf);
        nCorrectConf(d,3) = sum(data.correct & data.oneTargConf);
        
        deltas_day{d,1} = unique(data.delta);
        hasDelta(d) = any(data.delta~=0);
        
        % now remove all the one-targ conf trials for parsing the gross
        % behavioral metrics
        oneTargConf = logical(data.oneTargConf);
        for F = 1:length(fnames)
            data.(fnames{F})(oneTargConf) = [];
        end
        
        mods_day{d} = unique(data.modality);
        cohs_day{d} = unique(data.coherence);
        hdgs_day{d} = unique(data.heading);
        
        parsedData(d) = dots3DMP_parseData(data,mods_day{d},cohs_day{d},deltas_day{d},hdgs_day{d},conftask,RTtask); % basic parsing of data, logistic fits
        
        gfit(d) = dots3DMP_fit_cgauss(data,mods,cohs,deltas_day{d},conftask,RTtask); % gaussian fitting
%         wves(d) = dots3DMP_cueWeights(gfit(d),cohs,deltas_day{d});
       
            
%     end

end

%% plot each day, for visually assessment

% d = 11;
% mods_day{d}

for d = 20:66
dots3DMP_plots(parsedData(d),mods_day{d},cohs_day{d},deltas_day{d},hdgs_day{d},conftask,RTtask)
fprintf('%d Total number of trials: %d\n', days(d), sum(parsedData(d).n(:)))
pause
clf
end

%% plot example day

day = 3; % 
dots3DMP_plots_cgauss_func(gfit(day),parsedData(day),mods,cohs,deltas_day{day},hdgs,conftask,RTtask)

%% summary statistics over days

% loop over days to extract necessary values into a simpler format for
% plotting

accuracy = nCorrect./nTrials;
accuracyConf = nCorrectConf./nTrialsConf;

[wvesPred,wvesEmp] = deal(nan(nDays,length(cohs)));

[choiceBias,confBias,RTBias,choiceSigma,confSigma,RTSigma] = ...
    deal(nan(length(mods),length(cohs),length(deltas)+1,nDays));

clear *Bias

for d=1:nDays
    
    % thresholds

    % biases/shifts
    if length(deltas_day{d})>1
        choiceBias(:,:,1:length(deltas)+1,d) = gfit(d).choice.mu;
        confBias(:,:,1:length(deltas)+1,d)   = gfit(d).conf.mu;
        RTBias(:,:,1:length(deltas)+1,d)     = gfit(d).RT.mu;
        
        choiceSigma(:,:,1:length(deltas)+1,d) = gfit(d).choice.sigma;
        confSigma(:,:,1:length(deltas)+1,d)   = gfit(d).conf.sigma;
        RTSigma(:,:,1:length(deltas)+1,d)     = gfit(d).RT.sigma;
        
    else
        choiceBias(:,:,length(deltas)+1,d)   = gfit(d).choice.mu(:,:,end);
        confBias(:,:,length(deltas)+1,d)   = gfit(d).conf.mu(:,:,end);
        RTBias(:,:,length(deltas)+1,d)   = gfit(d).RT.mu(:,:,end);
        
        choiceSigma(:,:,length(deltas)+1,d) = gfit(d).choice.sigma(:,:,end);
        confSigma(:,:,length(deltas)+1,d)   = gfit(d).conf.sigma(:,:,end);
        RTSigma(:,:,length(deltas)+1,d)     = gfit(d).RT.sigma(:,:,end);
        
    end
    
    % wVes
    wvesPred(d,:) = wves(d).choice.pred;

    if length(deltas_day{d})>1
        wvesEmp(d,:)  = wves(d).choice.emp;

    end    
end


%%

% ACCURACY
figure; subplot(211); hold on
plot(1:nDays,accuracy,'k',1:nDays,accuracyConf(:,1),'r',1:nDays,accuracyConf(:,2),'r--')
legend('Mean acc','High Bet Acc','Low Bet Acc')
xlabel('Session no.'); ylabel('Accuracy');
subplot(212); hold on
hr=refline(1,0); hr.Color = 'r';
scatter(accuracyConf(:,2),accuracyConf(:,1));
axis([0.5 0.9 0.5 0.9])
xlabel('Low Bet Acc'); ylabel('High Bet Acc');
%%
% VESTIBULAR WEIGHT
figure; hold on;
plot(1:nDays,wvesPred(:,1),'k--',1:nDays,wvesPred(:,2),'k-','linew',2)
plot(1:nDays,wvesEmp(:,1),'r--',1:nDays,wvesEmp(:,2),'r-','linew',2)


%% BIASES

c = 2; % coh

figure('position',[500 500 500 700]); 
for m=1:length(mods)
    subplot(3,1,m); hold on
    if m==3,d=2; else d=4; end
    plot(1:nDays,squeeze(choiceBias(m,c,d,:)),'k','linew',1.5);
%     plot(1:nDays,squeeze(confBias(m,c,d,:)),'r','linew',1.5);
%     plot(1:nDays,squeeze(RTBias(m,c,d,:)),'b','linew',1.5);

%     ylim([-5 5]);
end
plot(1:nDays,squeeze(choiceBias(3,c,1,:)),'c','linew',1.5);
plot(1:nDays,squeeze(choiceBias(3,c,3,:)),'g','linew',1.5);
% legend('Choice','Conf','RT')

if 0
corrcoef(squeeze(choiceBias(m,c,4,:)),squeeze(confBias(m,c,4,:)))
corrcoef(squeeze(choiceBias(m,c,4,:)),squeeze(RTBias(m,c,4,:)))
corrcoef(squeeze(RTBias(m,c,4,:)),squeeze(confBias(m,c,4,:)))
end

% figure; hold on;
% scatter(squeeze(choiceBias(m,c,4,:)),squeeze(confBias(m,c,4,:)),100,'k','.')
% lsline;
% scatter(squeeze(choiceBias(m,c,4,:)),squeeze(RTBias(m,c,4,:)),30,'r','.')
% scatter(squeeze(RTBias(m,c,4,:)),squeeze(confBias(m,c,4,:)),30,'b','.')

%% SHIFTS UNDER CONFLICT CONDITION



