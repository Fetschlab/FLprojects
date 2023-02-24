% HUMAN DATA


% 03-2022

clear; clc; close all
cd /Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/dataStructs
addpath(genpath('/Users/stevenjerjian/Desktop/FetschLab/Analysis/codes/'))

%% select subject, load the data

subject = 'human';

conftask = 1;
RTtask   = 0; % change this to select RT or non-RT dataset


% load('human_20200213-20220113_RT_clean_Jan2022.mat') % human RT, Jan 2022
load('human_20200213-20220317_RT_clean_Mar2022.mat') % human RT, Mar 2022

RTlims = [0.25 2.5];

fnames = fieldnames(data);


if RTtask
    removethese = data.RT < RTlims(1) | data.RT > RTlims(2);
    
    for f=1:length(fnames)
        data.(fnames{f})(removethese) = [];
    end
end

% if strcmp(subject,'human') && RTtask
%     % not enough good data, so let's just remove for now?
%     removethese = data.heading==0;
%     fnames = fieldnames(data);
%     for f=1:length(fnames)
%         data.(fnames{f})(removethese) = [];
%     end
% end

mods   = unique(data.modality); 
cohs   = unique(data.coherence); 
deltas = unique(data.delta);
% deltas = [-3 3];
hdgs   = unique(data.heading);


%% parse data for individual subjects, and basic check on which subjects might merit individual analysis

clear parsedData* gfit* allCounts nTotalTrs
subjs = unique(data.subj);

for s = 1:length(subjs)
    sdata = data;
    removethese = ~strcmp(data.subj,subjs{s});
    
    for f=1:length(fnames)
        sdata.(fnames{f})(removethese) = [];
    end
    
    parsedDataIndivSubj(s) = dots3DMP_parseData(sdata,mods,cohs,deltas,hdgs,conftask,RTtask); 
    gfitIndivSubj(s) = dots3DMP_fit_cgauss(sdata,mods,cohs,deltas,conftask,RTtask); 

    allCounts(:,:,:,:,s) = parsedDataIndivSubj(s).n;
    nTotalTrs(s) = length(sdata.heading);
end

allCounts(1,2,:,:,:) = 0; % vestibular trials are double counted in n!
subjsForIndivAnalysis = nTotalTrs>750;
subjs2keep = subjs(subjsForIndivAnalysis);

modCounts = permute(allCounts(:,:,1:3,:,:),[1 5 2 3 4]);

modCounts = reshape(modCounts,size(modCounts,1),size(modCounts,2),size(modCounts,3)*size(modCounts,4)*size(modCounts,5));
modCounts = sum(modCounts,3);

if 0
clear cols
cols{1} = 'kmc';
cols{2} = 'krb';
figure;
hb=bar(modCounts','stacked');
for hh=1:length(hb)
    hb(hh).FaceColor = cols{2}(hh); hb(hh).EdgeColor = cols{2}(hh);
end
set(gca,'xticklabel',subjs)

figure;
for m=1:3
    subplot(3,1,m); hold on;
    for c = 1:2
        temp = squeeze(allCounts(m,c,2,:,:));
        plot(hdgs,temp,'color',cols{c}(m));
        
        for s = 1:length(subjs)
            text(hdgs(end)-1,temp(end,s)*1.05,subjs{s});
        end
            
    end
end
end


%% plot individual subject logistic fits, assess data quality

if 0
for s = 1%:length(subjs)
    fh = dots3DMP_plots(parsedData(s),mods,cohs,deltas,hdgs,conftask,RTtask);
    savefig(fh(1),sprintf('dots3DMP_logisticfits_delta0_%s',subjs{s}));
    savefig(fh(2),sprintf('dots3DMP_logisticfits_deltas_%s',subjs{s}));
    close all;
end
end

%%
clear ChoiceShift*
for s = 1:length(subjs)
        
    ChoiceShift(s,:,1) = gfit(s).choice.mu(3,:,1)-gfit(s).choice.mu(3,:,2);
    ChoiceShift(s,:,2) = gfit(s).choice.mu(3,:,3)-gfit(s).choice.mu(3,:,2);
    
    ChoiceShiftSE(s,:,:) = gfit(s).choice.muSE(3,:,[1 3]);
end

ChoiceShift = ChoiceShift(subjsForIndivAnalysis,:,:);
ChoiceShiftSE = ChoiceShiftSE(subjsForIndivAnalysis,:,:);

figure;
for c = 1:length(cohs)
    for d = 1:2
        subplot(2,2,d+2*(c-1)); hold on
        plot(ChoiceShift(:,c,d),'k');
        plot(ChoiceShiftSE(:,c,d),'r');
    end
end

% need more data from SJJ and DJB...
