
% count trials for a particular data struct, mainly for human aggregated
% data

fnames = fieldnames(data);

mods = unique(data.modality);
cohs = unique(data.coherence);
hdgs = unique(data.heading);
hdgs(hdgs==0) = [];
deltas = unique(data.delta);


nTrs = nan(length(mods),length(cohs),length(deltas)+1,length(hdgs),length(subjs));

for s = 1:length(subjs)

    temp = data;

    removethese = ~strcmp(data.subj,subjs{s});
    for f = 1:length(fnames)
        temp.(fnames{f})(removethese) = [];
    end

    disp([s length(unique(temp.subjDate))])


    for m = 1:length(mods)
        for c = 1:length(cohs)
            for d = 1:length(deltas)+1 % add extra column for all trials irrespective of delta

                for h = 1:length(hdgs)
                    
                    I = data.modality==mods(m) & data.coherence==cohs(c) & data.heading==hdgs(h) & strcmp(data.subj,subjs{s});

                    if d<length(deltas)+1
                        I = I & data.delta == deltas(d);
                    end

                    nTrs(m,c,d,h,s) = sum(I);

                end
            end
        end
    end


%     disp(length(temp.heading))

end


summaryTr = [];
summaryLabels = {'min','max','mean','median','std'};
nTrs(nTrs==0) = NaN;
for s = 1:length(subjs)
    theseSubj_trials = nTrs(:,:,2,:,s);
    summaryTr(s,:) = [min(theseSubj_trials(:)), max(theseSubj_trials(:)), nanmean(theseSubj_trials(:)), nanmedian(theseSubj_trials(:)) nanstd(theseSubj_trials(:))];
end

