function dots3DMP_dprime_measures(data,mods,cohs,deltas)

% dprime and meta-dprime for individual subjects, and for each condition

subjs = unique(data.subj);
fnames = fieldnames(data);
% across all subjs

meta_d_subjCond = nan(length(mods),length(cohs),length(deltas)+1,length(subjs));
meta_crit_subjCond = nan(length(mods),length(cohs),length(deltas)+1,length(subjs));

for s = 1:length(subjs)+1
    if s == length(subjs)+1
        
        [dp(s), crit(s)] = dprime(double(data.heading>0),data.choice-1);

        [meta_d(s), meta_crit(s)] = dprime(double(data.correct), double(data.conf>=median(data.conf)));

    else
        temp = data;
        for f = 1:length(fnames)
            temp.(fnames{f})(~strcmp(data.subj,subjs{s})) = [];
        end
        
        % confidence should be set to median withi each subject and within
        % each condition - i.e. criterion for that condition
%         medianconf = median(temp.conf); 
        
        [meta_d(s), meta_crit(s)] = dprime(double(temp.correct), double(temp.conf>=median(temp.conf)));
        [dp(s), crit(s)] = dprime(double(temp.heading>0),temp.choice-1);

        
        for m=1:length(mods)
            for c = 1:length(cohs)
                for d = 1:length(deltas)+1
                    if d == length(deltas)+1
                        J = temp.modality==mods(m) & temp.coherence==cohs(c);
                    else
                        J = temp.modality==mods(m) & temp.coherence==cohs(c) & temp.delta==deltas(d);
                    end
                    
                    [meta_d_subjCond(m,c,d,s), meta_crit_subjCond(m,c,d,s)] = ...
                        dprime(double(temp.correct(J)), double(temp.conf(J)>=median(temp.conf(J))));
                    
                    [d_subjCond(m,c,d,s), crit_subjCond(m,c,d,s)] = ...
                        dprime(double(temp.heading(J)>0),temp.choice(J)-1);

                    
                end
            end
        end


    end
end
