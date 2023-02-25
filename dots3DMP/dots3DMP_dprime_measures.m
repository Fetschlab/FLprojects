function [d, meta_d] = dots3DMP_dprime_measures(data,mods,cohs,deltas)

% dprime and meta-dprime across subjs, for individual subjects, and for each condition

subjs = unique(data.subj);
nsubj = length(subjs);

[d.dprime{1},meta_d.dprime{1},d.crit{1},meta_d.crit{1}] = deal(nan(nsubj,1));

N = nan(length(mods),length(cohs),length(deltas)+1,nsubj);

d.dprime{2} = N;
d.crit{2}   = N;
meta_d.dprime{2} = N;
meta_d.crit{2} = N;

for s = 1:length(subjs)+1

    I = ~isnan(data.choice) & ~isnan(data.conf);

    if s < length(subjs)+1
        I = I & strcmp(data.subj,subjs{s});
    end

    [d.dprime{1}(s), d.crit{1}(s)] = dprime(double(data.heading(I)>0),data.choice(I)-1);
    [meta_d.dprime{1}(s), meta_d.crit{1}(s)] = dprime(double(data.correct(I)), double(data.conf(I)>=median(data.conf(I))));

    % repeat for individual conditions
    for m=1:length(mods)
        for c = 1:length(cohs)
            for d = 1:length(deltas)+1
                if d == length(deltas)+1
                    J = temp.modality==mods(m) & temp.coherence==cohs(c);
                else
                    J = temp.modality==mods(m) & temp.coherence==cohs(c) & temp.delta==deltas(d);
                end

                if s < length(subjs)+1
                    J = J & strcmp(data.subj,subjs{s});
                end


                [d.dprime{2}(s), d.crit{2}(s)] = ...
                    dprime(double(temp.heading(J)>0),temp.choice(J)-1);

                [meta_d.dprime{2}(s), meta_d.crit{2}(s)] = ...
                    dprime(double(temp.correct(J)), double(temp.conf(J)>=median(temp.conf(J))));

            end
        end
    end


end
      
