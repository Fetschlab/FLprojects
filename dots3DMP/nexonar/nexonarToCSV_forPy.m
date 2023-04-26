

% SJ 02-15-2023
% temp for saving nexonar test data into csvs, for working with in Python

cd /Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/test/nexonar_tests-2023_02_16

nexfiles = dir('*_nexonar.mat');

newFnames = {'t1','t2','X','Y','Z'};




clear allNexData filename

overallTr = 1;
currTr = 0;
for f = 1:length(nexfiles)


    load(nexfiles(f).name)

    filetime = split(nexfiles(f).name,'_nexonar.mat');
    filetime = num2str(filetime{1}(end-3:end));

    numdatacols = size(nex.nexdata{t},2);
    pldapscols  = length(fields(nex.pldaps));
    condscols   = length(fields(nex.conditions));


    for t = 1:length(nex.nexdata)

        num_trs = size(nex.nexdata{t},1);

%         filename{currTr+1:currTr+num_trs,1} = nexfiles(f).name;

        allNexData(currTr+1:currTr+num_trs,1:numdatacols) = nex.nexdata{t};

        fnames = fieldnames(nex.pldaps);
        for fn = 1:length(fnames)
            allNexData(currTr+1:currTr+num_trs,numdatacols+fn) = nex.pldaps.(fnames{fn})(t);
            newFnames{numdatacols+fn} = fnames{fn};
        end

        fnames = fieldnames(nex.conditions);
        for fn = 1:length(fnames)
            allNexData(currTr+1:currTr+num_trs,numdatacols+pldapscols+fn) = nex.conditions.(fnames{fn})(t);
            newFnames{numdatacols+pldapscols+fn} = fnames{fn};
        end

        allNexData(currTr+1:currTr+num_trs,numdatacols+pldapscols+condscols+1) = overallTr;
        allNexData(currTr+1:currTr+num_trs,numdatacols+pldapscols+condscols+2) = f;
        currTr = size(allNexData,1);
        overallTr = overallTr+1;
    end


end

newFnames{end+1} = 'overall_trial_number';
newFnames{end+1} = 'file_number';


for f = 1:length(newFnames)
    S.(newFnames{f}) = allNexData(:,f);
end
struct2csv(S,'nex_20230216.csv')


%%
X.amps = [160 100 160 160 100 200 50 50]';
X.pad  = [0 0 0.9 0.9 0 0 0 0]';
X.dur  = [2 1 1.3 1.3 1.3 2 1 1]';
X.sigma = 0.14*ones(size(X.dur));
X.start_back = [80 50 80 100 80 100 25 50]';
struct2csv(X,'nex_20230216_motionsettings.csv')