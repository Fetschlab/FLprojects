function [muPMFboot,sigmaPMFboot,wvesEmpboot,wvesPredboot] = dots3DMP_cgauss_bootstrap_NN(data,mods,cohs,deltas,nboots)
% run cgauss fitting on data, bootstrapped with replacement to generate CIs

options = optimset('Display','off');

muPMFboot = cell(nboots,1);
sigmaPMFboot = cell(nboots,1);
wvesPredboot = nan(nboots,length(cohs));
wvesEmpboot = nan(nboots,length(cohs));

% define anonymous functions for fitting:

cgauss = @(b,hdg) 1/2 * ( 1 + erf( (hdg-b(1))./(b(2)*sqrt(2)) ) );
    % for probabilities, error is negative log likelihood of observing the data, which is
    % [ log(Pright(hdg)) + log(1-(~Pright(hdg))) ]
cgauss_err = @(param,choice,hdg) -(sum(log(cgauss(param,hdg(choice))))+sum(log(1-cgauss(param,hdg(~choice))))); 

unc = 0; % saves biases from fminunc instead of fminsearch (SEs always are fminunc, and plots are always fminsearch)

if nboots==0
    return
end

for n = 1:nboots
    
    if mod(n,10)==0
        fprintf('\nprogress: %d out of %d\n', n, nboots);
    end

    muPMFboot{n} = nan(length(mods),length(cohs),length(deltas));
    sigmaPMFboot{n} = nan(length(mods),length(cohs),length(deltas));

    % first the single-cues, and comb with either:
    % all trials irrespective of delta
    D = length(deltas)+1; % (the extra column we made for pooling across deltas)
    % OR select just delta=0:
    % D = find(deltas==0);
    for c = 1:length(cohs)
        for m = 1:length(mods)     % m c d h
            if m==1
                I = data.modality==mods(m);
            else
                if D==length(deltas)+1
                    I = data.modality==mods(m) & data.coherence==cohs(c); % all trials irrespective of delta
                else
                    I = data.modality==mods(m) & data.coherence==cohs(c) & data.delta==deltas(D);
                end
            end
            
            bootI = randsample(find(I),sum(I),true);

            if unc
                [beta,~,~,~,~,~] = fminunc(@(x) cgauss_err(x,data.choice(bootI)==1,data.heading(bootI)), [0 3], options);
            else
                beta = fminsearch(@(x) cgauss_err(x,data.choice(bootI)==1,data.heading(bootI)), [0 3], options);
            end
            muPMFboot{n}(m,c,D) = beta(1);
            sigmaPMFboot{n}(m,c,D) = beta(2);
            
        end
    end

    % then the comb only, separated by delta, to get biases
    for c = 1:length(cohs)
        for d = 1:length(deltas)     % m c d h
            I = data.modality==3 & data.coherence==cohs(c) & data.delta==deltas(d);

            bootI = randsample(find(I),sum(I),true);
            
            if unc
                [beta,~,~,~,~,~] = fminunc(@(x) cgauss_err(x,data.choice(bootI)==1,data.heading(bootI)), [0 3], options);
            else
                beta = fminsearch(@(x) cgauss_err(x,data.choice(bootI)==1,data.heading(bootI)), [0 3], options);
            end
            muPMFboot{n}(3,c,d) = beta(1);
            sigmaPMFboot{n}(3,c,d) = beta(2);
            
        end
    end
    
    % lastly compute the weights
    for c = 1:length(cohs)      % m c d
        wvesPredboot(n,c) = (1/sigmaPMFboot{n}(1,1,D)^2) / ((1/sigmaPMFboot{n}(1,1,D)^2) + (1/sigmaPMFboot{n}(2,c,D)^2));
                         % m c d
        actual(1) = (muPMFboot{n}(3,c,1)-muPMFboot{n}(3,c,2)+(deltas(1)/2)) / deltas(1);
        actual(2) = (muPMFboot{n}(3,c,3)-muPMFboot{n}(3,c,2)+(deltas(3)/2)) / deltas(3);    
        wvesEmpboot(n,c) = mean(actual);

    end

end
