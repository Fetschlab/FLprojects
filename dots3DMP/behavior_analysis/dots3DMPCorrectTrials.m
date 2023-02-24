function correct = dots3DMPCorrectTrials(choice,heading,delta)
%
% 04/2022 SJ realized bug in code for zero heading rewards...
% because heading == 0 is recoded as eps in gui to set heading vector, all
% zero headings are technically tiny rightward, and therefore right choices
% were always rewarded on zero heading, and left was never rewarded!

% for human task this is not a big deal since they don't get trial-by-trial reward

% for monkey task, any data before 04/07/2022 inclusive could have this issue!

% either way, this function post-hoc fixes the data.correct variable based
% on heading, delta, and choice

% reward randomly if abs heading is less than delta, or heading is 0

% why fix rng? for reproducibility
rng(1);

temp           = (abs(heading)<=abs(delta) | heading==0);
randCorr       = false(size(temp));
randCorr(temp) = rand(sum(temp),1) < 0.5;

rightCorr = heading>0 & choice==2;
leftCorr  = heading<0 & choice==1;

correct = randCorr | rightCorr | leftCorr;

