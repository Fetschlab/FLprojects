function createBinaryFiles(subject,date)
% SJ 04-2022 

% to run in Command prompt

% insert subject and date as arguments to createBinaryFiles in the
% following call

% > start matlab -nosplash -nodesktop -nojvm -minimize -r "cd('C:\Users\fetschlab\Documents\MATLAB\dots3DMP\offline_processing\'), createBinaryFiles('lucio',20220526)" > out.txt 2> out.err &

addpath(genpath('C:\Program Files (x86)\Ripple\Trellis\Tools'));
addpath(genpath('C:\Users\fetschlab\Documents\MATLAB\dots3DMP'));

fpath = ['\\172.30.3.33\homes\fetschlab\data\' subject '\' subject '_neuro\' num2str(date) '\'];

load(fullfile(fpath,sprintf('%s%ddots3DMP_info.mat',subject,date)),'info');
processTrellisData(info,1);
exit;