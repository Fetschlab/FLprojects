% RF Mapping

clear;clc;close all
addpath(genpath('/Users/stevenjerjian/Desktop/FetschLab/Analysis/codes/'))

% Load in the data

subject   = 'lucio';
dateRange = 20220307:20220804;

dataPath = '/Users/stevenjerjian/Desktop/FetschLab/Analysis/data/';
dataFileName = sprintf('%s_%d-%d_neuralData.mat',subject,dateRange(1),dateRange(end));

load(fullfile(dataPath,dataFileName));

% keep only those units which have sufficient RFMapping data
parSelect  = {'RFMapping'}; 
minRate    = 5;
minTrs     = 3;
dataStruct = dots3DMP_NeuralStruct_runCleanUp(dataStruct,parSelect,minRate,minTrs);

% CHECK THAT WE HAVE THE TARGET R AND THETA FOR NON_MOUSE POS
%%

% WORK IN PROGRESS

