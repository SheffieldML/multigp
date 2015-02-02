% SCHOOLGGDTC Runs the SCHOOL DATA EXPERIMENT with DTC and GG kernel

% MULTIGP

clear
clc
rand('twister', 1e6)
randn('state', 1e6)
dataSetName = 'schoolData2';
% Configuration of options
options = multigpOptions('dtc');
options.kernType = 'gg';
options.optimiser = 'scg';
options.nlf = 1;
options.initialInducingPositionMethod =  'kmeansHeterotopic';
options.beta = 1e-2;
options.noiseOpt = 1;
options.tieOptions.selectMethod = 'nofree';
options.isArd = false;
options.fixInducing = false;

numActive = [5 20 50];

% Options for the batch
typeOfInit = 0;
isScaled = 1;
numberSeed = 1e3;
isBatch = 1;

display = 1;
iters = 5;
totFolds = 5;

[totalError, elapsed_time_train, totalErrorSMSE, totalErrorSMLL] =  schoolSparseCore(dataSetName, options, ... 
       numActive, display, iters, totFolds, typeOfInit, isScaled, numberSeed, isBatch);
   
save('schoolGgDTC.mat')   
