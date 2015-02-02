clear
clc
rand('twister', 1e6)
randn('state', 1e6)
addToolboxes(0,1)
file = 'Cd';
isScaled = 1;
isBatch  = 1;
numberSeed = 1e3;
typeOfInit = 1;
numActive = [50 100 200];
numFolds = 5;
iters = 5;
options = multigpOptions('pitc');
options.kernType = 'gg';
options.optimiser = 'scg';
options.nlf = 2;
options.initialInducingPositionMethod = 'kmeansHeterotopic';
options.isArd = 1;
options.tieInducing = 1;
options.fixInducing = 0;

[maerror, elapsed_time] = demSpmgpJuraBatch(file, options, numActive, ...
    iters, numFolds, typeOfInit, isScaled, numberSeed, isBatch);

save('scriptBatchJuraGgSpmgpCd1', 'maerror', 'elapsed_time')
