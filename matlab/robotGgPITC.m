clear
clc
rand('twister', 1e6)
randn('state', 1e6)
addToolboxes(0,1)
dataSetName ='robotWireless';
nlf =[1 2 3 4];
iters = 1000;
options = multigpOptions('dtcvar');
options.kernType = 'gg';
options.optimiser = 'scg';
options.nlf = 1;
options.initialInducingPositionMethod = 'kmeansHeterotopic';
options.isArd = false;
options.fixInducing = 0;
options.numActive = 50;
options.beta = 1e-1;
options.normalise = false;

[maerror, elapsed_time] = robotSparseCore(dataSetName, options, nlf, iters);

save([dataSetName upper(options.approx)], 'maerror', 'elapsed_time')
