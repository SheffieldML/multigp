clear
clc
rand('twister', 1e6)
randn('state', 1e6)
addToolboxes(0,1)
file = 'Cd';
nlf = 1;
numFolds = 10;
iters = 100;
options = multigpOptions('ftc');
options.kernType = 'gg';
options.optimiser = 'scg';
options.nlf = 1;
options.beta = 1e3;
options.isArd = false;
[maerror, elapsed_time] = demJuraBatch(file, options, numFolds, iters);

save('scriptBatchJuraGgFullCdIsotopic' , 'maerror', 'elapsed_time')