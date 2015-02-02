clear
clc
addToolboxes(0,1)
file = 'Cd';
nlf = 1;
numFolds = 10;
iters = 1000;
options = multigpOptions('ftc');
options.kernType = 'ggwhite';
options.optimiser = 'scg';
options.nlf = 1;
options.beta = 1e3;
options.isArd = 1;
[maerror, elapsed_time] = demJuraBatch(file, options, numFolds, iters);

save('scriptBatchJuraGgwhiteFullCd' , 'maerror', 'elapsed_time')