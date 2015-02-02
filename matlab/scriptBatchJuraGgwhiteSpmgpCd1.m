clear
clc
addToolboxes(0,1)
file = 'Cd';
numActive =[50 100 200 500 600];
approx = 'dtc';
nlf = 1;
numFolds = 10;
iters = 1000;
options = multigpOptions('dtc');
options.kernType = 'ggwhite';
options.optimiser = 'scg';
options.nlf = 1;
options.initialInducingPositionMethod = 'randomDataIsotopic';
options.isArd = 1;
options.nVIKs = 1;
options.flagVIKs = false;
options.fixInducing = 0;
[maerror, elapsed_time] = demSpmgpJuraBatch(file, options, numActive, iters, numFolds);
save('scriptBatchJuraGgwhiteSpmgpCd1', 'maerror', 'elapsed_time')
