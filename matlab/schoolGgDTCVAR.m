clear
clc
rand('twister', 1e6)
randn('state', 1e6)
addToolboxes(0,1)
dataSetName = 'schoolData2';
% Configuration of options
options = multigpOptions('dtcvar'); 
options.kernType = 'gg';
options.optimiser = 'scg';
options.nlf = 1;
options.initialInducingPositionMethod =  'kmeansHeterotopic';
options.beta = 1e3;
options.noiseOpt = 1;
options.tieOptions.selectMethod = 'nofree';
options.isArd = true;
options.fixInducing = false;

numActive = [ 5 20 50 100];


display = 0;
iters = 200;
totFolds = 10;

[totalError, elapsed_time_train] =  schoolSparseCore(dataSetName, options, ... 
       numActive, display, iters, totFolds);
   
save('schoolGgDTCVAR.mat')   