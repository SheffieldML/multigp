clc
clear

addToolboxes(0,1)
dataSetName = 'compilerData';

options = multigpOptions('fitc');
options.kernType = 'gg';
options.optimiser = 'scg';
options.nlf = 1;
options.initialInducingPositionMethod = 'randomDataIsotopic';
options.beta = 1e3;
options.tieOptions.selectMethod = 'free';
options.isArd = true;
options.nVIKs = 1;
options.fixInducing = false;

numTraining = [16 32 64 128];

numActive{1} = [ 8 16];
numActive{2} = [16 32];
numActive{3} = [16 32 64];
numActive{4} = [32 64 128];

display = 0;
iters = 1000;
totFolds = 10;

[totalError, elapsed_time_train, elapsed_time_test] =  ... 
    compilerSparseCore(dataSetName, options, numTraining, ... 
       numActive, display, iters, totFolds);
   
save('compilerGgFITC.mat')   