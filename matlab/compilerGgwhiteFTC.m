% COMPILERGGWHITEFTC Runs the COMPILER DATA EXPERIMENT with full gp

% MULTIGP

clc
clear

dataSetName = 'compilerData';

options = multigpOptions('ftc');
options.kernType = 'ggwhite';
options.optimiser = 'scg';
options.nlf = 1;
options.tieOptions.selectMethod = 'free';
options.noise = 1e3;
options.isArd = true;

numTraining = [16 32 64 128];

display = 0;
iters = 1000;
totFolds = 10;

[totalError, elapsed_time_train, elapsed_time_test] =  ... 
    compilerFullCore(dataSetName, options, numTraining, ... 
        display, iters, totFolds);
   
save('compilerGgwhiteFTC.mat')   
  
