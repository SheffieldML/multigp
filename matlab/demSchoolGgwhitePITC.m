% DEMSCHOOLGGWHITEPITC SCHOOL DATA EXPERIMENT with PITC and Ggwhite kernel

% MULTIGP

clear
clc
rand('twister', 1e7)
randn('state', 1e7)

dataSetName = 'schoolData2';
% Configuration of options
options = multigpOptions('pitc');
options.kernType = 'ggwhite';
options.optimiser = 'scg';
options.nlf = 1;
options.initialInducingPositionMethod =  'kmeansHeterotopic';
options.beta = 1e-3;
options.noiseOpt = 1;
options.tieOptions.selectMethod = 'nofree';
options.isArd = true;
options.fixInducing = false;
options.nVIKs = 1;


numActive = [ 5 20 50];


display = 0;
iters = 500;
totFolds = 10;

[totalError, elapsed_time_train] =  schoolSparseCore(dataSetName, options, ... 
       numActive, display, iters, totFolds);
   
save('schoolGgwhitePITC.mat')   
