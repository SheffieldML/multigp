% SCHOOLGGDTC Runs the SCHOOL DATA EXPERIMENT with DTC and IMC kernel

% MULTIGP

clear
clc
rand('twister', 1e6)
randn('state', 1e6)
dataSetName = 'schoolData2';
% Configuration of options
options = multigpOptions('fitc');
options.kernType = 'lmc';
options.optimiser = 'scg';
options.nlf = 1;

options.initialInducingPositionMethod =  'kmeansHeterotopic';
options.beta = 1e-3;
options.noiseOpt = 1;
options.tieOptions.selectMethod = 'nofree';
options.isArd = true;
options.kern.isArd = options.isArd;
options.fixInducing = false;
options.includeNoise = false;
options.gamma = exp(-2);
options.kern.nout = 139;

numActive = [ 5 20 50];
rankOpts = [1 2 3 5 10];

display = 0;
iters = 200;
totFolds = 10;
totalErrorNlf = cell(length(rankOpts),1);
elapsedTime = cell(length(rankOpts),1);
for rank=1:length(rankOpts)
    options.rankCorregMatrix = rankOpts(rank);
    options.kern.rankCorregMatrix = options.rankCorregMatrix;
    [totalError, elapsed_time_train] =  schoolSparseCore(dataSetName, options, ...
        numActive, display, iters, totFolds);
    totalErrorNlf{rank} =  totalError;
    elapsedTime{rank} = elapsed_time_train;
    save('schoolGgFITC_ICM.mat')
end