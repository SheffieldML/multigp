% SCHOOLGGDTC Runs the SCHOOL DATA EXPERIMENT with DTC and IMC kernel

% MULTIGP

clear
clc
rand('twister', 1e6)
randn('state', 1e6)
dataSetName = 'schoolData2';
% Configuration of options
options = multigpOptions('dtc');
options.kernType = 'lmc';
options.optimiser = 'scg';
options.nlf = 1;
options.rankCorregMatrix = 1;
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
options.kern.rankCorregMatrix = options.rankCorregMatrix;

numActive = [ 5 20 50];
nlfOpts = [1 2 3 5 10];

display = 0;
iters = 200;
totFolds = 10;
totalErrorNlf = cell(length(nlfOpts),1);
elapsedTime = cell(length(nlfOpts),1);
for nlf=1:length(nlfOpts)
    fprintf('NUMBER OF LATENT FORCES: %d \n',nlfOpts(nlf));
    options.nlf = nlfOpts(nlf);
    [totalError, elapsed_time_train] =  schoolSparseCore(dataSetName, options, ...
        numActive, display, iters, totFolds);
    totalErrorNlf{nlf} =  totalError;
    elapsedTime{nlf} = elapsed_time_train;
    save('schoolGgDTC_SLFM.mat')
end