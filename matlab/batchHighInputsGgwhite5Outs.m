% BATCHHIGHINPUTSGGWHITE2OUTS Computes DTC VAR bound and full gp likelihood 
% for two outputs

% MULTIGP 


clc
clear
rand('twister',1e3);
randn('state',1e3);
nOuts = 5;
inverseWidth = 5 + 10*rand(1,nOuts);
sensitivity = 1 + 5*rand(1,nOuts);
inputDim = [2 5 10];
nInducing = [50 100 200 500];
iters = 5000;
paramInit = 1e1;
methodInit = 'randomComplete';
methodInitMGP = false;
folds = 10;
nTrainingPoints = 100;
llApprox = cell(length(inputDim),1);
llFull = zeros(length(inputDim),1);
for i = 1:length(inputDim)
    [llFull(i), X, y, Kout, noisePerOutput] = generateDataset('ggwhite', inputDim(i), ...
        nOuts, nTrainingPoints, inverseWidth, sensitivity);
    llApprox{i} = highInputsGgwhite(X, y, nInducing, nOuts, ...
        iters, noisePerOutput, inverseWidth, sensitivity,...
        paramInit, methodInit, methodInitMGP, folds);
end
