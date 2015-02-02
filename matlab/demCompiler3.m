% demCompiler3.m

% Demo of Sparse Multi Output Gaussian Process using FITC. In this demo, we use the
% Gaussian Kernel for all the covariances (or Kernels) involved and only one hidden function.
% When changing the kernel, the fix values in this code and the indeces in
% the tieParam vector in multigpCreate must also be changed.

% Mauricio Alvarez 2008
clc
clear
% rand('seed',1e5);
% randn('seed',1e5);
%
dataSetName = 'data_compiler_org.mat';
experimentNo = 1;
ntrainX =128;
tryNumber = 1;
ntrainX2 =32;
iters =3000;
approx = 'pitc';
saveFigures=1;
% load data
data = multigpLoadData(dataSetName);
data.nin = 1;  % Number of latent functions
missingData = cell(data.nout,1);


options.isSparse = 1;  % Indicates if the scheme is sparse (1) or if it is full (0)
options.stdBiasX = 0;  % Indicates if the inputs are to be standarized
options.stdBiasY = 1;  % Indicates if the outputs are to be standarized
options.outKernelName = 'gg'; % Indicates the name of the output kernel: lfm, sim or gg

% Setup model

[options, MXtrain, Xtrain, MYtrain, Ytrain, Xtest, Ytest] = multigpOptionsCompiler(options, data, ntrainX, ntrainX2, ...
    approx, missingData, tryNumber);

options.optimiser = 'scg';

% Creates the model
model = multigpCreate(Ytrain, Xtrain, options);

% Trains the model and counts the training time
ini_t = cputime;
[model, params] = multigpOptimise(model,1,iters);
total_t = cputime - ini_t;

% This part is to compute the prediction error

maxTest = length(Ytest{1});
numPart = 1000;
step = floor(maxTest/numPart);
remain = rem(maxTest,numPart);
XtestMod = cell(model.nout,1);
mserror = zeros(maxTest,1);
for j =1:(numPart+1),
    if (j==numPart+1)
        indexes = indexes(end)+1:maxTest;
    else
        indexes = (j-1)*step+1:j*step;
    end
    for k= 1:model.nout,
        XtestMod{k} = Xtest{k}(indexes,:); 
    end
    [mu, varsigma] = multigpPredictionMeanVar(model, XtestMod, options);
    for k= 1:model.nout,
        mserror(indexes,k) = abs((Ytest{k}(indexes) - mu{k}));
    end
end


% [mu, varsigma] = multigpPredictionMeanVar(model, Xtest, Options);
% desv_sserror = zeros(model.nout,1);
% mean_sserror = zeros(model.nout,1);
% desv_logprob = zeros(model.nout,1);
% mean_logprob = zeros(model.nout,1);
% 
% % Standarized performance measures
% for k= 1:model.nout,
%     mserror = abs((Ytest{k} - mu{k}));
%     mean_sserror(k) = mean(mserror);
%     desv_sserror(k) =  std(mserror);
% end
% 
% 
