function [total_t, maerror] = batch_demJura(options, data, ntrainX, ntrainX2, approx, missingData, iters)

% BATCH_DEMJURA 

% MULTIGP
% Setup model
[options, MXtrain, Xtrain, MYtrain, Ytrain, Xtest, Ytest] = multigpOptionsJura(options, data, ntrainX,missingData, ntrainX2, approx);
options.optimiser = 'scg';
options.includeInd = 0;
% Create the model
model = multigpCreate(q, data.nout, MYtrain, MXtrain, options);
% Change the initial parameters
params = multigpExtractParam(model);
params = params*rand + rand;
model = multigpExpandParam(model, params);
% Train the model
ini_t = cputime;
model = multigpOptimise(model,1,iters);
total_t = cputime - ini_t;
% Compute the error
mu = multigpPredictionMeanVar(model, Xtest, options);
maerror = mean(abs((Ytest{1} - mu{1})));



