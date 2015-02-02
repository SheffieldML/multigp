% DEMSPMGPGGTOY3KL Measures the KL divergence between the full Multigp with
% GG kernel against the approximation using PITC

% MULTIGP

clc
clear
addToolboxes(0,1)
rand('twister',1e6);
randn('state',1e6);
%
dataSetName = 'ggwhiteToyMissing';
experimentNo = 1;

try
    capName = dataSetName;
    capName(1) = upper(capName(1));
    load(['dem' capName num2str(experimentNo) '.mat'])
catch
    [XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData(dataSetName);
    options = multigpOptions('ftc');
    options.kernType = 'ggwhite';
    options.optimiser = 'scg';
    options.nlf = 1;
    X = cell(size(yTemp, 2)+options.nlf,1);
    y = cell(size(yTemp, 2)+options.nlf,1);
    for j=1:options.nlf
        y{j} = [];
        X{j} = 1;
    end
    for i = 1:size(yTemp, 2)
        y{i+options.nlf} = yTemp{i};
        X{i+options.nlf} = XTemp{i};
    end
    q = 1;
    d = size(yTemp, 2) + options.nlf;
    % Creates the model
    model = multigpCreate(q, d, X, y, options);
    params = modelExtractParam(model);
    index = paramNameRegularExpressionLookup(model, 'multi .* inverse .*');
    params(index) = log(100);
    model = modelExpandParam(model, params);
    display = 1;
    iters = 2000;
    % Trains the model
    model = multigpOptimise(model, display, iters);
    % Save the results.
    capName = dataSetName;
    capName(1) = upper(capName(1));
    save(['dem' capName num2str(experimentNo) '.mat'], 'model');
end
llFull = modelLogLikelihood(model);
% Trains the sparse model
nFolds = 10;
display = 1;
iters = 2000;
options = multigpOptions('dtc');
options.kernType = 'ggwhite';
options.optimiser = 'scg';
options.nlf = 1;
options.initialInducingPositionMethod = 'randomComplete';
options.numActive = 30;
options.fixInducing = 0;
options.beta = 1e-3;
numActive = [10 20 30 40 50 100 200];
llApprox = zeros(length(numActive), nFolds);
elapsed_time = zeros(length(numActive), nFolds);
models = cell(length(numActive), nFolds);
for j=1:length(numActive)
    options.numActive = numActive(j);
    for k = 1: nFolds,
        [XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData(dataSetName);
        X = cell(size(yTemp, 2),1);
        y = cell(size(yTemp, 2),1);
        for i = 1:size(yTemp, 2)
            y{i} = yTemp{i};
            X{i} = XTemp{i};
        end
        q = 1;
        d = size(yTemp, 2);
        % Creates the model
        model = multigpCreate(q, d, X, y, options);
        params = modelExtractParam(model);
        index = paramNameRegularExpressionLookup(model, 'multi .* inverse .*');
        params(index) = log(100);
        model = modelExpandParam(model, params);
        rand('twister',10^k);
        randn('state',10^k);
        % Change the variance (the amplitude of the kernel)
        if options.fixInducing == 0
            params = modelExtractParam(model);
            %for i = 1:model.nout,
            paramInd = paramNameRegularExpressionLookup(model, 'X_u .*');
            initialLoc = 0.05*randn(1,length(model.X_u));
            params(paramInd) = initialLoc;
            %end
            model = modelExpandParam(model, params);
        end
        % Train the model
        tic;
        model = multigpOptimise(model, display, iters);
        elapsed_time(j,k) = toc;
        llApprox(j,k) = modelLogLikelihood(model);
        models{j,k} = model;
    end
end

save('demSpmgpNoiseToy1KL.mat', 'models', 'llFull', 'llApprox', 'elapsed_time');







