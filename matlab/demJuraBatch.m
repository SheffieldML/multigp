function [maerror, elapsed_time] = demJuraBatch(file, options, numFolds, iters)

% DEMJURABATCH Demonstrate convolution models on JURA data.

% MULTIGP

testHeterotopic = true;
dataSetName = ['juraData' file];
display = 1;
maerror =  zeros(numFolds,1);
elapsed_time =  zeros(numFolds,1);
[XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData(dataSetName);
scaleVal = zeros(1,size(yTemp, 2));
biasVal = zeros(1,size(yTemp, 2));
for k =1:size(yTemp, 2),
    biasVal(k) = mean(yTemp{k});
    scaleVal(k) = sqrt(var(yTemp{k}));
end
options.bias =  [zeros(1, options.nlf) biasVal];
options.scale = [zeros(1, options.nlf) scaleVal];

X = cell(size(yTemp, 2)+options.nlf,1);
y = cell(size(yTemp, 2)+options.nlf,1);
for j=1:options.nlf
    y{j} = [];
    X{j} = zeros(1, 2);
end
for i = 1:size(yTemp, 2)
    y{i+options.nlf} = yTemp{i};
    X{i+options.nlf} = XTemp{i};
end
XTest = cell(size(yTemp, 2)+options.nlf,1);
for j=1:options.nlf
    XTest{j} = ones(1, 2);
end
for i = 1:size(yTemp, 2)
    XTest{i+options.nlf} = XTestTemp{i};
end

q = size(X{end},2);
d = size(yTemp, 2) + options.nlf;
rand('twister', 1e6)
randn('state', 1e6)
for fold =1: numFolds
    % Creates the model
    fprintf('Creating MODEL %d DATASET %s APPROX %s\n', fold, file, options.approx)
    model = multigpCreate(q, d, X, y, options);
    % Change initial conditions
    params = modelExtractParam(model);
    if strcmp(options.kernType, 'lmc')
        index = paramNameRegularExpressionLookup(model, 'multi .* inverse .*');
        params(index) = log(100 +  5*rand(1,length(index)));
        index = paramNameRegularExpressionLookup(model, 'multi .* variance .*');
        params(index) = log(50);
        index = paramNameRegularExpressionLookup(model, 'multi .* A(.*');
        params(index) = 1+rand(1, length(index));
    else
        index = paramNameRegularExpressionLookup(model, 'multi .* inverse width output .*');
        params(index) = log(100 +  5*rand(1,length(index)));
        index = paramNameRegularExpressionLookup(model, 'multi .* variance');
        params(index) = log(10);
        index = paramNameRegularExpressionLookup(model, 'multi .* sensitivity');
        params(index) = 1+ rand(1, length(index));
    end
    model = modelExpandParam(model, params);
    % Optimization procedure
    fprintf('Optimizing MODEL %d DATASET %s APPROX %s\n', fold, file, options.approx)
    tic;
    model = multigpOptimise(model, display, iters);
    elapsed_time(fold) = toc;
    fprintf('Elapsed time optimization: %f hours\n', elapsed_time(fold)/3600)
    % Likelihood
    loglik = modelLogLikelihood(model);
    fprintf('Model likelihood: %f\n', loglik);
    % Prediction
    % Create a new model but using prediction and validation data
    if testHeterotopic
        model2 = constructModelHeterotopic(file, options, q, d);
        params = modelExtractParam(model);
        model2 = modelExpandParam(model2, params);
        mu = multigpPosteriorMeanVar(model2, XTest);
    else
        mu = multigpPosteriorMeanVar(model, XTest);
    end
    maerror(fold) = mean(abs((yTestTemp{1} - mu{model.nlf + 1})));    
end

function model = constructModelHeterotopic(file, options, q, d)

dataSetName = ['juraData' file];
dataSetName(end-7:end) = [];
[XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData(dataSetName);
scaleVal = zeros(1,size(yTemp, 2));
biasVal = zeros(1,size(yTemp, 2));
for k =1:size(yTemp, 2),
    biasVal(k) = mean(yTemp{k});
    scaleVal(k) = sqrt(var(yTemp{k}));
end
options.bias =  [zeros(1, options.nlf) biasVal];
options.scale = [zeros(1, options.nlf) scaleVal];
X = cell(size(yTemp, 2)+options.nlf,1);
y = cell(size(yTemp, 2)+options.nlf,1);
for j=1:options.nlf
    y{j} = [];
    X{j} = zeros(1, 2);
end
for i = 1:size(yTemp, 2)
    y{i+options.nlf} = yTemp{i};
    X{i+options.nlf} = XTemp{i};
end
XTest = cell(size(yTemp, 2)+options.nlf,1);
for j=1:options.nlf
    XTest{j} = ones(1, 2);
end
for i = 1:size(yTemp, 2)
    XTest{i+options.nlf} = XTestTemp{i};
end
q = size(X{end},2);
d = size(yTemp, 2) + options.nlf;
model = multigpCreate(q, d, X, y, options);
