function [maerror, elapsed_time] = demSpmgpJuraBatch(file, options, ...
    numActive, iters, numFolds, typeOfInit, isScaled, numberSeed, isBatch)

% DEMSPMGPJURABATCH Demonstrate sparse convolution models on JURA data.

% MULTIGP
if nargin < 9
    isBatch = true;
    if nargin < 8
        numberSeed = 10^6;
        if nargin < 7
            isScaled = true;
            if nargin < 6
                typeOfInit = 1;
            end
        end
    end
end

dataSetName = ['juraData' file];
display = 1;
maerror =  zeros(numFolds,length(numActive));
elapsed_time =  zeros(numFolds,length(numActive));

[XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData(dataSetName);
scaleVal = zeros(1,size(yTemp, 2));
biasVal = zeros(1,size(yTemp, 2));
for k =1:size(yTemp, 2),
    biasVal(k) = mean(yTemp{k});
    scaleVal(k) = sqrt(var(yTemp{k}));
end

options.bias = biasVal;
if isScaled
    options.scale = scaleVal;
    options.beta = 10;
else
    options.scale = ones(1, length(biasVal));
    options.beta = 1./((0.15*options.scale).^2);
end

X = cell(1,size(yTemp, 2));
y = cell(1,size(yTemp, 2));

for i = 1:size(yTemp, 2)
    y{i} = yTemp{i};
    X{i} = XTemp{i};
end

q = 2;
d = size(yTemp, 2);

if isBatch
    rand('twister',numberSeed);
    SEEDS = rand('twister');
    SEEDS = SEEDS(1:numFolds);
    for seudoPoints =1: length(numActive)
        options.numActive = numActive(seudoPoints);
        if strcmp(options.approx, 'dtcvar')
            if options.flagVIKs
                options.nVIKs = options.numActive;
            end
        end
        for fold = 1:numFolds
            rand('twister', double(SEEDS(fold)));
            fprintf('Creating MODEL %d DATASET %s APPROX %s NUMACTIVE %d \n',...
                fold, file, options.approx, options.numActive);
            model = multigpCreate(q, d, X, y, options);
            % Initialize parameters
            % Change initial conditions
            if typeOfInit
                params = modelExtractParam(model);
                index = paramNameRegularExpressionLookup(model, 'multi .* inverse .*');
                params(index) = log(100 +  10*rand(1,length(index)));
                model = modelExpandParam(model, params);
            end
            % Optimization procedure
            fprintf('Optimizing MODEL %d DATASET %s APPROX %s NUMACTIVE %d \n',...
                fold, file, options.approx, options.numActive)
            tic
            model = multigpOptimise(model, display, iters);
            elapsed_time(fold,seudoPoints) = toc;
            fprintf('Elapsed time optimization: %f hours\n', elapsed_time(fold,seudoPoints)/3600)
            % Likelihood
            loglik = modelLogLikelihood(model);
            fprintf('Model likelihood: %f\n', loglik);
            % Prediction
            mu = multigpPosteriorMeanVar(model, XTestTemp);
            maerror(fold, seudoPoints) = mean(abs((yTestTemp{1} - mu{model.nlf + 1})));
        end
    end
else
    options.numActive = numActive;
    if strcmp(options.approx, 'dtcvar')
        if options.flagVIKs
            options.nVIKs = options.numActive;
        end
    end
    rand('twister',numberSeed);
    fprintf('Creating MODEL DATASET %s APPROX %s NUMACTIVE %d \n',...
        file, options.approx, options.numActive);
    model = multigpCreate(q, d, X, y, options);
    % Initialize parameters
    % Change initial conditions
    if typeOfInit
        params = modelExtractParam(model);
        index = paramNameRegularExpressionLookup(model, 'multi .* inverse .*');
        params(index) = log(100 +  10*rand(1,length(index)));
        model = modelExpandParam(model, params);
    end
    % Optimization procedure
    fprintf('Optimizing MODEL DATASET %s APPROX %s NUMACTIVE %d \n', ...
        file, options.approx, options.numActive)
    tic
    model = multigpOptimise(model, display, iters);
    elapsed_time = toc;
    fprintf('Elapsed time optimization: %f hours\n', elapsed_time/3600)
    % Likelihood
    loglik = modelLogLikelihood(model);
    fprintf('Model likelihood: %f\n', loglik);
    % Prediction
    mu = multigpPosteriorMeanVar(model, XTestTemp);
    maerror = mean(abs((yTestTemp{1} - mu{model.nlf + 1})));
end



