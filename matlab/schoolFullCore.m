function [totalError, elapsed_time] = schoolFullCore(dataSetName, ...
    options, display, iters, totFolds, numberSeed, isBatch)

% SCHOOLFULLCORE description. 

% MULTIGP


addToolboxes(0,1)
warning('off', 'multiKernParamInit:noCrossKernel')

if nargin < 7
    isBatch = true;
    if nargin < 6
        numberSeed = 1e6;
    end
end

[XTemp, yTemp] = mapLoadData(dataSetName);

% XTemp = XTemp(1:5);
% yTemp = yTemp(1:5);


X = cell(1,size(yTemp, 2)+options.nlf);
y = cell(1,size(yTemp, 2)+options.nlf);
XTest = cell(1,size(yTemp, 2)+options.nlf);
yTest = cell(1,size(yTemp, 2));
q = size(XTemp{1}, 2)-1;
d = size(yTemp, 2)+options.nlf;

rand('twister',numberSeed);
SEEDS = rand('twister');
SEEDS = SEEDS(1:totFolds);

if isBatch
    totalError = zeros(totFolds,1);
    elapsed_time = zeros(totFolds,1);    
    for nFolds =1:totFolds
        rand('twister', double(SEEDS(nFolds)));
        fprintf('FOLD NUMBER: %d \n', nFolds);
        for i=1:options.nlf
            y{i} = [];
            X{i} = ones(1,q);
            options.bias(i) = 0;
            options.scale(i) = 1;
        end
        for i = 1:size(yTemp, 2)
            maxSamples = size(yTemp{i}, 1);
            randIndex = randperm(maxSamples);
            samplesTraining = floor(0.75*maxSamples);
            indexTraining = randIndex(1:samplesTraining);
            indexTesting = randIndex(samplesTraining+1:end);
            y{i+options.nlf} = yTemp{i}(indexTraining,:);
            X{i+options.nlf} = XTemp{i}(indexTraining,1:q);
            options.nRepeats{i} = XTemp{i}(indexTraining,q+1);
            yTest{i} = yTemp{i}(indexTesting,:);
            XTest{i+options.nlf} = XTemp{i}(indexTesting,1:q);
            options.bias(i+options.nlf) = mean(y{i+options.nlf});
            options.scale(i+options.nlf) = std(y{i+options.nlf});
        end
        % Creates the model
        model = multigpCreate(q, d, X, y, options);
        % Change initial parameters
        params = modelExtractParam(model);
        index = paramNameRegularExpressionLookup(model, '.* white .* variance');
        params(index) = log(10);
        model = modelExpandParam(model, params);
        %Trains the model
        timeL = cputime;
        model = multigpOptimise(model, display, iters);
        elapsed_time(nFolds) = cputime - timeL;
        % Save the results.
        XTest{1} = XTest{2};
        CC = schoolResults(model, XTest, yTest, y(options.nlf+1:end));
        totalError(nFolds) = mean(CC);
        %save(['school' upper(options.kernType(1)) options.kernType(2:end) upper(options.approx)])
    end
else
    rand('twister', numberSeed);
    fprintf('ONLY ONE FOLD. \n');
    for i=1:options.nlf
        y{i} = [];
        X{i} = ones(1,q);
        options.bias(i) = 0;
        options.scale(i) = 1;
    end
    for i = 1:size(yTemp, 2)
        maxSamples = size(yTemp{i}, 1);
        randIndex = randperm(maxSamples);
        samplesTraining = floor(0.75*maxSamples);
        indexTraining = randIndex(1:samplesTraining);
        indexTesting = randIndex(samplesTraining+1:end);
        y{i+options.nlf} = yTemp{i}(indexTraining,:);
        X{i+options.nlf} = XTemp{i}(indexTraining,1:q);
        options.nRepeats{i} = XTemp{i}(indexTraining,q+1);
        yTest{i} = yTemp{i}(indexTesting,:);
        XTest{i+options.nlf} = XTemp{i}(indexTesting,1:q);
        options.bias(i+options.nlf) = mean(y{i+options.nlf});
        options.scale(i+options.nlf) = std(y{i+options.nlf});
    end
    % Creates the model
    model = multigpCreate(q, d, X, y, options);
    % Change initial parameters
    params = modelExtractParam(model);
    index = paramNameRegularExpressionLookup(model, 'Beta .*');
    params(index) = log(1/(10)^2);
    model = modelExpandParam(model, params);
    %Trains the model
    timeL = cputime;
    model = multigpOptimise(model, display, iters);
    elapsed_time = cputime - timeL;
    % Save the results.
    XTest{1} = XTest{2};
    CC = schoolResults(model, XTest, yTest, y(options.nlf+1:end));
    totalError = mean(CC);
end

