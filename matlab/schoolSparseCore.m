function [totalError, elapsed_time_samples, totalErrorSMSE, totalErrorMSLL] = ...
    schoolSparseCore(dataSetName, options, numActive, display, iters, ...
    totFolds, typeOfInit, isScaled, numberSeed, isBatch)

% SCHOOLSPARSECORE description.

% MULTIGP
if nargin < 10
    isBatch = 1;
    if nargin < 9
        numberSeed = 1e6;
        if nargin < 8
            isScaled = 1;
            if nargin < 7
                typeOfInit = 0;
            end
        end
    end
end

[XTemp, yTemp] = mapLoadData(dataSetName);

X = cell(1,size(yTemp, 2));
y = cell(1,size(yTemp, 2));
XTest = cell(1,size(yTemp, 2));
yTest = cell(1,size(yTemp, 2));
q = size(XTemp{1}, 2)-1;
d = size(yTemp, 2);

warning('off', 'multiKernParamInit:noCrossKernel');

if isBatch
    rand('twister',numberSeed);
    SEEDS = rand('twister');
    SEEDS = SEEDS(1:totFolds);
    elapsed_time_samples = zeros(totFolds, length(numActive));
    totalError = zeros(totFolds, length(numActive));
    scaleVal = zeros(1, size(yTemp, 2));
    if nargout > 2
        totalErrorSMSE = zeros(totFolds, length(numActive));
        totalErrorMSLL = zeros(totFolds, length(numActive));
    end
    for j =1:length(numActive),
        fprintf('NUMBER OF ACTIVE POINTS: %d \n',numActive(j));
        options.numActive = numActive(j);        
        for nFolds =1:totFolds
            rand('twister', double(SEEDS(nFolds)));
            fprintf('FOLD NUMBER: %d \n', nFolds);
            for i = 1:size(yTemp, 2)
                maxSamples = size(yTemp{i}, 1);
                randIndex = randperm(maxSamples);
                samplesTraining = floor(0.75*maxSamples);
                indexTraining = randIndex(1:samplesTraining);
                indexTesting = randIndex(samplesTraining+1:end);
                y{i} = yTemp{i}(indexTraining);
                X{i} = XTemp{i}(indexTraining,1:q);
                options.nRepeats{i} = XTemp{i}(indexTraining,q+1);
                yTest{i} = yTemp{i}(indexTesting,:);
                XTest{i} = XTemp{i}(indexTesting,1:q);
                options.bias(i) = mean(y{i});
                scaleVal(i) = std(y{i});  
            end
            if isScaled
                options.scale = scaleVal;
                options.beta = 10;
            else
                options.scale = ones(1, length(scaleVal));
                options.beta = 1./((0.15*options.scale).^2);
            end
            % Creates the model
            model = multigpCreate(q, d, X, y, options);
            % Change initial conditions
            if typeOfInit
                params = modelExtractParam(model);
                index = paramNameRegularExpressionLookup(model, 'multi .* inverse .*');
                params(index) = log(100 +  10*rand(1,length(index)));
                model = modelExpandParam(model, params);
            end
            %Trains the model
            tic;
            model = multigpOptimise(model, display, iters);
            elapsed_time_samples(nFolds, j) = toc;
            % Save the results.
            if nargout > 2
                [CC, smse, msll] = schoolResults(model, XTest, yTest, y);
                totalErrorSMSE(nFolds, j) = mean(smse);
                totalErrorMSLL(nFolds, j) = mean(msll);
            else
                CC = schoolResults(model, XTest, yTest);
            end
            totalError(nFolds, j) = mean(CC);
        end
        fprintf('TOTAL TIME: %f \n',sum(elapsed_time_samples(:, j)));
        save(['school' upper(options.kernType(1)) options.kernType(2:end) upper(options.approx)])
    end    
else
    fprintf('NUMBER OF ACTIVE POINTS: %d \n',numActive);
    options.numActive = numActive;
    rand('twister',numberSeed);
    fprintf('ONLY ONE FOLD');
    scaleVal = zeros(1, size(yTemp, 2));
    for i = 1:size(yTemp, 2)
        maxSamples = size(yTemp{i}, 1);
        randIndex = randperm(maxSamples);
        samplesTraining = floor(0.75*maxSamples);
        indexTraining = randIndex(1:samplesTraining);
        indexTesting = randIndex(samplesTraining+1:end);
        y{i} = yTemp{i}(indexTraining);
        X{i} = XTemp{i}(indexTraining,1:q);
        options.nRepeats{i} = XTemp{i}(indexTraining,q+1);
        yTest{i} = yTemp{i}(indexTesting,:);
        XTest{i} = XTemp{i}(indexTesting,1:q);
        options.bias(i) = mean(y{i});        
        scaleVal(i) = std(y{i});        
    end
    if isScaled
        options.scale = scaleVal;
        options.beta = 10;
    else
        options.scale = ones(1, length(scaleVal));
        %options.beta = 1./((0.15*options.scale).^2);
    end
    % Creates the model
    model = multigpCreate(q, d, X, y, options);
    % Change initial conditions
    if typeOfInit
        params = modelExtractParam(model);
        index = paramNameRegularExpressionLookup(model, 'multi .* inverse .*');
        params(index) = log(100 +  10*rand(1,length(index)));
        model = modelExpandParam(model, params);
    end
    %Trains the model
    tic;
    model = multigpOptimise(model, display, iters);
    elapsed_time_samples = toc;
    % Save the results.
    if nargout > 2
        [CC, smse, msll] = schoolResults(model, XTest, yTest, y);
        totalErrorSMSE = mean(smse);
        totalErrorMSLL = mean(msll);
    else
        CC = schoolResults(model, XTest, yTest);
    end
    totalError = mean(CC);
    save(['school' upper(options.kernType(1)) options.kernType(2:end) upper(options.approx) 'OneRun'])    
end