function [totalError, elapsed_time_train, elapsed_time_test] = ...
    compilerFullCore(dataSetName, options, numTraining, ...
    display, iters, totFolds)

% COMPILERFULLCORE description. 

% MULTIGP

[XTemp, yTemp] = mapLoadData(dataSetName);

load('all_samples.mat');

X = cell(1,size(yTemp, 2)+options.nlf);
y = cell(1,size(yTemp, 2)+options.nlf);
XTest = cell(1,size(yTemp, 2));
yTest = cell(1,size(yTemp, 2));
q = size(XTemp{1}, 2);
d = size(yTemp, 2)+options.nlf;

totalError = cell(length(numTraining),1);
elapsed_time_train = cell(length(numTraining),1);
elapsed_time_test = cell(length(numTraining),1);

for ns = 1:length(numTraining)
    nSamplesTraining = numTraining(ns);
    totalErrorSamples = zeros(totFolds, size(yTemp, 2));
    elapsed_time_samples = zeros(totFolds,1);
    elapsed_time_samples_test = zeros(totFolds, 1);
    fprintf('NUMBER OF TRAINING POINTS: %d \n',numTraining(ns));
    rand('twister',10^7);
    randn('state',10^7);            
    for nFolds =1:totFolds
        fprintf('FOLD NUMBER: %d \n', nFolds);
        for i=1:options.nlf
            y{i} = [];
            X{i} = ones(1,q);
            options.bias(i) = 0;
            options.scale(i) = 1;
        end
        for i = 1:size(yTemp, 2)
            y{i+options.nlf} = yTemp{i}(all_samples(1:nSamplesTraining,nFolds),:);
            X{i+options.nlf} = XTemp{i}(all_samples(1:nSamplesTraining,nFolds),:);
            yTest{i} = yTemp{i};
            yTest{i}(all_samples(1:nSamplesTraining,nFolds),:) = [];
            XTest{i} = XTemp{i};
            XTest{i}(all_samples(1:nSamplesTraining,nFolds),:) = [];
            options.bias(i+options.nlf) = mean(y{i+options.nlf});
            options.scale(i+options.nlf) = std(y{i+options.nlf});
        end
        % Creates the model
        model = multigpCreate(q, d, X, y, options);
        % Change initial parameters of noise
        params = modelExtractParam(model);
        index = paramNameRegularExpressionLookup(model, '.* white .* variance');
        params(index(model.nlf+1:end)) = log(options.noise);
        model = modelExpandParam(model, params);
        %Trains the model
        tic;
        model = multigpOptimise(model, display, iters);
        elapsed_time_samples(nFolds) = toc;
        % Save the results.
        [totalErrorSamples(nFolds,:), ...
            elapsed_time_samples_test(nFolds)] = compilerResults(model, XTest, yTest);
    end
    fprintf('TOTAL TIME: %f \n',sum(elapsed_time_samples_test) ...
        + sum(elapsed_time_samples));
    totalError{ns} = totalErrorSamples;
    elapsed_time_train{ns} = elapsed_time_samples;
    elapsed_time_test{ns} = elapsed_time_samples_test;
    save(['compiler' upper(options.kernType(1)) options.kernType(2:end) upper(options.approx)]) 
end

