clc
clear

rand('twister',1e6);
randn('state',1e6);
%
dataSetName = 'compilerData';
experimentNo = 11;

[XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData(dataSetName);

% Configuration of options

options = multigpOptions('dtc');
options.kernType = 'ggwhite';
options.optimiser = 'scg';
options.nlf = 1;
options.initialInducingPositionMethod = 'randomData';
options.numActive = 128;
options.beta = 1e1*ones(1, size(yTemp, 2));
options.tieOptions.selectMethod = 'free';
options.nVIKs = 100;
options.fixInducing = false;
options.includeInd = false;

% Configuration of data

nSamplesTraining = 128;
displayGradchek = 1;
iters = 1000;
load('all_samples.mat');
X = cell(1,size(yTemp, 2));
y = cell(1,size(yTemp, 2));
XTest = cell(1,size(yTemp, 2));
yTest = cell(1,size(yTemp, 2));
totalError = zeros(10, size(yTemp, 2));
models = cell(size(yTemp, 2),1);

for nFolds =1:10,
    for i = 1:size(yTemp, 2)
        y{i} = yTemp{i}(all_samples(1:nSamplesTraining,nFolds),:);
        X{i} = XTemp{i}(all_samples(1:nSamplesTraining,nFolds),:);
        yTest{i} = yTemp{i};
        yTest{i}(all_samples(1:nSamplesTraining,nFolds),:) = [];
        XTest{i} = XTemp{i};
        XTest{i}(all_samples(1:nSamplesTraining,nFolds),:) = [];
        options.bias(i) = mean(y{i});
        options.scale(i) = std(y{i});
    end

    % Configuration of parameters
    q = 13;
    d = size(yTemp, 2);

    % Creates the model
    model = multigpCreate(q, d, X, y, options);
    
    %Trains the model
    init_time = cputime;
    model = multigpOptimise(model, displayGradchek, iters);
    elapsed_time = cputime - init_time;

    % Save the results.
    capName = dataSetName;
    capName(1) = upper(capName(1));
    totalError(nFolds, :) = compilerResults(dataSetName, experimentNo, XTest, yTest);
    models{nFolds} = model;
end


save(['batch' capName num2str(experimentNo) '.mat'], 'model');
