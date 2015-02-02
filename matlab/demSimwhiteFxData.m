% DEMSIMWHITEFXDATA Demonstrate latent force model on exchange rates data.

% MULTIGP


close all
clear
clc

rand('seed', 1e6);
randn('seed', 1e6);

dataSetName = 'ShortFxData2009-6outputs';
% dataSetName = 'ShortFxData2008-2009-6outputs';
% dataSetName = 'LongFxData2008-2009-6outputs';
% dataSetName = 'LongFxData2007-2008';
experimentNo = 2;

% load data

if exist('numSamplesTrainingSet', 'var')
    data = loadFxData(dataSetName, numSamplesTrainingSet);
else
    data = loadFxData(dataSetName);
end

% Scale the ouputs

meanVal = zeros(size(data.y));
scaleVal = ones(size(data.y));

% Set the Options 
options.type = 'multigp';
options.numModels = 1;
options.compOptions = multigpOptions('ftc');
% options.compOptions.includeNoise = false;
options.compOptions.optimiser = 'conjgrad';
% options.compOptions.optimiser = 'scg';
options.compOptions.kernType = 'simwhite';
options.compOptions.nlf = 2;
options.compOptions.meanFunction = true;
options.compOptions.meanFunctionOptions.type = options.compOptions.kernType;
options.compOptions.meanFunctionOptions.nlf  = options.compOptions.nlf;
options.separate = [];

% Set the missing data for the demo

missingData = cell(1, size(data.y, 2));

% Set the inputs and outputs in the correct format
X = cell(1, options.numModels);
y = cell(1, options.numModels);
XTest = cell(1, options.numModels);
yTest = cell(1, options.numModels);
if (options.numModels == 1)
    for i = 1:options.compOptions.nlf
        X{1}{i} = 0;
        y{1}{i} = 0;
    end
    count = (1:size(data.X, 1))';
    for i = 1:size(data.y, 2)
        if ~isempty(missingData{i})
            data.y{i}(missingData{i}) = 0.0;
        end
        ind = find(data.y{i} ~= 0.0);
        scaleVal(i) = std(data.y{i}(ind));
%         meanVal(i) = mean(data.y{i}(ind));
%         meanVal(i) = data.y{i}(ind(1));
        y{1}{i+options.compOptions.nlf} = (data.y{i} - meanVal(i))/scaleVal(i);
        X{1}{i+options.compOptions.nlf} = count(ind);
    end
else
    error('The current demo admits only one model');
end

% Set the input and ouput dimensions
q = 1;
d = size(y{1}, 2);

% Creates the model
model = multimodelCreate(q, d, X, y, options);
model.optimiser = model.comp{1}.optimiser;

% Make a better initilization of the basal rate parameter
params = modelExtractParam(model);
for i = 1:size(y{1}, 2)-options.compOptions.nlf
    paramInd = paramNameRegularExpressionLookup(model, ['. simwhite ' num2str(i)  ' basal']);
    params(paramInd) = log(mean(y{1}{i+options.compOptions.nlf}));
end
model = modelExpandParam(model, params);

display = 1;
iters = 1000;

% Trains the model and counts the training time
trainingTime = cputime;
model = modelOptimise(model, [], [], display, iters);
trainingTime = cputime - trainingTime;

% Save and plot the results.
capName = dataSetName;
capName(1) = upper(capName(1));
save(['demSimwhite' capName 'Exp' num2str(experimentNo) '.mat'], ...
    'model', 'missingData', 'meanVal', 'scaleVal');

if exist('numSamplesTrainingSet', 'var')
    [rmse, nmse] = fxDataResults(dataSetName, 'simwhite', experimentNo, numSamplesTrainingSet, [1 1 1]);
else
    [rmse, nmse] = fxDataResults(dataSetName, 'simwhite', experimentNo, [], [1 1 1]);
end
