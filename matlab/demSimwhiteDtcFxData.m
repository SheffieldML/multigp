% DEMSIMWHITEDTCFXDATA Demonstrate latent force model on exchange rates
% data using the SIM-WHITE model.
%
% COPYRIGHT : David Luengo, 2009
%
% SEEALSO : demSimDtcFxData, loadFxData, fxDataResults

% MULTIGP


close all
clear
clc

rand('seed', 1e6);
randn('seed', 1e6);

% Data set used

% dataSetName = 'ShortFxData2009-6outputs';
dataSetName = 'LongFxData2007-13outputs';
% dataSetName = 'ShortFxData2008-2009-6outputs';
% dataSetName = 'LongFxData2008-2009-6outputs';
% dataSetName = 'LongFxData2007-2008';

% baseDirResults = './';
baseDirResults = './results/October2009/';

% Experiment and parameters of the demo

numLatentFunc = 4;
numInducingPoints = 100;
approx = 'dtcvar';
% isStationary = false;
isStationary = true;
% isNormalised = false;
isNormalised = true;
useMeanFunctions = false;
% useMeanFunctions = true;
optionMeanVal = 'sampleMean';
% optionMeanVal = 'initVal';
% useIndGPInit = true;
useIndGPInit = false;
offsetX = 1;

display = 1;
iters = 2000;

experimentNo = str2num(strcat(num2str(iters), num2str(numInducingPoints), ...
    '0', num2str(numLatentFunc)));

optionDisplayResults = true;
% optionDisplayResults = false;
optionPlotFigures = true;
% optionPlotFigures = false;
% optionSaveFigures = true;
optionSaveFigures = true;
optionShowError = true;
% optionShowError = false;

% Load data

data = loadFxData(dataSetName);

% Set the Options 
options.type = 'multigp';
options.numModels = 1;
options.compOptions = multigpOptions(approx);
options.compOptions.numActive = numInducingPoints;
options.compOptions.initialInducingPositionMethod = 'espacedInRange';
% options.compOptions.optimiser = 'conjgrad';
options.compOptions.optimiser = 'scg';
options.compOptions.kernType = 'simwhite';
options.compOptions.nlf = numLatentFunc;
options.compOptions.tieOptions.selectMethod = 'typeLf';
options.separate = [];
if useMeanFunctions == true
    options.compOptions.meanFunction = true;
    options.compOptions.meanFunctionOptions.type = options.compOptions.kernType;
    options.compOptions.meanFunctionOptions.nlf  = options.compOptions.nlf;
end

% Set the missing data for the demo

missingData = cell(1, size(data.y, 2));
missingData{4} = 50:100;
missingData{6} = 100:150;
missingData{9} = 150:200;

% Set the inputs and outputs in the correct format
X = cell(1, options.numModels);
y = cell(1, options.numModels);
XTest = cell(1, options.numModels);
yTest = cell(1, options.numModels);
meanVal = zeros(size(data.y));
scaleVal = ones(size(data.y));
if (options.numModels == 1)
    count = offsetX + (0:size(data.X, 1)-1)';
    for i = 1:size(data.y, 2)
        if ~isempty(missingData{i})
            data.y{i}(missingData{i}) = 0.0;
        end
        ind = find(data.y{i} ~= 0.0);
        scaleVal(i) = std(data.y{i}(ind));
        if strcmp(optionMeanVal, 'sampleMean')
            meanVal(i) = mean(data.y{i}(ind));
        elseif strcmp(optionMeanVal, 'initVal')
            meanVal(i) = data.y{i}(ind(1));
        else
            error('Mean value option not defined');
        end
        y{1}{i} = (data.y{i}(ind) - meanVal(i))/scaleVal(i);
        X{1}{i} = count(ind);
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
model.comp{1}.X_u = linspace(count(1), count(end), model.comp{1}.k(1))';

% Forcing the kernels to be stationary

if isStationary
    for i=1:model.comp{1}.nlf
        for j=1:model.comp{1}.nout
            model.comp{1}.kern.comp{i}.comp{model.comp{1}.nlf+j}.isStationary = true;
        end
    end
end

% Make a better initilization of the basal rate parameter when the mean
% functions are being used

if (useMeanFunctions == true)
    params = modelExtractParam(model);
    for i = 1:size(y{1}, 2)
        paramInd = paramNameRegularExpressionLookup(model, ['. simwhite ' num2str(i)  ' basal']);
        params(paramInd) = log(mean(y{1}{i}));
    end
    model = modelExpandParam(model, params);
end

% Initialisation of the model's parameters using a set of independent GPs
% (if desired). Otherwise the parameters associated with the inverse widths
% of the inducing kernels are simply set to appropriate different initial
% values (defined by the user) for symmetry breaking purposes.

if (useIndGPInit == true)
    % Try to load the file with the appropriate parameters according to the
    % experiment's number. Otherwise generate it.
    fileName = ['paramsInitSimwhiteExp' num2str(experimentNo)];
    try
        load(fileName);
    catch
        paramsInit = ...
            initDemSimDtcFxData(X, y, count, options.compOptions, [isNormalised isStationary]);
        save(fileName);
    end
    % Set the number of smooth and non-smooth latent forces.
    if isfield(model.comp{1}, 'typeLf')
        NumLfType1 = model.comp{1}.typeLf(1);
        NumLfType2 = model.comp{1}.typeLf(2);
    else
        NumLfType1 = model.comp{1}.nlf;
        NumLfType2 = 0;
    end
    % Get the current parameters of the model and insert the new ones.
    params = modelExtractParam(model);
    switch options.compOptions.tieOptions.selectMethod
        case 'typeLf'
            for i = 1:NumLfType1
                paramInd = paramNameRegularExpressionLookup(model, ['multi ' num2str(i) ...
                    ' .* inverse width']);
                params(paramInd) = paramsInit(1, 1, i);
                for j = 1:model.comp{1}.nout
                    paramInd = paramNameRegularExpressionLookup(model, ['multi ' num2str(i) ...
                        model.comp{1}.kernType ' ' num2str(j) ' decay']);
                    params(paramInd) = log(mean(exp(paramsInit(j, 3, :))));
                    if strcmp(model.comp{1}.kernType, 'sim')
                        paramInd = paramNameRegularExpressionLookup(model, ['multi ' num2str(i) ...
                            ' sim ' num2str(j) ' variance']);
                        params(paramInd) = paramsInit(j, 4, i);
                    elseif strcmp(model.comp{1}.kernType, 'simwhite')
                        paramInd = paramNameRegularExpressionLookup(model, ['multi ' num2str(i) ...
                            ' simwhite ' num2str(j) ' sensitivity']);
                        params(paramInd) = paramsInit(j, 4, i);
                    else
                        error('Kernel type not supported yet');
                    end
                end
            end
            for i = NumLfType1+1:NumLfType2
                paramInd = paramNameRegularExpressionLookup(model, ['multi ' num2str(i) ...
                    ' .* inverse width']);
                params(paramInd) = paramsInit2(1, 1, i-NumLfType1);
                for j = 1:model.comp{1}.nout
                    paramInd = paramNameRegularExpressionLookup(model, ['multi ' num2str(i) ...
                        ' simwhite ' num2str(j) ' decay']);
                    params(paramInd) = log(mean(exp(paramsInit2(j, 3, :))));
                    paramInd = paramNameRegularExpressionLookup(model, ['multi ' num2str(i) ...
                        ' simwhite ' num2str(j) ' sensitivity']);
                    params(paramInd) = paramsInit2(j, 4, i-NumLfType1);
                end
            end
    end
    model = modelExpandParam(model, params);
else
    % No initialisation using independent GPs.
    invWidth0 = 1./([1 1.5 2 2.5 3 3.5 4].^2);
    params = modelExtractParam(model);
    for i = 1:options.compOptions.nlf
        paramInd = paramNameRegularExpressionLookup(model, ['multi ' num2str(i) ...
            ' .* inverse width']);
        params(paramInd) = expTransform(invWidth0(i), 'xtoa');
    end
    model = modelExpandParam(model, params);
end

% Trains the model and counts the training time

tic;
model = modelOptimise(model, [], [], display, iters);
trainingTime = toc;

% Save and plot the results.
capName = dataSetName;
capName(1) = upper(capName(1));
save([baseDirResults 'demSimwhite' upper(approx(1)) approx(2:end) capName ...
    'Exp' num2str(experimentNo) '.mat'], 'model', 'missingData', ...
    'meanVal', 'scaleVal', 'trainingTime');

if optionDisplayResults == true
    [rmse, nmse] = fxDataResults(dataSetName, model.comp{1}.kernType, model.comp{1}.approx, ...
        experimentNo, baseDirResults, [optionPlotFigures optionSaveFigures optionShowError]);
end
