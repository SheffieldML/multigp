function [rmse, nmse] = fxDataResults(dataSetName, kernType, approx, experimentNo, varargin);

% FXDATARESULTS Show the results of any given demo for the foreign exchange
% rates data set.
%
% FORMAT
% DESC Shows the results of the given demo for the FX data set.
% RETURN rmse : Root mean square error for each output.
% RETURN nmse : Normalised root mean square error for each output.
% ARG dataSetName : The data set to load.
% ARG kernType : The kernel used in the simulation ('sim' or 'simwhite').
% ARG approx : The approximation used ('ftc', 'dtc' or 'dtcvar').
% ARG experimentNo : Experiment number used in the demo.
%
% FORMAT Does the same as above but allows to specify a base directory for
% loading the result's file and flags regarding the saving and display of
% the results.
% RETURN rmse : Root mean square error for each output.
% RETURN nmse : Normalised root mean square error for each output.
% ARG dataSetName : The data set to load.
% ARG kernType : The kernel used in the simulation ('sim' or 'simwhite').
% ARG approx : The approximation used ('ftc', 'dtc' or 'dtcvar').
% ARG experimentNo : Experiment number used in the demo.
% ARG baseDirResults : Base directory for the result's file.
% ARG flags : Binary vector specifying: (1) whether the results should be
% plotted on the screen or not. (2) whether the figures should be saved.
% (3) whether a summary of the errors and likelihood of the model is shown
% at the end.
%
% COPYRIGHT : David Luengo, 2009
%
% SEEALSO : demSimDtcFxData, demSimwhiteDtcFxData, loadFxData,
% fxDataResultsPpca

% MULTIGP


if (nargin >= 5) & ~isempty(varargin{1})
    baseDirResults = varargin{1};
else
    baseDirResults = './';
end
if (nargin == 6)
    plotFigures = varargin{end}(1);
    saveFigures = varargin{end}(2);
    showError = varargin{end}(3);
else
    plotFigures = 0;
    saveFigures = 1;
    showError = 0;
end

offsetX = 1;

capName = dataSetName;
capName(1) = upper(capName(1));
kernName = kernType;
kernName(1) = upper(kernType(1));
approxName = approx;
approxName(1) = upper(approxName(1));

load([baseDirResults 'dem' kernName approxName capName 'Exp' num2str(experimentNo) '.mat'], ...
    'model', 'missingData', 'meanVal', 'scaleVal');

data = loadFxData(dataSetName);

% Setting the flags regarding the Stationarity and Normalisation of the
% kernels in the model (assuming that all of them are either stationary or
% not, normalised or not)

if isfield(model.comp{1}.kern.comp{1}.comp{model.comp{1}.nlf+1}, 'isStationary')
    isStationary = model.comp{1}.kern.comp{1}.comp{model.comp{1}.nlf+1}.isStationary;
else
    isStationary = false;
end
if isfield(model.comp{1}.kern.comp{1}.comp{model.comp{1}.nlf+1}, 'isNormalised')
    isNormalised = model.comp{1}.kern.comp{1}.comp{model.comp{1}.nlf+1}.isNormalised;
else
    isNormalised = false;
end

% Set the Options
options.type = 'multigp';
options.numModels = 1; % Number of models for test, not for training
options.compOptions = multigpOptions(model.comp{1}.approx);
if ~strcmp(model.comp{1}.approx, 'ftc')
    options.compOptions.numActive = model.comp{1}.k;
end
options.compOptions.initialInducingPositionMethod = 'espacedInRange';
options.compOptions.kernType = model.comp{1}.kernType;
options.compOptions.nlf = model.comp{1}.nlf;
if isfield(model.comp{1}, 'typeLf')
    options.compOptions.typeLf = model.comp{1}.typeLf;
end
options.compOptions.tieOptions.selectMethod = 'typeLf';
options.separate = [];
if isfield(model.comp{1}, 'meanFunction')
    options.compOptions.meanFunction = true;
    options.compOptions.meanFunctionOptions.type = options.compOptions.kernType;
    options.compOptions.meanFunctionOptions.nlf  = options.compOptions.nlf;
end

% Preparing the test data that are going to be plotted later
XTestPlot = cell(1, options.numModels);
yTestPlot = cell(1, options.numModels);
count = offsetX + (0:size(data.XTest, 1)-1)';
if (strcmp(options.compOptions.approx, 'ftc')) % format for full GP model
    for i = 1:options.compOptions.nlf
        XTestPlot{1}{i} = 0;
        yTestPlot{1}{i} = 0;
    end
    for i = 1:model.comp{1}.nout
        ind = find(data.yTest{i} ~= 0.0);
        yTestPlot{1}{i+options.compOptions.nlf} = data.yTest{i}(ind);
        XTestPlot{1}{i+options.compOptions.nlf} = count(ind);
    end
else
    for i = 1:model.comp{1}.nout
        ind = find(data.yTest{i} ~= 0.0);
        yTestPlot{1}{i} = data.yTest{i}(ind);
        XTestPlot{1}{i} = count(ind);
        ind2 = find(data.yTest{i} == 0.0);
        yTestPlot2{1}{i} = data.yTest{i};
        yTestPlot2{1}{i}(ind2) = NaN;
        XTestPlot2{1}{i} = count;
    end
end

% Set the inputs and outputs in the correct format
XTest = cell(1, options.numModels);
yTest = cell(1, options.numModels);
if (strcmp(options.compOptions.approx, 'ftc')) % format for full GP model
    for i = 1:options.compOptions.nlf
        XTest{1}{i} = 0;
        yTest{1}{i} = 0;
    end
    for i = 1:model.comp{1}.nout
        yTest{1}{i+options.compOptions.nlf} = data.yTest{i};
        XTest{1}{i+options.compOptions.nlf} = offsetX + (0:size(data.yTest{i}, 1)-1)';
    end
else
    for i = 1:model.comp{1}.nout
        if ~isempty(missingData{i})
            data.yTest{i}(missingData{i}) = 0.0;
        end
        ind = find(data.yTest{i} ~= 0.0);
        yTest{1}{i} = (data.yTest{i}(ind) - meanVal(i))/scaleVal(i);
        XTest{1}{i} = count(ind);
    end
end

% Set the input and ouput dimensions and create the model

q = 1;
d = model.outputDim;
testModel = multimodelCreate(q, d, XTest, yTest, options);

% Forcing the kernels to be stationary (if that was the case in the demo)

if (isStationary == true)
    for i=1:model.comp{1}.nlf
        for j=1:model.comp{1}.nout
            testModel.comp{1}.kern.comp{i}.comp{testModel.comp{1}.nlf+j}.isStationary = true;
        end
    end
end

% Forcing the SIM and SIM-WHITE kernels to be normalised (if that was the
% case in the demo)

if isNormalised
    for i=1:model.comp{1}.nlf
        for j=1:model.comp{1}.nlf+model.comp{1}.nout
            testModel.comp{1}.kern.comp{i}.comp{j}.isNormalised = true;
        end
    end
end

% Setting the parameters of the test model to the values learnt for model

param = modelExtractParam(model);
testModel = modelExpandParam(testModel, param);
testModel.comp{1}.X_u = model.comp{1}.X_u;

% Prediction

if strcmp(options.compOptions.approx, 'ftc')
    Xt = linspace(count(1), count(end), 500)';
    [mu, varsigma] = multigpPosteriorMeanVar(testModel.comp{1}, Xt);
    [predVal, varsigma2] = multigpPosteriorMeanVar(testModel.comp{1}, XTestPlot{1});
elseif strcmp(options.compOptions.approx, 'dtc') || strcmp(options.compOptions.approx, 'dtcvar')
    Xt = linspace(count(1), count(end), 500)';
    [mu, varsigma] = multigpPosteriorMeanVar(testModel.comp{1}, Xt);
    [predVal, varsigma2] = multigpPosteriorMeanVar(testModel.comp{1}, XTestPlot{1});
else
    % TBD
    error('Approximations different to FTC and DTC are currently not supported');
end

% Undoing the mean value substraction and scaling

for k = testModel.comp{1}.nlf+1:length(predVal)
    mu{k} = mu{k}*scaleVal(k-testModel.comp{1}.nlf) + meanVal(k-testModel.comp{1}.nlf);
    varsigma{k} = varsigma{k}*(scaleVal(k-testModel.comp{1}.nlf)^2);
    predVal{k} = predVal{k}*scaleVal(k-testModel.comp{1}.nlf) + meanVal(k-testModel.comp{1}.nlf);
    varsigma2{k} = varsigma2{k}*(scaleVal(k-testModel.comp{1}.nlf)^2);
end

if exist(['demPpca' dataSetName 'Nlf' num2str(options.compOptions.nlf) '.mat'], 'file')
    eval(['load demPpca' dataSetName 'Nlf' num2str(options.compOptions.nlf)]);
end

% Plotting the results

close all
xlim = [min(XTest{1}{end}) max(XTest{1}{end})];
for k = 1:length(predVal)
    if plotFigures
        figure
        hold on
        if strcmp(options.compOptions.approx, 'ftc')
            f = [(mu{k}+2*real(sqrt(varsigma{k}))); flipdim((mu{k}-2*real(sqrt(varsigma{k}))),1)];
            a = fill([Xt; flipdim(Xt,1)], f, [7 7 7]/8, 'EdgeColor', [7 7 7]/8);
            a = [a plot(Xt, mu{k},'k-')];
            set(a, 'lineWidth', 2);
            % Plotting the true data
            if k > testModel.comp{1}.nlf
                c = plot(XTest{1}{k}, yTest{1}{k}*scaleVal(k)+meanVal(k), 'k.');
                d = plot(XTestPlot{1}{k}, yTestPlot{1}{k},'k--');
                set(c, 'markerSize', 20);
                set(d, 'lineWidth', 2)
            end
        elseif strcmp(options.compOptions.approx, 'dtc') || strcmp(options.compOptions.approx, 'dtcvar')
            if k <= testModel.comp{1}.nlf
                f = [(mu{k}+2*real(sqrt(varsigma{k}))); flipdim((mu{k}-2*real(sqrt(varsigma{k}))),1)];
                a = fill([Xt; flipdim(Xt,1)], f, [7 7 7]/8, 'EdgeColor', [7 7 7]/8);
                a = [a plot(Xt, mu{k},'k-')];
                set(a, 'lineWidth', 2);
            else
                f = [(mu{k}+2*real(sqrt(varsigma{k}))); flipdim((mu{k}-2*real(sqrt(varsigma{k}))),1)];
                a = fill([Xt; flipdim(Xt, 1)], f, [7 7 7]/8, 'EdgeColor', [7 7 7]/8);
                a =[a plot(Xt, mu{k}, 'k-')];
                set(a, 'lineWidth', 2);
                % Plotting the inducing points
                b = plot(model.comp{1}.X_u, min(f)-0.1*(max(f)-min(f)), 'kx');
                set(b, 'lineWidth', 2)
                set(b, 'markerSize', 10);
                % Plotting the true data
                c = plot(XTest{1}{k-testModel.comp{1}.nlf}, ...
                    yTest{1}{k-testModel.comp{1}.nlf}*scaleVal(k-testModel.comp{1}.nlf)+meanVal(k-testModel.comp{1}.nlf), 'k.');
                d = plot(XTestPlot2{1}{k-testModel.comp{1}.nlf}, yTestPlot2{1}{k-testModel.comp{1}.nlf}, 'k--');
%                 e = plot(count, yEst(:, k-testModel.comp{1}.nlf), 'r');
                set(c, 'markerSize', 20);
                set(d, 'lineWidth', 2)
            end
        else
            % TBD
        end
        set(gca, 'fontname', 'arial', 'fontsize', 15, 'xlim', xlim, 'Color', 'none')
        box on
        if saveFigures==1
            fileName = ['dem' kernName '_' capName '_Exp' num2str(experimentNo) '_Fig' num2str(k)];
            print('-depsc', [baseDirResults fileName]);
            saveas(gcf,[baseDirResults fileName],'fig');
            pos = get(gcf, 'paperposition');
            origpos = pos;
            pos(3) = pos(3)/2;
            pos(4) = pos(4)/2;
            set(gcf, 'paperposition', pos);
            lineWidth = get(gca, 'lineWidth');
            set(gca, 'lineWidth', lineWidth);
            print('-dpng', [baseDirResults fileName])
            set(gca, 'lineWidth', lineWidth);
            set(gcf, 'paperposition', origpos);
        end
    end
    % Compute Error over test inputs
    if k > testModel.comp{1}.nlf 
        if strcmp(options.compOptions.approx, 'ftc')
            rmse(k-testModel.comp{1}.nlf) = std(predVal{k} - yTestPlot{1}{k});
        elseif strcmp(options.compOptions.approx, 'dtc') || strcmp(options.compOptions.approx, 'dtcvar')
            rmse(k-testModel.comp{1}.nlf) = std(predVal{k} - yTestPlot{1}{k-testModel.comp{1}.nlf});
        else
            % TBD
        end
    end
end

vary = cellfun(@var, yTestPlot{1});
nmse = (rmse.^2)./vary;

numMissingData = 0;
for i=1:length(missingData)
    if exist('yEst', 'var')
        nmsePpca(i) = var(yEst(XTestPlot{1}{i}, i)-yTestPlot{1}{i});
    end
    if isempty(missingData{i})
        nmseTest(i) = 0;
        nmseTrain(i) = nmse(i);
        nmseTestPpca(i) = 0;
        nmseTrainPpca(i) = nmsePpca(i) / vary(i);
    else
        ind1 = missingData{i}(1) - 1;
        ind2 = missingData{i}(end) + 1;
        yTrain = [data.y{i}(1:ind1); data.y{i}(ind2:end)];
        yPredTrain = [predVal{model.comp{1}.nlf+i}(1:ind1); predVal{model.comp{1}.nlf+i}(ind2:end)];
        varyTrain = var(yTrain);
        nmseTrain(i) = var(yPredTrain - yTrain) / varyTrain;
        varyTest = var(data.y{i}(missingData{i}));
        nmseTest(i) = var(data.y{i}(missingData{i})-predVal{model.comp{1}.nlf+i}(missingData{i}))/varyTest;
        if exist('yEst', 'var')
            yPredTrainPpca = [yEst(1:ind1, i); yEst(ind2:end, i)];
            nmseTrainPpca(i) = var(yPredTrainPpca - yTrain) / varyTrain;
            nmseTestPpca(i) = var(data.y{i}(missingData{i})-yEst(missingData{i}, i))./varyTest;
        end
        numMissingData = numMissingData + 1;
    end
end
if exist('yEst', 'var')
    nmsePpca = nmsePpca./vary;
end

if showError
    disp(['Root Mean Square Error (RMSE)                             : ' num2str(rmse)]);
    disp(['Standardised MSE (SMSE) - Test and training data          : ' num2str(nmse)]);
    disp(['Standardised MSE (SMSE) - Training data only              : ' num2str(nmseTrain)]);
    disp(['Standardised MSE (SMSE) - Test data only                  : ' num2str(nmseTest)]);
    if exist('yEst', 'var')
        disp(['Standardised MSE (SMSE) for PPCA - Test and training data : ' num2str(nmsePpca)]);
        disp(['Standardised MSE (SMSE) for PPCA - Training data only     : ' num2str(nmseTrainPpca)]);
        disp(['Standardised MSE (SMSE) for PPCA - Test data only         : ' num2str(nmseTestPpca)]);
    end
    disp(['Model -log-likelihood                                     : ' num2str(-modelLogLikelihood(model))]);
    disp(['Average SMSE (test and training data)                     : ' num2str(100*mean(nmse)) ' %']);
    disp(['Average SMSE (training data only)                         : ' num2str(100*mean(nmseTrain)) ' %']);
    disp(['Average SMSE (test data only)                             : ' num2str(100*sum(nmseTest)/numMissingData) ' %']);
    if exist('yEst', 'var')
        disp(['Average SMSE for PPCA (test and training data)            : ' num2str(100*mean(nmsePpca)) ' %']);
        disp(['Average SMSE for PPCA (training data only)                : ' num2str(100*mean(nmseTrainPpca)) ' %']);
        disp(['Average SMSE for PPCA (test data only)                    : ' num2str(100*sum(nmseTestPpca)/numMissingData) ' %']);
    end
end

return