function [rmse, nmse] = fxDataResultsPpca(dataSetName, nlf, varargin);

% FXDATARESULTSPPCA Show the results of the PPCA demo for the foreign
% exchange rates data set.
%
% FORMAT
% DESC Shows the results of the PPCA demo for the FX data set.
% RETURN rmse : Root mean square error for each output.
% RETURN nmse : Normalised root mean square error for each output.
% ARG dataSetName : The data set to load.
% ARG nlf : Dimension of the latent space, i.e. number of principal
% components used for PPCA.
%
% FORMAT Does the same as above but allows to specify a base directory for
% loading the result's file and flags regarding the saving and display of
% the results.
% RETURN rmse : Root mean square error for each output.
% RETURN nmse : Normalised root mean square error for each output.
% ARG dataSetName : The data set to load.
% ARG nlf : Dimension of the latent space, i.e. number of principal
% components used for PPCA.
% ARG flags : Binary vector specifying: (1) whether the results should be
% plotted on the screen or not. (2) whether the figures should be saved.
% (3) whether a summary of the errors and likelihood of the model is shown
% at the end.
%
% COPYRIGHT : David Luengo, 2009
%
% SEEALSO : fxDataResults

% MULTIGP


if (nargin == 3)
    plotFigures = varargin{1}(1);
    saveFigures = varargin{1}(2);
    showError = varargin{1}(3);
else
    plotFigures = 0;
    saveFigures = 1;
    showError = 0;
end

baseDirResults = './';

offsetX = 1;

data = loadFxData(dataSetName);

% Preparing the test data that are going to be plotted later

nout = length(data.y);
XTestPlot = cell(1, 1);
yTestPlot = cell(1, 1);
dataSet = cell(1, 1);
count = offsetX + (0:size(data.XTest, 1)-1)';
for i = 1:nout
    ind = find(data.yTest{i} ~= 0.0);
    dataSet{1}{i} = ind;
    yTestPlot{1}{i} = data.yTest{i}(ind);
    XTestPlot{1}{i} = count(ind);
    ind2 = find(data.yTest{i} == 0.0);
    yTestPlot2{1}{i} = data.yTest{i};
    yTestPlot2{1}{i}(ind2) = NaN;
    XTestPlot2{1}{i} = count;
end

% Set the missing data for the demo

missingData = cell(1, size(data.y, 2));
missingData{4} = 50:100;
missingData{6} = 100:150;
missingData{9} = 150:200;

% Set the inputs and outputs in the correct format for the model

XTest = cell(1, 1);
yTest = cell(1, 1);
yTestMissing = cell(1, 1);
trainingSet = cell(1, 1);
testSet = cell(1, 1);
for i = 1:nout
    if ~isempty(missingData{i})
        yTestMissing{1}{i} = data.yTest{i}(missingData{i});
        data.yTest{i}(missingData{i}) = 0.0;
    else
        yTestMissing{1}{i} = [];
    end
    ind = find(data.yTest{i} ~= 0.0);
    trainingSet{1}{i} = ind;
    meanVal(i) = 0;
    scaleVal(i) = 1;
    yTest{1}{i} = (data.yTest{i}(ind) - meanVal(i))/scaleVal(i);
    XTest{1}{i} = count(ind);
    ind2 = find(data.yTest{i} == 0.0);
    testSet{1}{i} = missingData{i};
end

% Loading the PPCA results (if the simulation exists) or performing the
% simulation (if it doesn't)

try
    load(strcat('demPpca', dataSetName, 'Nlf', num2str(nlf)));
catch
    [yEst, z, W] = demPpcaFxData(dataSetName, nlf);
end

% Obtaining a version without the missing data for evaluating the error

yEst2 = cell(1, 1);
for i=1:nout
    ind = find(~isnan(yTestPlot2{1}{i}));
    yEst2{1}{i} = yEst(ind, i);
end

% Plotting the results

close all
xlim = [min(XTest{1}{end}) max(XTest{1}{end})];
for k = 1:length(data.y)
    if plotFigures
        figure
        hold on
        % Plotting the true data
        c = plot(XTest{1}{k}, yTest{1}{k}*scaleVal(k)+meanVal(k), 'k.');
        d = plot(XTestPlot2{1}{k}, yTestPlot2{1}{k}, 'k--');
        e = plot(count, yEst(:, k), 'k');
        set(c, 'markerSize', 20);
        set(d, 'lineWidth', 2);
        set(e, 'lineWidth', 2);
        set(gca, 'fontname', 'arial', 'fontsize', 15, 'xlim', xlim, 'Color', 'none')
        box on
        if saveFigures==1
            fileName = ['dem_' dataSetName '_Ppca_Nlf' num2str(nlf) '_Fig' num2str(k)];
            print('-depsc', [baseDirResults 'results/' fileName]);
            saveas(gcf,[baseDirResults 'results/' fileName],'fig');
            pos = get(gcf, 'paperposition');
            origpos = pos;
            pos(3) = pos(3)/2;
            pos(4) = pos(4)/2;
            set(gcf, 'paperposition', pos);
            lineWidth = get(gca, 'lineWidth');
            set(gca, 'lineWidth', lineWidth);
            print('-dpng', [baseDirResults 'results/' fileName])
            set(gca, 'lineWidth', lineWidth);
            set(gcf, 'paperposition', origpos);
        end
    end
    % Compute Error over test inputs
    rmse(k) = std(mean(yEst2{1}{k} - yTestPlot{1}{k}));
end

vary = cellfun(@var, yTestPlot{1});
nmse = (rmse.^2)./vary;

numMissingData = 0;
for i=1:length(missingData)
    nmsePpca(i) = var(yEst(dataSet{1}{i}, i)-yTestPlot{1}{i}) / vary(i);
    nmseTrainPpca(i) = var(yEst(trainingSet{1}{i}, i) - yTest{1}{i}) / var(yTest{1}{i});
    if ~isempty(testSet{1}{i})
        nmseTestPpca(i) = var(yEst(testSet{1}{i}, i) - yTestMissing{1}{i}) / var(yTestMissing{1}{i});
        numMissingData = numMissingData + 1;
    else
        nmseTestPpca(i) = 0;
    end
end

if showError
    disp(['Standardised MSE (SMSE) for PPCA - Test and training data : ' num2str(nmsePpca)]);
    disp(['Standardised MSE (SMSE) for PPCA - Training data only     : ' num2str(nmseTrainPpca)]);
    disp(['Standardised MSE (SMSE) for PPCA - Test data only         : ' num2str(nmseTestPpca)]);
    disp(['Average SMSE for PPCA (test and training data)            : ' num2str(100*mean(nmsePpca)) ' %']);
    disp(['Average SMSE for PPCA (training data only)                : ' num2str(100*mean(nmseTrainPpca)) ' %']);
    disp(['Average SMSE for PPCA (test data only)                    : ' num2str(100*sum(nmseTestPpca)/numMissingData) ' %']);
end
