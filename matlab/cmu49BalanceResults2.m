function [rmse, nmse] = cmu49BalanceResults2(dataSetName, kernType, experimentNo)

% CMU49BALANCERESULTS2 Display the graphics (latent forces, test data and
% prediction) as well as the root mean square error (RMSE) and normalised
% mean square error (NMSE) for the CMU49 balance experiment, both using a
% convolutional model and independen GPs (if the simulation exists) for
% each output.
%
% FORMAT
% DESC Shows the results obtained with demCmu49BalanceArm.
% ARG dataSetName : string with the name of the data set used.
% ARG experimentNo : number of experiment loaded.
% ARG kernType : string with the type of kernel used (typically 'lfm' for
% demCmu49BalanceArm).
% RETURN err : RMS error for each test output signal.
%
% SEEALSO : demCmu49BalanceArm
%
% COPYRIGHT : Mauricio Alvarez and Neil D. Lawrence, 2009
%
% MODIFICATIONS : David Luengo, 2009

% MULTIGP


close all

saveFigures = 0;
capName = dataSetName;
capName(1) = upper(capName(1));
kernName = kernType;
kernName(1) = upper(kernType(1));

% Loading the IGP results (if they exist)
fileNameIgp = ['demIgp' capName num2str(experimentNo) '.mat'];
flagIgp = false;
if exist(fileNameIgp, 'file')
    flagIgp = true;
    load(fileNameIgp, 'model');
    model_igp = model;
end

% Loading the convolutional model results
load(['dem' kernName capName num2str(experimentNo) '.mat'], 'model');
% load(['dem' capName num2str(experimentNo) '.mat'], 'model');

[yTemp, void, yTestTemp, void] = lvmLoadData(dataSetName);

% Get the time index.
fps = 120/32;

% Scale the ouputs

scaleVal = sqrt(sum(var(yTemp)));
yTemp = yTemp/scaleVal;
yTestTemp = yTestTemp/scaleVal;

% Set the Options
options.type = 'multigp';
options.numModels = 1; % Number of models for test, not for training
options.compOptions = multigpOptions(model.comp{1}.approx);
options.compOptions.initialInducingPositionMethod = 'espacedInRange';
options.compOptions.optimiser = 'conjgrad';
options.compOptions.kernType = model.comp{1}.kernType;
options.compOptions.nlf = model.comp{1}.nlf;
options.compOptions.meanFunction = true;
options.compOptions.meanFunctionOptions.type = model.comp{1}.kernType;
options.compOptions.meanFunctionOptions.nlf  = options.compOptions.nlf;
options.separate = [];

% Set the inputs and outputs in the correct format
XTest = cell(1, options.numModels);
yTest = cell(1, options.numModels);
if (strcmp(options.compOptions.approx, 'ftc')) % format for full GP model
    for i = 1:options.compOptions.nlf
        XTest{1}{i} = 0;
        yTest{1}{i} = [];
    end
    for i = 1:size(yTemp, 2)
        yTest{1}{i+options.compOptions.nlf} = yTestTemp(:, i);
        XTest{1}{i+options.compOptions.nlf} = (1:size(yTestTemp, 1))'/fps;
    end
else
    for i = 1:size(yTemp, 2)
        yTest{1}{i} = yTestTemp(:, i);
        XTest{1}{i} = (1:size(yTestTemp, 1))'/fps;
    end
end

% Set the input and ouput dimensions
q = 1;
if strcmp(options.compOptions.approx, 'ftc')
    d = size(yTemp, 2) + options.compOptions.nlf;
else
    d = size(yTemp, 2);
end

testModel = multimodelCreate(q, d, {XTest{1}}, {yTest{1}}, options);
param = modelExtractParam(model);
testModel = modelExpandParam(testModel, param);
% Since each realization has its own set of parameters B, we apply the relation B_q = mean_q * D_q
for i = 1:size(yTemp, 2)
    paramInd = paramNameRegularExpressionLookup(model, ['. lfm ' num2str(i)  ' basal']);
    param(paramInd) = yTestTemp(1,i)*testModel.comp{1}.mu_D(i);
end
testModel = modelExpandParam(testModel, param);

Xt = linspace(0, 9, 100)';
[mu, varsigma] = multigpPosteriorMeanVar(testModel.comp{1}, Xt);
predVal = multigpPosteriorMeanVar(testModel.comp{1}, XTest{1}{options.compOptions.nlf+1});

if (flagIgp == true)
    rmse = zeros(2, testModel.comp{1}.nout);
    nmse = zeros(2, testModel.comp{1}.nout);
else
    rmse = zeros(1, testModel.comp{1}.nout);
    nmse = zeros(1, testModel.comp{1}.nout);
end
xlim = [min(XTest{1}{options.compOptions.nlf+1}) max(XTest{1}{options.compOptions.nlf+1})];
for k = 1:size(yTestTemp,2)+testModel.comp{1}.nlf
    figure
    hold on
    f = [(mu{k}+2*real(sqrt(varsigma{k})))*scaleVal;flipdim((mu{k}-2*real(sqrt(varsigma{k})))*scaleVal,1)];
    a = fill([Xt; flipdim(Xt,1)], f, [7 7 7]/8, 'EdgeColor', [7 7 7]/8);
    a =[ a plot(Xt, mu{k}*scaleVal,'k-')];
    set(a,   'lineWidth', 2);
    if isfield(model, 'X_u') && ~isempty(model.X_u);
        b = plot(model.X_u, -8, 'kx');
        set(b, 'linewidth', 2)
        set(b, 'markersize', 10);
    end
    if k > testModel.comp{1}.nlf
        c =plot(XTest{1}{testModel.comp{1}.nlf+1}, yTestTemp(:, k-testModel.comp{1}.nlf)*scaleVal,'k.');
        set(c,   'markersize', 20);        
        rmse(1, k-testModel.comp{1}.nlf) ...
            = sqrt(mean((predVal{k-testModel.comp{1}.nlf} - yTestTemp(:, k - testModel.comp{1}.nlf)).^2))*scaleVal;
        if flagIgp
            [mu_igp{k-testModel.comp{1}.nlf}, varsigma_igp{k-testModel.comp{1}.nlf}] ...
                = gpPosteriorMeanVar(model_igp{k-testModel.comp{1}.nlf}, Xt);
            predVal_igp{k-testModel.comp{1}.nlf} ...
                = gpPosteriorMeanVar(model_igp{k-testModel.comp{1}.nlf}, XTest{1}{options.compOptions.nlf+1});
            d = errorbar(Xt, mu_igp{k-testModel.comp{1}.nlf}*scaleVal, sqrt(varsigma_igp{k-testModel.comp{1}.nlf})*scaleVal, 'x');
            set(d, 'linewidth', 2)
            rmse(2, k-testModel.comp{1}.nlf) ...
                = sqrt(mean((predVal_igp{k-testModel.comp{1}.nlf} - yTestTemp(:, k - testModel.comp{1}.nlf )).^2))*scaleVal;
        end
    end
    set(gca, 'fontname', 'arial', 'fontsize', 15, 'xlim', xlim, 'Color', 'none')
    box on
    if saveFigures==1
        fileName = ['Toy_prediction' 'full' num2str(k)];
        print('-depsc', ['./results/' fileName]);
        saveas(gcf,['./results/' fileName],'fig');
        pos = get(gcf, 'paperposition');
        origpos = pos;
        pos(3) = pos(3)/2;
        pos(4) = pos(4)/2;
        set(gcf, 'paperposition', pos);
        lineWidth = get(gca, 'lineWidth');
        set(gca, 'lineWidth', lineWidth);
        print('-dpng', ['./results/' fileName])
        set(gca, 'lineWidth', lineWidth);
        set(gcf, 'paperposition', origpos);
    end
end

vary = var(yTestTemp)*(scaleVal^2);
if (flagIgp == true)
    vary = repmat(vary, 2, 1);
end
nmse = (rmse.^2)./vary;
ind = find(nmse > 1e10);
nmse(ind) = inf;
