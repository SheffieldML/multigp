% DROS

% MULTIGP

rand('seed', 1e6);
randn('seed', 1e6);

addToolboxes(0,1)
dataSetName = 'demDrosMelInf';

% load data
[xTemp, yTemp, xTestTemp, yTestTemp] = mapLoadData(dataSetName);

options = multigpOptions('ftc');
options.optimiser = 'scg';
options.kernType = 'heat';
options.kern.nTerms = 5;
options.kern.includeIC = true;
options.nlf = 1;

% Use only the genes Hb, Kr, Gt, Kni
% xTemp = xTemp(3:6);
% yTemp = yTemp(3:6);
% xTestTemp = xTestTemp(3:6);
% yTestTemp = yTestTemp(3:6);

xTemp = xTemp(3:4);
yTemp = yTemp(3:4);
xTestTemp = xTestTemp(3:4);
yTestTemp = yTestTemp(3:4);


X = cell(length(yTemp),1);
y = cell(length(yTemp),1);
XTest = cell(length(yTemp),1);
yTest = cell(length(yTemp),1);

for i = 1:options.nlf
    XTest{i} = 0;
    X{i} = zeros(1,2);
    y{i} = [];
end
for i = 1:length(yTemp)
    y{i+options.nlf} = yTemp{i};
    X{i+options.nlf} = xTemp{i};
    yTest{i+options.nlf} = yTestTemp{i};
    XTest{i+options.nlf} = xTestTemp{i};
    options.bias(i+options.nlf) = mean(yTemp{i});
    options.scale(i+options.nlf) = std(yTemp{i});
end

% Set the input and ouput dimensions
q = 2;
d = length(yTemp) + options.nlf;

warning('off','multiKernParamInit:noCrossKernel');
warning('off','multigp:FixedParam');

% Creates the model with equal parameters for all components
%profile on
initt = cputime;
model = multigpCreate(q, d, X, y, options);
finaltt = cputime - initt; 
%profile viewer

if options.kern.includeIC || options.nlf > 1
    params = modelExtractParam(model);
    index = paramNameRegularExpressionLookup(model, '.* inverse width space\.');
    params(index) = log(1000 + 100*rand(1, length(index)));
    index = paramNameRegularExpressionLookup(model, '.* inverse width space IC\.');
    params(index) = log(1000 + 100*rand(1, length(index)));
    model = modelExpandParam(model, params);
end

display = 1;
iters = 10;

% Trains the model and counts the training time
trainingTime = cputime;
model = multigpOptimise(model, display, iters);
trainingTime = cputime - trainingTime;

[mu, varsigma] = multigpPosteriorMeanVar(model,  XTest{model.nlf+1});
[mae, mse, smse, msll] = multigpErrorMeasures(y(model.nlf+1:end), yTest(model.nlf+1:end), mu(model.nlf+1:end), ...
    varsigma(model.nlf+1:end), model.nout);

% Save the results.
dataSet = dataSetName(4:end);
dataSet(1) = lower(dataSet(1));
approx = upper(options.approx);
kernel = upper(options.kernType);
forces = num2str(options.nlf);
nTerms = num2str(options.kern.nTerms);
if options.kern.includeIC
    includeIC = '1';
else
    includeIC = '0';
end

namesOfThisFile = [dataSet approx kernel 'nTerms' nTerms 'includeIC' includeIC];

save([namesOfThisFile '.mat'], 'model', 'trainingTime', 'smse', 'msll');
