% DEMDROSMELICMGENEHB Demo of ICM over Dros Melanogaster data

% MULTIGP

clc
clear
rand('twister',1e5);
randn('state',1e5);
%
dataSetName = 'demDrosMel';

[xTemp, yTemp, xTestTemp, yTestTemp] = mapLoadData(dataSetName);

% Use only the genes Hb, Kr, Gt, Kni
number = 3;
xTemp2 = xTemp{number};
yTemp2 = yTemp{number};
xTestTemp2 = xTestTemp{number};
yTestTemp2 = yTestTemp{number};

% Organize the dataset as we need it.
base = 1;
ut = unique(xTemp2(:,1));
us = unique(xTemp2(:,2));
XTemp = cell(1, length(ut) - base);
yTemp = cell(1, length(ut) - base);
XTempTest = cell(1, length(ut) - base);
yTempTest = cell(1, length(ut) - base);

trainSpace = 1:40;


for i=(base+1):length(ut)
   index = find(xTemp2(:,1)==ut(i));
   indexTrain = index(trainSpace);
   indexTest = index;
   indexTest(trainSpace) = [];
   XTemp{i-base} = xTemp2(indexTrain,2); 
   XTempTest{i-base} = xTemp2(indexTest,2); 
   yTemp{i-base} = yTemp2(indexTrain,1); 
   yTempTest{i-base} = yTemp2(indexTest,1); 
end

options = multigpOptions('ftc');
options.kernType = 'lmc';
options.optimiser = 'scg';
options.nlf = 1;
options.rankCorregMatrix = 7;
options.kern.nout = size(yTemp, 2);
options.kern.rankCorregMatrix = options.rankCorregMatrix;

q = 1; % Input dimension
d = size(yTemp, 2) + options.nlf;

X = cell(size(yTemp, 2)+options.nlf,1);
y = cell(size(yTemp, 2)+options.nlf,1);
XTest = cell(size(yTemp, 2)+options.nlf,1);
yTest = cell(size(yTemp, 2)+options.nlf,1);

mean_global = mean(cell2mat(yTemp'));
std_global =  std(cell2mat(yTemp'));

for j=1:options.nlf
   y{j} = [];
   X{j} = zeros(1, q);
   XTest{j} = zeros(1, q);
%    options.bias(j) = 0;
%    options.scale(j) = 1;
end
for i = 1:size(yTemp, 2)
  y{i+options.nlf} = yTemp{i};
  X{i+options.nlf} = XTemp{i};
  yTest{i+options.nlf} = yTempTest{i};
  XTest{i+options.nlf} = XTempTest{i};
%   options.bias(i+options.nlf) = mean(y{i+options.nlf});%mean_global;
%   options.scale(i+options.nlf) = std(y{i+options.nlf});%std_global;
end

% Creates the model
warning('off','multiKernParamInit:noCrossKernel')

model = multigpCreate(q, d, X, y, options);

% params = modelExtractParam(model);
% index = paramNameRegularExpressionLookup(model, 'multi .* inverse .*');
% params(index) = log(100);
% model = modelExpandParam(model, params);

display = 1;
iters = 100;

% Trains the model 
init_time = cputime;
model = multigpOptimise(model, display, iters);
elapsed_time = cputime - init_time;

[mu, varsigma] = multigpPosteriorMeanVar(model,  XTest{model.nlf+1});
[mae, mse, smse, msll, coefR] = multigpErrorMeasures(y(model.nlf+1:end), yTest(model.nlf+1:end), mu(model.nlf+1:end), ...
    varsigma(model.nlf+1:end), model.nout);


% % Save the results.
% capName = dataSetName;
% capName(1) = upper(capName(1));
% save(['dem' capName num2str(experimentNo) '.mat'], 'model');



