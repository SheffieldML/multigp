% DEMTOY1DLMCFTC1 Demo of full multi output GP with LMC kernel.
% FORMAT
% DESC Demo of Full Multi Output Gaussian Process with LMC kernel 

% MULTIGP

clc
clear
rand('twister',1e5);
randn('state',1e5);
%
%dataSetName = 'ggToyTrainTest';
dataSetName = 'ggToyMissing';
experimentNo = 11;

[XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData(dataSetName);

options = multigpOptions('ftc');
options.kernType = 'lmc';
options.optimiser = 'scg';
options.nlf = 1;
options.rankCorregMatrix = 2;
options.kern.nout = size(yTemp, 2);
options.kern.rankCorregMatrix = options.rankCorregMatrix;


q = 1; % Input dimension
d = size(yTemp, 2) + options.nlf;

X = cell(size(yTemp, 2)+options.nlf,1);
y = cell(size(yTemp, 2)+options.nlf,1);

for j=1:options.nlf
   y{j} = [];
   X{j} = zeros(1, q);  
end
for i = 1:size(yTemp, 2)
  y{i+options.nlf} = yTemp{i};
  X{i+options.nlf} = XTemp{i};
end

% Creates the model
warning('off','multiKernParamInit:noCrossKernel')

model = multigpCreate(q, d, X, y, options);

params = modelExtractParam(model);
index = paramNameRegularExpressionLookup(model, 'multi .* inverse .*');
params(index) = log(100);
model = modelExpandParam(model, params);

display = 2;
iters = 200;

% Trains the model 
init_time = cputime;
model = multigpOptimise(model, display, iters);
elapsed_time = cputime - init_time;

% Save the results.
capName = dataSetName;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');

% Load complete data to plot the ground truth

[XGT, void, void, fGT] = mapLoadData('ggToy');

ggToyResults(dataSetName, experimentNo, XTemp, yTemp, XGT, fGT);

