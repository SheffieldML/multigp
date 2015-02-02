% DEMGGWHITETOY1 Demo of full multi output GP with missing data and GGWHITE
% kernel
% FORMAT
% DESC Demo of Full Multi Output Gaussian Process with missing data and
% GGWHITE kernel

% MULTIGP

clc
clear
rand('twister',1e6);
randn('state',1e6);
%
dataSetName = 'ggwhiteToyMissing';
experimentNo = 1;

[XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData(dataSetName);

options = multigpOptions('ftc');
options.kernType = 'ggwhite';
options.optimiser = 'scg';
options.nlf = 1;

X = cell(size(yTemp, 2)+options.nlf,1);
y = cell(size(yTemp, 2)+options.nlf,1);

for j=1:options.nlf
   y{j} = [];
   X{j} = 1;
end
for i = 1:size(yTemp, 2)
  y{i+options.nlf} = yTemp{i};
  X{i+options.nlf} = XTemp{i};
end

q = 1;
d = size(yTemp, 2) + options.nlf;

% Creates the model
model = multigpCreate(q, d, X, y, options);

params = modelExtractParam(model);
index = paramNameRegularExpressionLookup(model, 'multi .* inverse .*');
params(index) = log(100); 
model = modelExpandParam(model, params);

display = 2;
iters = 1000;

% Trains the model 
init_time = cputime;
model = multigpOptimise(model, display, iters);
elapsed_time = cputime - init_time;

% Save the results.
capName = dataSetName;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');

ggToyResults(dataSetName, experimentNo, XTemp, yTemp);

