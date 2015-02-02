% DEMGGTOY1 Demo of full multi output GP with missing data.
% FORMAT
% DESC Demo of Full Multi Output Gaussian Process with missing data. In this demo, we use the
% Gaussian Kernel for all the covariances (or Kernels) involved and only one hidden function.

% MULTIGP

clc
clear
rand('twister',1e5);
randn('state',1e5);
%
dataSetName = 'ggToyMissing';
experimentNo = 1;

[XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData(dataSetName);

options = multigpOptions('ftc');
options.kernType = 'gg';
options.optimiser = 'scg';
options.nlf = 1;

q = 1; % Input dimension
d = size(yTemp, 2) + options.nlf;

X = cell(size(yTemp, 2)+options.nlf,1);
y = cell(size(yTemp, 2)+options.nlf,1);

% When we want to include the structure of the latent force kernel within
% the whole kernel structure, and we don't have access to any data from the
% latent force, we just put ones in the vector X and empty in the vector y.

for j=1:options.nlf
   y{j} = [];
   X{j} = zeros(1, q);  
end
for i = 1:size(yTemp, 2)
  y{i+options.nlf} = yTemp{i};
  X{i+options.nlf} = XTemp{i};
end

% Creates the model
model = multigpCreate(q, d, X, y, options);

params = modelExtractParam(model);
index = paramNameRegularExpressionLookup(model, 'multi .* inverse .*');
params(index) = log(100+10*rand(1,length(index)));
model = modelExpandParam(model, params);

display = 1;
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

