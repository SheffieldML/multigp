% TOY1DGGFTCEXAMPLE Demo of full multi output GP with Gaussian kernel.
% FORMAT
% DESC Demo of Full Multi Output Gaussian Process. 

% MULTIGP

clc
clear
rand('twister',1e5);
randn('state',1e5);
%
dataSetName = 'ggToyICM';
experimentNo = 1;

[XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData(dataSetName);

options = multigpOptions('ftc');
options.kernType = 'gg';
options.optimiser = 'scg';
options.nlf = 5;

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
  options.bias(i+options.nlf)  = mean(y{i+options.nlf});
  options.scale(i+options.nlf) = std(y{i+options.nlf});
end

% Creates the model
model = multigpCreate(q, d, X, y, options);

params = modelExtractParam(model);
index = paramNameRegularExpressionLookup(model, 'multi .* inverse width latent .*');
params(index) = log(1+ rand(1,length(index)));
index = paramNameRegularExpressionLookup(model, 'multi .* sensitivity');
params(index) = 1+ rand(1,length(index));
model = modelExpandParam(model, params);
 
display = 1;
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

[XGT, void, void, fGT] = mapLoadData(dataSetName);

ggToyResults(dataSetName, experimentNo, XTemp, yTemp, XGT, fGT);

