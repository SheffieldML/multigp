% TOY1DGGINDMISSING Independent GP with missing Data

% MULTIGP

clc
clear
rand('twister',1e5);
randn('state',1e5);
%
dataSetName = 'ggToyMissing';
experimentNo = 15;

[XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData(dataSetName);

%Just the fourth output

XTemp = XTemp(4);
yTemp = yTemp(4);


options = multigpOptions('ftc');
options.kernType = 'gg';
options.optimiser = 'scg';
options.nlf = 1;

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
model = multigpCreate(q, d, X, y, options);

params = modelExtractParam(model);
index = paramNameRegularExpressionLookup(model, 'multi .* inverse .*');
params(index) = log(100);
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

[XGT, void, void, fGT] = mapLoadData('ggToy');

ggToyResults(dataSetName, experimentNo, XTemp, yTemp, XGT(4), fGT(4))

