% TOY1DGGFITCMISSING Sparse multigp on TOY data using FITC with missing Data

% MULTIGP

rand('twister', 1e5);
randn('state', 1e5);

dataSetName = 'ggToyMissing';
experimentNo = 2;

[XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData(dataSetName);

options = multigpOptions('fitc');
options.kernType = 'gg';
options.optimiser = 'scg';
options.nlf = 1;
options.initialInducingPositionMethod = 'espaced';
options.numActive = 30;
options.beta = ones(1, size(yTemp, 2));
options.fixInducing = false;

X = cell(size(yTemp, 2),1);
y = cell(size(yTemp, 2),1);

for i = 1:size(yTemp, 2)
  y{i} = yTemp{i};
  X{i} = XTemp{i};
end

q = 1;
d = size(yTemp, 2);

% Creates the model
model = multigpCreate(q, d, X, y, options);

params = modelExtractParam(model);
index = paramNameRegularExpressionLookup(model, 'multi .* inverse .*');
params(index) = log(100); 
model = modelExpandParam(model, params);

display = 1;
iters = 2000;

% Train the model 
init_time = cputime; 
model = multigpOptimise(model, display, iters);
elapsed_time = cputime - init_time;

% Save the results.
capName = dataSetName;
capName(1) = upper(capName(1));
save(['demSpmgp' capName num2str(experimentNo) '.mat'], 'model');

[XGT, void, void, fGT] = mapLoadData('ggToy');

ggSpmgpToyResults(dataSetName, experimentNo, XTemp, yTemp, ...
    XGT, fGT);

