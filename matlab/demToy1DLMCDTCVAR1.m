% DEMTOY1DLMCPITC1 Demo of DTCVAR multi output GP with LMC kernel.
% FORMAT
% DESC Demo of DTVAR Multi Output Gaussian Process with LMC kernel 

% MULTIGP

rand('twister', 1e5);
randn('state', 1e5);

dataSetName = 'ggToyMissing';
experimentNo = 33;

[XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData(dataSetName);

options = multigpOptions('dtcvar');
options.kernType = 'lmc';
options.optimiser = 'scg';
options.nlf = 1;
options.initialInducingPositionMethod = 'espaced';
options.numActive = 30;
options.beta = 1e-1;
options.fixInducing = false;
options.rankCorregMatrix = 2;
options.kern.nout = size(yTemp, 2);
options.kern.rankCorregMatrix = options.rankCorregMatrix;
options.kern.isArd = 0;
options.includeNoise = false;
options.gamma = exp(-2);

warning('off', 'multiKernParamInit:noCrossKernel')

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


display = 2;
iters = 500;

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







