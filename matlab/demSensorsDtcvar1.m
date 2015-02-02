% DEMSENSORSDTCVAR1 Sparse VIKs multigp on temperature sensor data 

% MULTIGP

rand('twister', 1e5);
randn('state', 1e5);

dataSetName = 'sensorsTemperature';
experimentNo = 1;

[XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData(dataSetName);

options = multigpOptions('dtcvar');
options.kernType = 'simwhite';
options.optimiser = 'scg';
options.includeInd = true;
options.nlf = 1;
options.initialInducingPositionMethod = 'espaced';
options.numActive = 200;
options.beta = 1e-3*ones(1, size(yTemp, 2));
options.fixInducing = true;
options.isStationary = true;

X = cell(size(yTemp, 2),1);
y = cell(size(yTemp, 2),1);

for i = 1:size(yTemp, 2)
  y{i} = yTemp{i};
  X{i} = XTemp{i};
 options.bias(i) = mean(y{i});
 options.scale(i) = std(y{i});
end

q = 1;
d = size(yTemp, 2);

% Creates the model
model = multigpCreate(q, d, X, y, options);

display = 1;
iters = 50;

% Train the model 
init_time = cputime; 
model = multigpOptimise(model, display, iters);
elapsed_time = cputime - init_time;

% Save the results.
capName = dataSetName;
capName(1) = upper(capName(1));
save(['demSensors' capName num2str(experimentNo) '.mat'], 'model');







