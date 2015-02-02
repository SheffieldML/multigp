% DEMSPMGPNOISETOY3 Demonstrate latent force model on TOY data with noise
% input function

% In this demo, we use the
% Gaussian White Kernel for all the covariances (or Kernels) involved and only one hidden function.


% MULTIGP

rand('twister', 1e6);
randn('state', 1e6);

dataSetName = 'ggToyMissing';
experimentNo = 4;

[XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData(dataSetName);

options = multigpOptions('dtc');
options.kernType = 'ggwhite';
options.optimiser = 'scg';
options.nlf = 1;
options.initialInducingPositionMethod = 'randomComplete';
options.numActive = 30;
options.fixInducing = 0;
options.beta = 1e-3;

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

% Change the variance (the amplitude of the kernel)
if options.fixInducing == 0
    params = modelExtractParam(model);
    %for i = 1:model.nout,
    paramInd = paramNameRegularExpressionLookup(model, 'X_u .*');
    initialLoc = 0.05*randn(1,length(model.X_u));
    params(paramInd) = initialLoc;
    %end
    model = modelExpandParam(model, params);
end

display = 1;
iters = 1000;

% Trains the model 
init_time = cputime; 
model = multigpOptimise(model, display, iters);
elapsed_time = cputime - init_time;

% Save the results.
capName = dataSetName;
capName(1) = upper(capName(1));
save(['demSpmgp' capName num2str(experimentNo) '.mat'], 'model');


ggSpmgpToyResults(dataSetName, experimentNo, XTemp, yTemp, initialLoc);

