% DEMSARCOSGGDTC Demonstrate sparse convolution models on SARCOS data.

% MULTIGP

rand('twister', 1e6)
randn('state', 1e6)

dataSetName = 'sarcosMultiGP';
experimentNo = 2;


[XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData(dataSetName);


scaleVal = zeros(1,size(yTemp, 2));
biasVal = zeros(1,size(yTemp, 2));
for k =1:size(yTemp, 2),
    biasVal(k) = mean(yTemp{k});
    scaleVal(k) = sqrt(var(yTemp{k}));
end

options = multigpOptions('dtc');
options.kernType = 'gg';
options.optimiser = 'scg';
options.nlf = 1;
options.initialInducingPositionMethod = 'kmeansHeterotopic';
options.numActive = 50;
options.fixInducing = 0;
options.bias = biasVal;
options.scale = scaleVal;
options.beta = ones(1, size(yTemp, 2));


X = cell(size(yTemp, 2),1);
y = cell(size(yTemp, 2),1);

for i = 1:size(yTemp, 2)
    y{i} = yTemp{i};
    X{i} = XTemp{i};
end

q = 2;
d = size(yTemp, 2);

% Creates the model
model = multigpCreate(q, d, X', y, options);


params = modelExtractParam(model);
index = paramNameRegularExpressionLookup(model, '.* inverse .*');
params(index) = log(100);
model = modelExpandParam(model, params);

display = 1;
iters = 1000;
model = multigpOptimise(model, display, iters);

% Prediction
mu = multigpPosteriorMeanVar(model, XTestTemp);
maerror = mean(abs((yTestTemp{1} - mu{model.nlf + 1})));





