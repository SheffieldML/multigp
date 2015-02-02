% DEMGGJURA Demonstrate multigp convolution model on JURA data using
% the FULL covariance matrix.

% MULTIGP

rand('twister', 1e6);
randn('state', 1e6);

dataSetName = 'juraData';
experimentNo = 1;
file = 'Cd';


[XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData([dataSetName file]);

scaleVal = zeros(1,size(yTemp, 2));
biasVal = zeros(1,size(yTemp, 2));
for k =1:size(yTemp, 2),
    biasVal(k) = mean(yTemp{k});
    scaleVal(k) = sqrt(var(yTemp{k}));
end

options = multigpOptions('ftc');
options.kernType = 'gg';
options.optimiser = 'scg';
options.nlf = 1;
options.beta = ones(1, size(yTemp, 2));
options.bias =  [zeros(1, options.nlf) biasVal];
options.scale = [zeros(1, options.nlf) scaleVal];

q = 2;
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

XTest = cell(size(yTemp, 2)+options.nlf,1);

for j=1:options.nlf
    XTest{j} = ones(1, q);
end
for i = 1:size(yTemp, 2)
    XTest{i+options.nlf} = XTestTemp{i};
end


% Creates the model
model = multigpCreate(q, d, X, y, options);

display = 1;
iters = 10;
model = multigpOptimise(model, display, iters);

% Prediction
mu = multigpPosteriorMeanVar(model, XTest);
maerror = mean(abs((yTestTemp{1} - mu{model.nlf + 1})));





