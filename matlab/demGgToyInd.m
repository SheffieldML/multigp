clc
clear

rand('twister',1e6);
randn('state',1e6);
%
dataSetName = 'ggToyMissing';
experimentNo = 15;
[XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData(dataSetName);
% Configuration of data
q = size(XTemp{1},2);
d = size(yTemp, 2);
nout = size(yTemp, 2);

options.type = 'gp';
options.numModels = nout;
options.compOptions = gpOptions('ftc');
options.compOptions.optimiser = 'scg';
options.compOptions.kern = {'gg', 'white'};
options.separate = [];
options.optimiser = 'scg';
iters =100;
display = 1;

X = cell(1,size(yTemp, 2));
y = cell(1,size(yTemp, 2));

for i = 1:size(yTemp, 2)
    X{i}  = XTemp{i};
    y{i} = yTemp{i};
end
% Configuration of parameters

model = multimodelCreate(q, 1, X, y, options);
params = modelExtractParam(model);
options.separate = 1:length(params);
model = multimodelCreate(q, 1, X, y, options);

params = modelExtractParam(model);
index = paramNameRegularExpressionLookup(model, 'multimodel .* inverse .*');
params(index) = log(100);

model = modelExpandParam(model, params);

model = modelOptimise(model, [], [], display, iters);


