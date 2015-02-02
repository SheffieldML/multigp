% ROBOTGGIND ..



clc
clear

rand('twister', 1e5);
randn('state', 1e5);
addToolboxes(0,1)
dataSetName = 'robotWireless';
experimentNo = 3;

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
iters = 1000;
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

% Save the results.
capName = dataSetName;
capName(1) = upper(capName(1));
save(['demSpmgp' capName num2str(experimentNo) '.mat'], 'model');

