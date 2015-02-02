% DEMCMU49BALANCEARM4 Demonstrate latent force model on CMU data.

% MULTIGP

rand('seed', 1e6);
randn('seed', 1e6);

dataSetName = 'cmu49BalanceArm';
experimentNo = 4;


% load data
[y, void, yTest, void] = lvmLoadData(dataSetName);

% Get the time index.
fps = 120/32;
scaleVal = sqrt(sum(var(y)));
y = y/scaleVal;
X{1}{1} = [0];
Y{1}{1} = [0];
X{2}{1} = [0];
Y{2}{1} = [0];
for i = 1:size(y, 2)
  Y{1}{i+1} = y(1:35, i);
  Y{2}{i+1} = y(36:end, i);
  X{1}{i+1} = (1:35)'/fps;
  X{2}{i+1} = (1:30)'/fps;
end


% Get the time index.
fps = 120/32;
yTest = yTest/scaleVal;
Xtest{1} = [0];
Ytest{1} = [0];
for i = 1:3
  Ytest{i+1} = yTest(1:29, i);
  Xtest{i+1} = (1:29)'/fps;
end
for i = 4:9
  Ytest{i+1} = yTest(1, i);
  Xtest{i+1} = (1)'/fps;
end
options.type = 'multigp';
options.numModels = 2;
options.compOptions = multigpOptions('ftc');
options.compOptions.optimiser = 'conjgrad';
options.compOptions.kernType = 'lfm';
options.compOptions.nlf = 1;
options.compOptions.tieOptions.selectMethod = 'free';
options.separate = [];
q = 1;
d = size(y, 2)+1;

% Creates the model
model = multimodelCreate(q, d, {X{1}; X{2}}, {Y{1}; Y{2}}, options);

% Parameters 4 and 38 are inverse widths, set to different values for
% symmetry breaking.

param = modelExtractParam(model);
param(1) = param(1)+0.01*randn(1);
param(39) = param(39)+0.01*randn(1);

model = modelExpandParam(model, param);

display = 1;
iters = 10;

% Trains the model and counts the training time
model = modelOptimise(model, [], [], display, iters);

% Save the results.
capName = dataSetName;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');


%cmu49BalanceResults(dataSetName, experimentNo);