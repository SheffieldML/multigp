% DEMCMU49BALANCEARM1 Demonstrate latent force model on CMU data.

% MULTIGP

rand('seed', 1e6);
randn('seed', 1e6);

dataSetName = 'cmu49BalanceArm';
experimentNo = 1;


% load data
[y, void, yTest, void] = lvmLoadData(dataSetName);

% Get the time index.
fps = 120/32;
t = 1:size(y, 1);
t = t'/fps;
y = y/sqrt(sum(var(y)));

for i = 1:size(y, 2)
  Y{i} = y(:, i);
  X{i} = t;
end


% Set up model
options = multigpOptions('ftc');
options.optimiser = 'conjgrad';
options.kernType = 'lfm';

options.tieOptions.selectMethod = 'free';

q = 1;
d = size(y, 2);

% Creates the model
model = multigpCreate(q, d, X, Y, options);
%model.scale = repmat(sqrt(sum(var(y))/model.d), 1, model.d);

model.scale = repmat(1, 1, model.d);
model.m = multigpComputeM(model);


display = 2;
iters = 1000;

% Trains the model and counts the training time
model = multigpOptimise(model, display, iters);

% Save the results.
capName = dataSetName;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');


