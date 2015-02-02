% DEMPITCCMU49BALANCE Trains latent force model for motion 19, CMU dataset

% MULTIGP

rand('seed', 1e6);
randn('seed', 1e6);


dataSetName = 'cmu49Balance';
experimentNo = 13;


% load data
[yTemp, void, yTestTemp, void] = lvmLoadData(dataSetName);

% Get the time index.
fps = 120/32;

% Scale the ouputs

scaleVal = sqrt(sum(var(yTemp)));
yTemp = yTemp/scaleVal;
yTestTemp = yTestTemp/scaleVal;

%scaleVal = sqrt(var(yTemp));
%yTemp = yTemp./repmat(scaleVal, size(yTemp, 1),1);
%yTestTemp = yTestTemp./repmat(scaleVal, size(yTestTemp, 1),1);


% Set the Options 
options.type = 'multigp';
options.numModels = 2;
options.compOptions = multigpOptions('pitc');
options.compOptions.optimiser = 'conjgrad';
options.compOptions.kernType = 'lfm';
options.compOptions.nlf = 3;
options.compOptions.initialInducingPositionMethod = 'espacedInRange';
options.compOptions.numActive = 20;
options.compOptions.beta = 1e3;
options.compOptions.meanFunction = true;
options.compOptions.meanFunctionOptions.type = 'lfm';
options.compOptions.meanFunctionOptions.nlf  = options.compOptions.nlf;
options.separate = [];
options.optimiser = 'scg';

% Set the inputs and outputs in the correct format
X = cell(1, options.numModels);
y = cell(1, options.numModels);
XTest = cell(1, options.numModels);
yTest = cell(1, options.numModels);
for k = 1:options.numModels    
    for i = 1:size(yTemp, 2)
        if k == 1,
            y{k}{i} = yTemp(1:35, i);
            X{k}{i} = (1:35)'/fps;
        else
            y{k}{i} = yTemp(36:end, i);
            X{k}{i} =(1:30)'/fps;
        end
    end
end

% Set the input and ouput dimensions
q = 1;
d = size(yTemp, 2);

% Creates the model
model = multimodelCreate(q, d, {X{1}; X{2}}, {y{1}; y{2}}, options);

% Make a better initilization of the basal rate parameter
params = modelExtractParam(model);
for i = 1:size(yTemp, 2)
    paramInd = paramNameRegularExpressionLookup(model, ['. lfm ' num2str(i)  ' basal']);
    params(paramInd) = mean(yTemp(:,i));
end
model = modelExpandParam(model, params);

% Set parameters associated with the inverse widths to different values for
% symmetry breaking.

if options.compOptions.nlf>1,
   params = modelExtractParam(model);
   for i = 1:options.compOptions.nlf
       paramInd = paramNameRegularExpressionLookup(model, ['multi ' num2str(i) ...
           ' .* inverse width']);
       params(paramInd) = params(paramInd) + 0.01*randn;
   end
   model = modelExpandParam(model, params);
end
 

display = 1;
iters = 5000;

% Trains the model and counts the training time
trainingTime = cputime;
model = modelOptimise(model, [], [], display, iters);
trainingTime = cputime - trainingTime;

% Save the results.
capName = dataSetName;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model','trainingTime');

err = cmu49BalanceResults2(dataSetName, 'lfm', experimentNo);
