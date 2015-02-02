% DEMCMU49BALANCEARM5 Demonstrate latent force model on CMU data. In this
% experiment the data from the left arm during two motions (namely 18 and
% 19) from subject 49 in the data base is used to train the LFM.
%
% COPYRIGHT : Mauricio Alvarez and Neil D. Lawrence, 2008
%
% MODIFICATIONS : David Luengo, 2009

% MULTIGP


rand('seed', 1e6);
randn('seed', 1e6);

dataSetName = 'cmu49BalanceArm';
experimentNo = 2;

% load data
[yTemp, void, yTestTemp, void] = lvmLoadData(dataSetName);

% Get the time index.
fps = 120/32;

% Scale the ouputs

scaleVal = sqrt(sum(var(yTemp)));
yTemp = yTemp/scaleVal;
yTestTemp = yTestTemp/scaleVal;

% Set the Options 
options.type = 'multigp';
options.numModels = 2;
options.compOptions = multigpOptions('ftc');
% options.compOptions.optimiser = 'conjgrad';
options.compOptions.kernType = 'lfm';
options.compOptions.nlf = 2;
options.compOptions.meanFunction = true;
options.compOptions.meanFunctionOptions.type = 'lfm';
options.compOptions.meanFunctionOptions.nlf  = options.compOptions.nlf;
options.separate = [];

% Set the inputs and outputs in the correct format
X = cell(1, options.numModels);
y = cell(1, options.numModels);
XTest = cell(1, options.numModels);
yTest = cell(1, options.numModels);
for k = 1:options.numModels
    for i = 1:options.compOptions.nlf
        X{k}{i} = 0;
        y{k}{i} = 0;
    end
    for i = 1:size(yTemp, 2)
        if k == 1,
            y{k}{i+options.compOptions.nlf} = yTemp(1:35, i);
            X{k}{i+options.compOptions.nlf} = (1:35)'/fps;
        else
            y{k}{i+options.compOptions.nlf} = yTemp(36:end, i);
            X{k}{i+options.compOptions.nlf} =(1:30)'/fps;
        end
    end
end
for i = 1:options.compOptions.nlf
    XTest{i} = 0;
    yTest{i} = 0;
end
for i = 1:size(yTemp, 2)
    if i < 4 ,
        yTest{i+options.compOptions.nlf} = yTestTemp(1:29, i);
        XTest{i+options.compOptions.nlf} = (1:29)'/fps;
    else
        yTest{i+options.compOptions.nlf} = yTestTemp(1, i);
        XTest{i+options.compOptions.nlf} =(1)'/fps;
    end
end

% Set the input and ouput dimensions
q = 1;
d = size(yTemp, 2) + options.compOptions.nlf;

% Creates the model
model = multimodelCreate(q, d, {X{1}; X{2}}, {y{1}; y{2}}, options);
model.optimiser = model.comp{1}.optimiser;

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
iters = 1000;

% Trains the model and counts the training time
trainingTime = cputime;
model = modelOptimise(model, [], [], display, iters);
trainingTime = cputime - trainingTime;
disp(['Training time : ' num2str(trainingTime) ' s.']);

% Save the results.
capName = dataSetName;
capName(1) = upper(capName(1));
save(['demLfm' capName num2str(experimentNo) '.mat'], 'model');

err = cmu49BalanceResults2(dataSetName, 'lfm', experimentNo);