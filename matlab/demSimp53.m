% DEMSIMP53 Demonstrate latent force model on gene expression data.

% MULTIGP

rand('seed', 1e6);
randn('seed', 1e6);

dataSetName = 'demp53_5genes';
experimentNo = 1;


% load data
[xTemp, yTemp, xTestTemp, yTestTemp] = mapLoadData(dataSetName);


% Set the Options 
options.type = 'multigp';
options.numModels = 3;
options.compOptions = multigpOptions('ftc');
options.compOptions.optimiser = 'scg';
options.compOptions.kernType = 'sim';
options.compOptions.nlf = 2;
%options.compOptions.typeLf = [2 1];
options.compOptions.meanFunction = true;
options.compOptions.meanFunctionOptions.type = 'sim';
options.compOptions.meanFunctionOptions.nlf  = options.compOptions.nlf;
options.compOptions.includeInd = 1;
options.separate = [];
options.optimiser = 'scg';

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
        switch k
        case  1
            y{k}{i+options.compOptions.nlf} = yTemp(1:7, i);
            X{k}{i+options.compOptions.nlf} = xTemp;            
        case 2
            y{k}{i+options.compOptions.nlf} = yTemp(8:14, i);
            X{k}{i+options.compOptions.nlf} = xTemp;
        case 3
            y{k}{i+options.compOptions.nlf} = yTemp(15:21, i);
            X{k}{i+options.compOptions.nlf} = xTemp;
        end
    end
end

% Set the input and ouput dimensions
q = 1;
d = size(yTemp, 2) + options.compOptions.nlf;

% Creates the model
model = multimodelCreate(q, d, {X{1}; X{2}; X{3}}, {y{1}; y{2}; y{3}}, options);

% Creates the model with separate parameters for all components
numComp = options.compOptions.nlf + options.compOptions.includeInd + ...
    options.compOptions.includeNoise;
paramIndVarianceNoise =  paramNameRegularExpressionLookup(model, ['multi ' ...
    num2str(numComp) '.* variance']);
paramIndBasal = paramNameRegularExpressionLookup(model, '. sim .* basal');


options.separate = [paramIndVarianceNoise  paramIndBasal];
model = multimodelCreate(q, d, {X{1}; X{2};X{3}}, {y{1}; y{2};y{3}}, options);

% Make a better initilization of the basal rate parameter
params = modelExtractParam(model);

if ~isempty(options.separate)
    for k = 1:model.numModels
        for i = 1:size(yTemp, 2)
            paramInd = paramNameRegularExpressionLookup(model, ...
                ['multimodel ' num2str(k) '.* sim ' num2str(i)  ' basal']);
            params(paramInd) = log(mean(yTemp((k-1)*7+1:k*7,i)));
        end
    end
else
    for i = 1:size(yTemp, 2)
        paramInd = paramNameRegularExpressionLookup(model, ['.* sim ' num2str(i)  ' basal']);
        % We use the first set of outputs to initialize the basal mean and then these values
        % are shared among the three components. This is done like this for simplicity,
        % but in reality each individual output should have its own basal rate.
        params(paramInd) = log(mean(yTemp(1:7,i)));
    end
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
 

display = 2;
iters = 200;

% Trains the model and counts the training time
trainingTime = cputime;
model = modelOptimise(model, [], [], display, iters);
trainingTime = cputime - trainingTime;

% Save the results.
capName = dataSetName;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');

simResults(dataSetName, experimentNo, 1)


%err = cmu49BalanceResults2(dataSetName, experimentNo);