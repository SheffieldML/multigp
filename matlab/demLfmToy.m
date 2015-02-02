% DEMLFMTOY Demonstrate latent force model on TOY data.

% MULTIGP

rand('seed', 1e6);
randn('seed', 1e6);

dataSetName = 'lfmOde';
experimentNo = 1;


% load data
[XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData(dataSetName);

% Set equal variances 
scaleVal = sqrt(sum(var(yTemp)));
yTemp = yTemp/scaleVal;
yTestTemp = yTestTemp/scaleVal;

options = multigpOptions('ftc');
options.kernType = 'lfm';
options.optimiser = 'scg';
options.nlf = 1;
options.meanFunction = true;
options.meanFunctionOptions.type = 'lfm';
options.meanFunctionOptions.nlf  = options.nlf;

X = cell(size(yTemp, 2)+options.nlf,1);
y = cell(size(yTemp, 2)+options.nlf,1);
XTest = cell(size(yTemp, 2)+options.nlf,1);
yTest = cell(size(yTemp, 2)+options.nlf,1);

for i = 1:options.nlf
    X{i} = 0;
    y{i} = 0;
end
for i = 1:size(yTemp, 2)
  y{i+options.nlf} = yTemp(:, i);
  X{i+options.nlf} = XTemp;
end

for i = 1:options.nlf,
    XTest{i} = 0;
    yTest{i} = 0;
end
for i = 1:size(yTemp, 2)
    yTest{i+options.nlf} = yTestTemp(:, i);
    XTest{i+options.nlf} = XTestTemp;
end

q = 1;
d = size(yTemp, 2) + options.nlf;

% Creates the model
model = multigpCreate(q, d, X, y, options);

% Make a better initilization of the basal rate parameter

params = modelExtractParam(model);
for i = 1:size(yTemp, 2)
    paramInd = paramNameRegularExpressionLookup(model, ['. lfm ' num2str(i)  ' basal']);
    params(paramInd) = mean(yTemp(:,i));
end
model = modelExpandParam(model, params);

% Change the initialization of inverse widths 

if options.nlf>1,
   params = modelExtractParam(model);
   for i = 1:options.nlf
       paramInd = paramNameRegularExpressionLookup(model, ['multi ' num2str(i) ...
           ' .* inverse width']);
       params(paramInd) = params(paramInd) + 0.01*randn;
   end
   model = modelExpandParam(model, params);
end

display = 2;
iters = 5;

% Trains the model 
init_time = cputime;
model = multigpOptimise(model, display, iters);
elapsed_time = cputime - init_time;

% Save the results.
capName = dataSetName;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');


lfmToyResults(dataSetName, experimentNo);