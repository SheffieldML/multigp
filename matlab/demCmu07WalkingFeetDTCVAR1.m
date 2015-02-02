% DEMCMU07WALKINGLEGS1 Demonstrate latent force model on CMU data. In this
% experiment the data from both legs during two motions (namely 01 and
% 02) from subject 07 in the data base is used to train the LFM.
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% MULTIGP


rand('seed', 1e6);
randn('seed', 1e6);

dataSetName = 'cmu07WalkFeet';
experimentNo = 1;

% load data
[yTemp, void, yTestTemp, void] = lvmLoadData(dataSetName);

% yTemp = yTemp(:,[1 3 4 10]);


% Get the time index.
fps = 120/2;

% Scale the ouputs
isSpeedUp = 2;
% scaleVal = sqrt(sum(var(yTemp)));
% yTemp = yTemp/scaleVal;
% yTestTemp = yTestTemp/scaleVal;

% Set the Options 
options.type = 'multigp';
options.numModels = 2;
options.separate = [];
options.optimiser = 'scg';

% Set the inputs and outputs in the correct format
X = cell(1, options.numModels);
y = cell(1, options.numModels);
XTest = cell(1, options.numModels);
yTest = cell(1, options.numModels);
for k = 1:options.numModels
    options.compOptions{k} = multigpOptions('dtcvar');
    options.compOptions{k}.kernType = 'lfm';
    options.compOptions{k}.nlf = 2;
    options.compOptions{k}.numActive = 50; 
    options.compOptions{k}.initialInducingPositionMethod = 'espaced';
    options.compOptions{k}.beta = 1e-2;
    if isSpeedUp~=0        
        options.compOptions{k}.isSpeedUp = isSpeedUp;
        if isSpeedUp == 2
            options.compOptions{k}.gamma = 1e-1;            
            options.compOptions{k}.includeNoise = false;            
            options.compOptions{k}.kern.isMassFixed.state = true;
            options.compOptions{k}.kern.isMassFixed.massFixedValue = 0.1;
        end
    end
    options.compOptions{k}.useKernDiagGradient = true;
    for i = 1:size(yTemp, 2)
        if k == 1,
            y{k}{i} = yTemp(1:158, i);
            X{k}{i} = (1:158)'/fps;
        else
            y{k}{i} = yTemp(159:end, i);
            X{k}{i} =(1:165)'/fps;
        end
        options.compOptions{k}.bias(i) = mean(y{k}{i});
        options.compOptions{k}.scale(i) = std(y{k}{i});
    end
end

warning('off', 'multigpCreate:NoTyingFunctionGlobalKernel')
warning('off','multiKernParamInit:noCrossKernel')

for i = 1:size(yTemp, 2)
    %     if i < 4 ,
    yTest{i} = yTestTemp(1:208, i);
    XTest{i} = (1:208)'/fps;
    %     else
    %         yTest{i} = yTestTemp(1, i);
    %         XTest{i} =(1)'/fps;
    %     end
end

% Set the input and ouput dimensions
q = 1;
d = size(yTemp, 2);

% Creates the model
createTime = cputime;
model = multimodelCreate(q, d, {X{1}; X{2}}, {y{1}; y{2}}, options);
createTime = cputime - createTime;
model.optimiser = model.comp{1}.optimiser;

% Set parameters associated with the inverse widths to different values for
% symmetry breaking.

if options.compOptions{1}.nlf>1,    
   params = modelExtractParam(model);
   if isSpeedUp==2
       for i = 1:options.compOptions{k}.nlf
           paramInd = paramNameRegularExpressionLookup(model, ['inverse width ' num2str(i)]);
           params(paramInd) = params(paramInd) + 0.01*randn;
       end
   else
       for i = 1:options.compOptions{k}.nlf
           paramInd = paramNameRegularExpressionLookup(model, ['multi ' num2str(i) ...
               ' .* inverse width']);
           params(paramInd) = params(paramInd) + 0.01*randn;
       end
   end
   model = modelExpandParam(model, params);
end
 
display = 1;
iters = 500;

% Trains the model and counts the training time
trainingTime = cputime;
model = modelOptimise(model, [], [], display, iters);
trainingTime = cputime - trainingTime;
disp(['Training time : ' num2str(trainingTime) ' s.']);

[mu, varsigma] = multigpPosteriorMeanVar(model.comp{2},  XTest{1});
[mae, mse, smse, msll] = multigpErrorMeasures(y{2}, yTest, mu(model.comp{2}.nlf+1:end), ...
    varsigma(model.comp{2}.nlf+1:end), model.comp{2}.nout);

% Save the results.
forces = options.compOptions{1}.nlf;
approx = options.compOptions{1}.approx;
kernel = options.compOptions{1}.kernType;
kernel(1) = upper(kernel(1));

capName = dataSetName;
capName(1) = upper(capName(1));
save(['dem' kernel capName num2str(experimentNo) 'Forces' forces 'Approx' approx  '.mat'], 'model');

