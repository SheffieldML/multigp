% DEMSPMGPSIMYEASTSPELLMANPITC Latent force model on Yeast Spellman Data
% DESC
% Latent force model using a SIM kernel and applied to the cell cycle of
% Yeast data due to Spellman et al (1998) with a PITC approximation for 
% the covariance

% MULTIGP

randn('state', 1e6)
rand('twister', 1e6)

dataSetName = 'yeastSpellman';

experimentNo = 52;

% load data
[xTemp, yTemp, tfNames, connect] = mapLoadData(dataSetName);

% Set the Options 
options = multigpOptions('pitc');
options.type = 'multigp';
options.kernType = 'sim';
options.optimiser = 'scg';
options.connect = connect;
options.nlf = size(connect,2);
options.beta = 1e-3;
options.gamma = exp(-2);
options.initialInducingPositionMethod = 'espacedInRange';
options.numActive = 15;
options.isNormalised = true;
options.isNegativeS = true;
options.meanFunction = true;
options.meanFunctionOptions.type = 'sim';
options.meanFunctionOptions.nlf  = options.nlf;

yTemp2 = cell(size(yTemp,1),1);
xTemp2 = cell(size(yTemp,1),1);
for k = 1:size(yTemp,1)
    yTemp2{k} = yTemp(k, :)';        
    xTemp2{k} = xTemp;
    yTemp2{k} = exp(yTemp2{k});  
end

        
% Set the input dimension
q = 1;
d = size(yTemp2,1);

if options.nlf > 1
   warning('off','multiKernParamInit:noCrossKernel');
end

% Creates the model with equal parameters for all components
creatingTime = cputime;
model = spmultimodelCreate(q, d, xTemp2, yTemp2, options);
creatingTime = cputime - creatingTime;

display = 1;
iters = 1000;

% Trains the model and counts the training time
trainingTime = cputime;
model = modelOptimise(model, [], [], display, iters);
trainingTime = cputime - trainingTime;

% % Save the results.
capName = dataSetName;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');

 
