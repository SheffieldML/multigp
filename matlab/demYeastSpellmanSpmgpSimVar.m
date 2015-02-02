% DEMYEASTSPELLMANSPMGPSIMVAR Variational LFM on Yeast Spellman Data 
% DESC
% Latent force model using a SIM kernel and applied to the cell cycle of
% Yeast data due to Spellman et al (1998) with a DTCVAR approximation for 
% the covariance and variational learning of the sensitivities.

% MULTIGP

randn('state', 1e6)
rand('twister', 1e6)

dataSetName = 'yeastSpellmanRed';

experimentNo = 34;

% load data
[xTemp, yTemp, tfNames, connect] = mapLoadData(dataSetName);

% Set the Options 
options = multigpOptions('dtcvar');
options.type = 'multigp';
options.kernType = 'sim';
options.optimiser = 'scg';
options.connect = connect;
options.nlf = size(connect,2);
options.beta = 1;
options.gamma = exp(-2);
options.initialInducingPositionMethod = 'espacedInRange';
options.numActive = 15;
options.meanFunction = true;
options.meanFunctionOptions.type = 'sim';
options.meanFunctionOptions.nlf  = options.nlf;
options.varS = 1;
options.kernOptions.rbf.isNormalised = true;
options.kernOptions.sim.isNormalised = true;
options.kernOptions.sim.isVarS = true;
options.kernOptions.sim.isNegativeS = true;

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
iters = 10;
niters = 20;

% Trains the model and counts the training time
trainingTime = cputime;
model = spmultimodelOptimiseVar(model, display, iters, niters);
trainingTime = cputime - trainingTime;

% % Save the results.
capName = dataSetName;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');

 
