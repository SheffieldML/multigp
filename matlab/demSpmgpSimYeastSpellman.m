% DEMSPMGPSIMYEASTSPELLMAN Demonstrate latent force model on gene network data using an
% sparse approximation for the covariance matrix

% MULTIGP

randn('state', 1e6)
rand('twister', 1e6)

addToolboxes(0,1)
dataSetName = 'yeastSpellmanRed';

experimentNo = 32;

% load data
[xTemp, yTemp, tfNames, connect] = mapLoadData(dataSetName);

% index = find(strcmp('ACE2', tfNames)==1);
% whichOutput = find(connect(:, index) == 1);
% yTemp = yTemp(whichOutput, :);
% connect = ones(length(whichOutput),2);
% 
%yTemp = yTemp(1:20,:);
%connect = connect(1:20,:);

% Set the Options 
options = multigpOptions('pitc');
options.type = 'multigp';
options.kernType = 'sim';
options.optimiser = 'scg';
options.connect = connect;
options.nlf = size(connect,2);
options.beta = 1;
options.gamma = exp(-2);
options.initialInducingPositionMethod = 'espacedInRange';
options.numActive = 15;
options.isNormalised = true;
options.meanFunction = true;
options.meanFunctionOptions.type = 'sim';
options.meanFunctionOptions.nlf  = options.nlf;

yTemp2 = cell(size(yTemp,1),1);
xTemp2 = cell(size(yTemp,1),1);
for k = 1:size(yTemp,1)
%    tem = floor(10 + 10*rand);     
%     yTemp2{k} = yTemp(k, 1:tem)';    
%     xTemp2{k} = xTemp(1:tem);
    yTemp2{k} = yTemp(k, :)';        
    xTemp2{k} = xTemp;
    yTemp2{k} = exp(yTemp2{k});  
    %yTemp2{k} = (yTemp2{k})/std(yTemp2{k});  
    %options.beta(k) = 1/(var(yTemp2{k}));
    %yTemp2{k} = yTemp2{k}/std(yTemp2{k});    
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
iters = 20;

% Trains the model and counts the training time
trainingTime = cputime;
model = modelOptimise(model, [], [], display, iters);
trainingTime = cputime - trainingTime;

% % Save the results.
capName = dataSetName;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');

 
% simResults(dataSetName, experimentNo, 1)
