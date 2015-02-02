% DEMSPMGPGGTOY3KL Measures the KL divergence between the full Multigp with
% GG kernel against the approximation using PITC

% MULTIGP

clc
clear
addToolboxes(0,1)
rand('twister',1e6);
randn('state',1e6);
%
dataSetName = 'ggToyMissing';
experimentNo = 1;

try
    capName = dataSetName;
    capName(1) = upper(capName(1));
    load(['dem' capName num2str(experimentNo) '.mat'])
catch
    [XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData(dataSetName);
    options = multigpOptions('ftc');
    options.kernType = 'gg';
    options.optimiser = 'scg';
    options.nlf = 1;
    q = 1; % Input dimension
    d = size(yTemp, 2) + options.nlf;
    X = cell(size(yTemp, 2)+options.nlf,1);
    y = cell(size(yTemp, 2)+options.nlf,1);
    for j=1:options.nlf
        y{j} = [];
        X{j} = zeros(1, q);
    end
    for i = 1:size(yTemp, 2)
        y{i+options.nlf} = yTemp{i};
        X{i+options.nlf} = XTemp{i};
    end
    % Creates the model
    model = multigpCreate(q, d, X, y, options);
    display = 1;
    iters = 1000;
    params = modelExtractParam(model);
    index = paramNameRegularExpressionLookup(model, 'multi .* inverse .*');
    params(index) = log(100);
    model = modelExpandParam(model, params);
    % Trains the model
    model = multigpOptimise(model, display, iters);    
    capName = dataSetName;
    capName(1) = upper(capName(1));
    save(['dem' capName num2str(experimentNo) '.mat'], 'model');
end

llFull = modelLogLikelihood(model);
modelFull = model;
% Trains the sparse model
nFolds = 10;
display = 0;
iters = 1000;
options = multigpOptions('pitc');
options.kernType = 'gg';
options.optimiser = 'scg';
options.nlf = 1;
options.initialInducingPositionMethod = 'randomComplete';
options.beta = 1e3;
options.fixInducing = false;
numActive = [10 20 30 40 50 100 200];
llApprox = zeros(length(numActive), nFolds);
elapsed_time = zeros(length(numActive), nFolds);
% models = cell(length(numActive), nFolds);
for j=1:length(numActive)
    options.numActive = numActive(j);
    for k = 1: nFolds,
        [XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData(dataSetName);
        X = cell(size(yTemp, 2),1);
        y = cell(size(yTemp, 2),1);
        for i = 1:size(yTemp, 2)
            y{i} = yTemp{i};
            X{i} = XTemp{i};
        end
        q = 1;
        d = size(yTemp, 2);
        % Creates the model
        rand('twister',5^(k));
        fprintf('Creating MODEL %d APPROX %s NUMACTIVE %d \n',...
        k, options.approx, options.numActive);
        model = multigpCreate(q, d, X, y, options);
        modelFullAux = modelFull;
        modelFullAux.paramGroups = speye(size(modelFull.paramGroups,1));
        modelAux = model;
        modelAux.paramGroups = speye(size(model.paramGroups,1));
        [paramsFull, namesFull] = modelExtractParam(modelFullAux);
        indexNoTransform =  paramNameRegularExpressionLookup(modelAux, ' .* sensitivity');
        indexToAugment = paramNameRegularExpressionLookup(modelFullAux, ['multi ' ...
            num2str(length(model.kern.comp)) ' white .* variance']);
        count = length(model.fix);
        count2 = 1;
        for i=1:length(namesFull)
            if any(i == indexToAugment)
                count2 = count2 + 1;
                nameWhite = ['multi ' num2str(length(model.kern.comp)) ' white ' ...
                    num2str(count2) ' variance'];
                paramInd = paramNameRegularExpressionLookup(modelAux, nameWhite);
            else
                paramInd = paramNameRegularExpressionLookup(modelAux, namesFull{i});
            end
            if ~isempty(paramInd)
                count = count + 1;
                if any(paramInd == indexNoTransform)
                    model.fix(count).index = paramInd;
                    model.fix(count).value = paramsFull(i);
                else
                    model.fix(count).index = paramInd;
                    model.fix(count).value = expTransform(exp(paramsFull(i)), 'xtoa');
                end
            end
        end
        % Train the model
        fprintf('Optimizing MODEL %d  APPROX %s NUMACTIVE %d \n',...
            k, options.approx, options.numActive);
        tic;
        model = multigpOptimise(model, display, iters);
        elapsed_time(j,k) = toc;
        llApprox(j,k) = modelLogLikelihood(model);
        %models{j,k} = model;
    end
end

save('demSpmgpGgToy3KL.mat', 'llFull', 'llApprox', 'elapsed_time');







