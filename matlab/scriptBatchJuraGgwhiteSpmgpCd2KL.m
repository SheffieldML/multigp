clear
clc
addToolboxes(0,1)

file = 'Cd';
nameFull = 'juraFullGgWhiteKLCd.mat';
nameThisFile = 'scriptBatchJuraGgwhiteSpmgpCd2KL.mat';

rand('twister', 1e6)
randn('state', 1e6)

if exist(nameFull, 'file')
    load(nameFull);
else    
    iters = 1000;
    options = multigpOptions('ftc');
    options.kernType = 'ggwhite';
    options.optimiser = 'scg';
    options.nlf = 1;
    options.beta = 1e3;
    options.isArd = 1;
    dataSetName = ['juraData' file];
    display = 1;
    [XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData(dataSetName);
    scaleVal = zeros(1,size(yTemp, 2));
    biasVal = zeros(1,size(yTemp, 2));
    for k =1:size(yTemp, 2),
        biasVal(k) = mean(yTemp{k});
        scaleVal(k) = sqrt(var(yTemp{k}));
    end
    options.bias =  [zeros(1, options.nlf) biasVal];
    options.scale = [zeros(1, options.nlf) scaleVal];
    X = cell(size(yTemp, 2)+options.nlf,1);
    y = cell(size(yTemp, 2)+options.nlf,1);
    for j=1:options.nlf
        y{j} = [];
        X{j} = zeros(1, 2);
    end
    for i = 1:size(yTemp, 2)
        y{i+options.nlf} = yTemp{i};
        X{i+options.nlf} = XTemp{i};
    end
    XTest = cell(size(yTemp, 2)+options.nlf,1);
    for j=1:options.nlf
        XTest{j} = ones(1, 2);
    end
    for i = 1:size(yTemp, 2)
        XTest{i+options.nlf} = XTestTemp{i};
    end
    q = size(X{end},2);
    d = size(yTemp, 2) + options.nlf;
    % Creates the model
    model = multigpCreate(q, d, X, y, options);
    % Change initial conditions
    params = modelExtractParam(model);
    index = paramNameRegularExpressionLookup(model, 'multi .* inverse .*');
    params(index) = log(100);
    model = modelExpandParam(model, params);
    % Optimization procedure  
    model = multigpOptimise(model, display, iters);
    llFull = modelLogLikelihood(model);
    save(nameFull, 'model')
end

numActive =[50 100 200 500 600];
numFolds = 10;
iters = 1000;
options = multigpOptions('dtc');
options.kernType = 'ggwhite';
options.optimiser = 'scg';
options.nlf = 1;
options.initialInducingPositionMethod = 'randomDataIsotopic';
options.isArd = 1;
options.nVIKs = 1;
options.flagVIKs = true;
options.fixInducing = 0;

[bound, elapsed_time] = demSpmgpJuraBatchKL(file, options, numActive, iters, numFolds, model);

save(nameThisFile, 'bound', 'elapsed_time')
