clc
clear

rand('twister',1e6);
randn('state',1e6);
%
addToolboxes(0,1)
dataSetName = 'demDrosMel';
[xTemp, yTemp, xTestTemp, yTestTemp] = mapLoadData(dataSetName);

% Use only the genes Hb, Kr, Gt, Kni
number = 3;
xTemp = xTemp(number);
yTemp = yTemp(number);
xTestTemp = xTestTemp(number);
yTestTemp = yTestTemp(number);


X = cell(1,size(yTemp, 2));
y = cell(1,size(yTemp, 2));
for i = 1:size(yTemp, 2)
    y{i} = yTemp{i};
    X{i} = xTemp{i};
end

q = size(xTemp{1},2);
nout = size(yTemp, 2);
isGG = 0;
options = gpOptions('ftc');
options.optimiser = 'scg';
if isGG
    optionsk.isArd = true;
    kernType = 'gaussian';
    kernel = kernType;
    kernel(1) = upper(kernel(1));
    nameOfThisFile = ['drosMelNoTransfer' kernel 'GeneHb'];

else
    optionsk.nTerms = 30;
    optionsk.includeIC = false;
    optionsk.includeIndSens = false;
    optionsk.pde = 'sin';
    kernType = 'heat';
    kernel = kernType;
    kernel(1) = upper(kernel(1));    
    nameOfThisFile = ['drosMelNoTransfer' kernel 'nTerms' num2str(optionsk.nTerms) 'GeneHb']; 
end

options.kern = {'cmpnd', {'parametric', optionsk, kernType}, 'white'};
% options.scale2var1 = 1;

iters = 100;
display = 1;
totNfolds = 10;
numberSeed = 1e3;
typeOfInit = 0;



rand('twister',numberSeed);
SEEDS = rand('twister');
SEEDS = SEEDS(1:totNfolds);
SEEDS(1) = 1e6;

modelT = cell(totNfolds,length(yTemp));
elapsed_time = zeros(totNfolds,length(yTemp));
mae = zeros(totNfolds, length(yTemp));
mse = zeros(totNfolds, length(yTemp));
smse = zeros(totNfolds, length(yTemp));
msll = zeros(totNfolds, length(yTemp));
coefR = zeros(totNfolds, length(yTemp));

mu = cell(1);
varsigma = cell(1);
base = 1;

for nFolds = 1:totNfolds
    rand('twister', double(SEEDS(nFolds)));
    ut = unique(xTemp{1}(:,1));
    ut2 = length(ut);
    ut = ut(base+1:end);
    us = unique(xTemp{1}(:,2));
    tN = length(ut);
    sN = length(us);    
    trainSpace = 1:40;
    % Form the new yTemp and xTemp
    % Find the global indexes for training and for testing
    indxGTrain = zeros(length(ut)*length(trainSpace),1);
    startVal = 1;
    endVal = 0;
    for i=1:length(ut)
        endVal = endVal + length(trainSpace);
        tempoSpace = sN*(i-1)+1:sN*i;
        indxGTrain(startVal:endVal, 1) = tempoSpace(trainSpace);
        startVal = endVal + 1;
    end
    indxGTest = 1:tN*sN;
    indxGTest(indxGTrain) = [];
    indxGTrain = indxGTrain + base*sN;
    indxGTest = indxGTest + base*sN;
    X = cell(length(yTemp),1);
    y = cell(length(yTemp),1);
    XTest = cell(length(yTemp),1);
    yTest = cell(length(yTemp),1);        
    for i = 1:length(yTemp)
        y{i} = yTemp{i}(indxGTrain);
        X{i} = xTemp{i}(indxGTrain, :);
        yTest{i} = yTestTemp{i}(indxGTest);
        XTest{i} = xTestTemp{i}(indxGTest, :);
    end
    for j =1:length(y)
        model = gpCreate(q, 1, X{j}, y{j}, options);
        if ~isGG
            params = modelExtractParam(model);
            index = paramNameRegularExpressionLookup(model, '.* decay');
            params(index(1)) = log(0.1);
            index = paramNameRegularExpressionLookup(model, '.* inverse width time');
            %params(index(1)) = log(1e3 + 1e1*rand(1, length(index(1))));
            params(index(1)) = log(1e1 + 1e1*rand(1, length(index(1))));
            index = paramNameRegularExpressionLookup(model, '.* inverse width space\.');
            params(index(1)) = log(1e2 + 1e2*rand(1, length(index(1))));
            %params(index(1)) = log(1e1 + 1e1*rand(1, length(index(1))));
            model = modelExpandParam(model, params);
%             model.kern.comp{2}.variance = std(y{j});            
        end
        initt = cputime;
        model = modelOptimise(model, display, iters);
        elapsed_time(nFolds, j) = cputime - initt;
        modelT{nFolds, j} = model;
        [mu{1}, varsigma{1}] = gpPosteriorMeanVar(model, XTest{j});        
        [mae(nFolds, j), mse(nFolds, j), smse(nFolds, j), msll(nFolds, j), ...
            coefR(nFolds, j)] = multigpErrorMeasures(y(j), yTest(j), mu, varsigma, 1);
        % For comparison with other methods, we should partition the data
        % as needed for multi output in the ICM context
        numDataTrain = length(trainSpace);
        numDataTest = sN - numDataTrain;
        y2 = cell(1, ut2 -base);
        yTest2 = cell(1, ut2 -base);
        mu2 = cell(1, ut2 -base);
        varsigma2 = cell(1, ut2 -base);
        startVal1 = 1;
        endVal1 = 0;
        startVal2 = 1;
        endVal2 = 0;
        for i=base+1:ut2
            endVal1 = endVal1 + numDataTrain;
            endVal2 = endVal2 + numDataTest;
            y2{i-base} = y{j}(startVal1:endVal1);
            yTest2{i-base} = yTest{j}(startVal2:endVal2);
            mu2{i-base} = mu{j}(startVal2:endVal2);
            varsigma2{i-base} = varsigma{j}(startVal2:endVal2);
            startVal1 = endVal1 + 1;
            startVal2 = endVal2 + 1;
        end
        [mae2, mse2, smse2, msll2, coefR2] = multigpErrorMeasures(y2, yTest2, mu2, varsigma2, tN);
    end
end

save([nameOfThisFile '.mat'], 'modelT', 'mae', 'mse', 'smse', ...
    'msll', 'coefR', 'elapsed_time');
