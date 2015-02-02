clc
clear

rand('twister',1e6);
randn('state',1e6);
%
dataSetName = 'schoolData';
experimentNo = 12;
[XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData(dataSetName);
% Configuration of data
q = size(XTemp{1},2)-1;
d = size(yTemp, 2);
nout = size(yTemp, 2);


options.type = 'gp';
options.numModels = nout;
options.compOptions = gpOptions('ftc');
options.compOptions.optimiser = 'scg';
isRbf = 0;
if isRbf
    options.compOptions.kern = {'rbf', 'white'};    
else
%    options.compOptions.kern = {'ggwhite', 'white'};  % For Heskes approach  
   options.compOptions.kern = {'gg', 'whiteh'}; % For Bonilla approach
end
options.separate = [];
options.optimiser = 'scg';
options.compOptions.scale2var1 = 1;

iters = 100;
display = 1;
totNfolds = 5;
logmod = zeros(totNfolds, nout);
res = zeros(totNfolds, nout);
smse = zeros(totNfolds, nout);
msll = zeros(totNfolds, nout);
meanTrain = zeros(1, nout);
varTrain = zeros(1, nout);

X = cell(1,size(yTemp, 2));
y = cell(1,size(yTemp, 2));
XTest = cell(1, size(yTemp, 2));
yTest = cell(1, size(yTemp, 2));
X2 = cell(1,size(yTemp, 2));
y2 = cell(1,size(yTemp, 2));
XTest2 = cell(1, size(yTemp, 2));
yTest2 = cell(1, size(yTemp, 2));


for nFolds = 1:totNfolds  
%     rand('twister',5^(nFolds));
%     randn('state',5^(nFolds));
    for i = 1:size(yTemp, 2)
        maxSamples = size(yTemp{i}, 1);
        randIndex = randperm(maxSamples);
        samplesTraining = floor(0.75*maxSamples);
        indexTraining = randIndex(1:samplesTraining);
        indexTesting = randIndex(samplesTraining+1:end);
        y{i} = yTemp{i}(indexTraining,:);
        X{i} = XTemp{i}(indexTraining,:);
%         y2{i} = cell2mat(yTestTemp{i}(indexTraining));        
%         X2{i} = cell2mat(XTestTemp{i}(indexTraining));
        yTest{i} =  yTemp{i}(indexTesting,:);
        XTest{i} =  XTemp{i}(indexTesting,:);
        meanTrain(i) = mean(y{i});        
        varTrain(i) = var(y{i},1);
        %varTrain(i) = (y{i} - mean(y{i}))'*(y{i} - mean(y{i}));
%         yTest2{i} = yTestTemp{i}(indexTesting,:);
%         XTest2{i} = XTestTemp{i}(indexTesting,:);                
    end
    
%     % Configuration of parameters    
%     model = multimodelCreate(q, 1, X', y', options);    
%     params = modelExtractParam(model);
%     %options.separate = 1:(length(params)-1);
%     options.separate = 1:(length(params));
%     model = multimodelCreate(q, 1, X', y', options);           
%     for i =1:nout,
%         model.comp{i}.kern = kernSetIndex(model.comp{i}.kern, 1, 1:q);
%         model.comp{i}.kern = kernSetIndex(model.comp{i}.kern, 2, 1:q);
%         if ~isRbf
%             fhandle = str2func([options.compOptions.kern{1} 'KernParamInit']);
%             model.comp{i}.kern.comp{1} = fhandle(model.comp{i}.kern.comp{1}, 1);
% %             model.comp{i}.kern.comp{1} = ggwhiteKernParamInit...
% %                 (model.comp{i}.kern.comp{1}, 1);
%             model.comp{i}.kern.nParams     = model.comp{i}.kern.comp{1}.nParams + ...
%                 model.comp{i}.kern.comp{2}.nParams;
%             model.comp{i}.kern.paramGroups = speye(model.comp{i}.kern.nParams);
%             model.comp{i}.nParams = model.comp{i}.kern.nParams;
%             %model.separateIndices = 1:model.comp{i}.kern.nParams-1;
%             model.separateIndices = 1:model.comp{i}.kern.nParams;
%             model.numSep = length(model.separateIndices);
%             model.numParams = model.comp{1}.nParams;
%             model.sharedIndices = 1:model.numParams;
%             model.sharedIndices(model.separateIndices) = [];
%             model.numParams = model.numParams + (model.numModels-1)*model.numSep;
%         end
%     end


    model = multimodelCreate(q+1, 1, X', y', options);    
    params = modelExtractParam(model);
    options.separate = 1:(length(params)-1);
    model = multimodelCreate(q+1, 1, X', y', options);           
    for i =1:nout,
        model.comp{i}.kern = kernSetIndex(model.comp{i}.kern, 1, 1:q);
        model.comp{i}.kern = kernSetIndex(model.comp{i}.kern, 2, q+1);
        if ~isRbf
            fhandle = str2func([options.compOptions.kern{1} 'KernParamInit']);
            model.comp{i}.kern.comp{1} = fhandle(model.comp{i}.kern.comp{1}, 1);
%             model.comp{i}.kern.comp{1} = ggwhiteKernParamInit...
%                 (model.comp{i}.kern.comp{1}, 1);
            model.comp{i}.kern.nParams     = model.comp{i}.kern.comp{1}.nParams + ...
                model.comp{i}.kern.comp{2}.nParams;
            model.comp{i}.kern.paramGroups = speye(model.comp{i}.kern.nParams);
            model.comp{i}.nParams = model.comp{i}.kern.nParams;
            model.separateIndices = 1:model.comp{i}.kern.nParams-1;
            model.numSep = length(model.separateIndices);
            model.numParams = model.comp{1}.nParams;
            model.sharedIndices = 1:model.numParams;
            model.sharedIndices(model.separateIndices) = [];
            model.numParams = model.numParams + (model.numModels-1)*model.numSep;
        end
    end
    initt = cputime;
    model = modelOptimise(model, [], [], display, iters);   
    elapsed_time = cputime - initt;
    for k =1:nout,
        xTest = XTest{k};
        Ytest = yTest{k};
%         xTest = cell2mat(XTest{k});
%         xTest(:,q+1) = ones(size(xTest,1),1);
%         Ytest = cell2mat(yTest{k});        
        [kk, kkvar] = gpPosteriorMeanVar(model.comp{k}, xTest);       
        coef = corrcoef(kk,Ytest);
        res2(nFolds, k) = (coef(2,1))^2;
%         SStot = std(y{k})^2;
%         SSerr = (Ytest - kk)'*(Ytest - kk)/size(kk,1);
        %SStot = (Ytest - mean(Ytest))'*(Ytest - mean(Ytest));
        SStot = var(Ytest,1);
        SSerr = ((Ytest - kk)'* (Ytest - kk))/size(Ytest,1);
        LLerr = 0.5*log(2*pi*kkvar) + 0.5*((Ytest - kk).^2)./kkvar;
        LLtot = 0.5*log(2*pi*varTrain(k)) + 0.5*((Ytest - meanTrain(k)).^2)/(varTrain(k));
%         kk = gpPosteriorMeanVar(model.comp{k}, X{k});       
%         SStot = (y{k} - mean(y{k}))'*(y{k} - mean(y{k}));
%         SSerr = (y{k} - kk)'* (y{k} - kk);
        msll(nFolds, k) = mean(LLerr - LLtot);
        smse(nFolds, k) = SSerr/SStot;
        res(nFolds, k) = 1 - (SSerr/SStot);   
%         logmod(nFolds, k) = modelLogLikelihood(model.comp{k});
    end
%     clear model
    options.separate = [];
end