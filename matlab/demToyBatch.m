function [mean_sserror, desv_sserror, mean_logprob, desv_logprob, elapsed_time] = ...
    demToyBatch(options,iters, numFolds, dataSetName)

% MULTIGP

if nargin <4
    dataSetName = 'ggToy';
end
display = 1;
elapsed_time =  zeros(numFolds,1);


[XTemp, yTemp, XTestTemp, yTest] = mapLoadData(dataSetName);

desv_sserror = zeros(numFolds,size(yTemp, 2));
mean_sserror = zeros(numFolds,size(yTemp, 2));
desv_logprob = zeros(numFolds,size(yTemp, 2));
mean_logprob = zeros(numFolds,size(yTemp, 2));

q = 1; % Input dimension
d = size(yTemp, 2) + options.nlf;

X = cell(size(yTemp, 2)+options.nlf,1);
y = cell(size(yTemp, 2)+options.nlf,1);
modelBatch = cell(1, numFolds);

rand('twister',1e6);
randn('state',1e6);
for fold =1: numFolds    
    for j=1:options.nlf
        y{j} = [];
        X{j} = zeros(1, q);
    end
    maxSamples = size(yTemp{1}, 1);
    randIndex = randperm(maxSamples);
    samplesTraining = floor(0.4*maxSamples);
    indexTraining = randIndex(1:samplesTraining);
    indexTesting = randIndex(samplesTraining+1:end);
    for i = 1:size(yTemp, 2)
        y{i+options.nlf} = yTemp{i}(indexTraining);
        X{i+options.nlf} = XTemp{i}(indexTraining);
        yTest{i} = yTemp{i}(indexTesting);
    end
    XTestTemp = XTemp{1}(indexTesting);
    fprintf('Creating MODEL %d APPROX %s \n', fold, options.approx);
    model = multigpCreate(q, d, X, y, options);
    % Initialize parameters
    % Change initial conditions for inverse widths
    params = modelExtractParam(model);
    index = paramNameRegularExpressionLookup(model, 'multi .* inverse .*');
    params(index) = log(100);
    model = modelExpandParam(model, params);
    fprintf('Optimizing MODEL %d APPROX %s \n', fold, options.approx)
    tic
    model = multigpOptimise(model, display, iters);
    elapsed_time(fold) = toc;
    fprintf('Elapsed time optimization: %f hours\n', elapsed_time(fold)/3600)
    % Prediction
    [mu, varsigma] = multigpPosteriorMeanVar(model, XTestTemp);
    modelBatch{fold} = model;
    % Standarized performance measures
    for k= 1:model.nout,
        mserror = ((yTest{k} - mu{k+model.nlf}).^2)/var(yTest{k});
        sserror = ((yTest{k} - mu{k+model.nlf}).^2)./(2*varsigma{k+model.nlf});
        logprob = 0.5*log(2*pi.*varsigma{k+model.nlf}) + sserror - 0.5*log(2*pi.*var(y{k+model.nlf})) ...
            - ((yTest{k} - mean(y{k+model.nlf})).^2)./(2*var(y{k+model.nlf}));
        mean_logprob(fold,k) = mean(logprob);
        desv_logprob(fold,k) = std(logprob);
        mean_sserror(fold,k) = mean(mserror);
        desv_sserror(fold,k) =  std(mserror);
    end
    save([dataSetName upper(options.kernType(1)) options.kernType(2:end) upper(options.approx) ...
                    'NumberForces_' num2str(options.nlf)])
end



