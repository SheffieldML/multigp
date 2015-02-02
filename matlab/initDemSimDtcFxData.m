function [paramsInit, paramsInit2] = initDemSimDtcFxData(X, y, count, options, flags);

% INITDEMSIMDTCFXDATA initialisation of the parameters of the SIM and SIM-WHITE models using an independent GP for each output.
%
% COPYRIGHT : David Luengo, 2009
%
% SEEALSO : demSimDtcFxData, demSimwhiteDtcFxData

% MULTIGP

% Parameters of the function

niter = 25;
niter2 = 25;

isNormalised = flags(1);
isStationary = flags(2);

indGpOptions = options;
indGpOptions.nlf = 1;
if isfield(indGpOptions, 'typeLf')
    indGpOptions = rmfield(indGpOptions, 'typeLf');
end

% Dummy model to initialise the parameters matrix
XInd = X{1}(:, 1);
yInd = y{1}(:, 1);
indGpModel = multigpCreate(1, 1, XInd, yInd, indGpOptions);
params = zeros(size(y{1}, 2), size(modelExtractParam(indGpModel), 2));
% names = cell(size(y{1}, 2), size(modelExtractParam(indGpModel), 2));

% Training the independent models and collecting the parameters

for i=1:size(y{1}, 2)
    XInd = X{1}(:, i);
    yInd = y{1}(:, i);
    indGpModel = multigpCreate(1, 1, XInd, yInd, indGpOptions);
    indGpModel.X_u = linspace(count(1), count(end), indGpOptions.numActive)';
    indGpModel.kern.comp{1}.comp{1}.isNormalised = isNormalised;
    indGpModel.kern.comp{1}.comp{2}.isNormalised = isNormalised;
    indGpModel.kern.comp{1}.comp{2}.isStationary = isStationary;
    par = modelExtractParam(indGpModel);
    indGpModel = modelExpandParam(indGpModel, par);
    indGpModel = modelOptimise(indGpModel, 1, niter);
    modelDisplay(indGpModel);
    params(i, :) = modelExtractParam(indGpModel);
end

% If there are two types of latent forces then we have to initialise a
% second kernel

if isfield(options, 'typeLf')
    % Dummy model to initialise the parameters matrix
    XInd = X{1}(:, 1);
    yInd = y{1}(:, 1);
    indGpModel = multigpCreate(1, 1, XInd, yInd, indGpOptions);
    indGpOptions.kernType = 'simwhite';
    params2 = zeros(size(y{1}, 2), size(modelExtractParam(indGpModel), 2));
%     names2 = cell(size(y{1}, 2), size(modelExtractParam(indGpModel), 2));
    for i=1:size(y{1}, 2)
        XInd = X{1}(:, i);
        yInd = y{1}(:, i);
        indGpModel = multigpCreate(1, 1, XInd, yInd, indGpOptions);
        indGpModel.X_u = linspace(count(1), count(end), indGpOptions.numActive)';
        indGpModel.kern.comp{1}.comp{1}.isNormalised = isNormalised;
        indGpModel.kern.comp{1}.comp{2}.isNormalised = isNormalised;
        indGpModel.kern.comp{1}.comp{2}.isStationary = isStationary;
        par = modelExtractParam(indGpModel);
        indGpModel = modelExpandParam(indGpModel, par);
        indGpModel = modelOptimise(indGpModel, 1, niter);
        modelDisplay(indGpModel);
        params2(i, :) = modelExtractParam(indGpModel);
    end
end

if isfield(options, 'typeLf')
    numSimLf = options.typeLf(1);
    numSimwhiteLf = options.typeLf(2);
else
    numSimLf = options.nlf;
    numSimwhiteLf = 0;
end

params = exp(params);
if numSimLf==1
    invWidthSim = mean(params(:, 1));
else
    currentFolder = pwd;
    try
        cd([matlabroot '/toolbox/stats']);
        [void, invWidthSim] = kmeans(params(:, 1), numSimLf);
        cd(currentFolder);
    catch
        error('You need to hav the Statistics Toolbox installed in the default path');
    end
end

if numSimwhiteLf > 0
    params2(:,[1:3 5:7]) = exp(params2(:,[1:3 5:7]));
    if numSimwhiteLf==1
        invWidthSimwhite = mean(params2(:, 1));
    else
        currentFolder = pwd;
        cd([matlabroot '/toolbox/stats']);
        [void, invWidthSimwhite] = kmeans(params2(:, 1), numSimwhiteLf);
        cd(currentFolder);
    end
end

paramsInit = zeros(size(y{1}, 2), size(modelExtractParam(indGpModel), 2), numSimLf);
% namesInit = cell(size(y{1}, 2), size(modelExtractParam(indGpModel), 2), numSimLf);
for i = 1:numSimLf
    for j = 1:size(y{1}, 2)
        XInd = X{1}(:, j);
        yInd = y{1}(:, j);
        indGpOptions.kernType = options.kernType;
        indGpModel = multigpCreate(1, 1, XInd, yInd, indGpOptions);
        indGpModel.X_u = linspace(count(1), count(end), indGpOptions.numActive)';
        indGpModel.kern.comp{1}.comp{1}.isNormalised = isNormalised;
        indGpModel.kern.comp{1}.comp{2}.isNormalised = isNormalised;
        indGpModel.kern.comp{1}.comp{2}.isStationary = isStationary;
        par = modelExtractParam(indGpModel);
        paramInd = paramNameRegularExpressionLookup(indGpModel, ['inverse width']);
        par(paramInd) = log(invWidthSim(i));
        if strcmp(indGpModel.kernType, 'sim')
            indGpModel.fix(end+1).index = 1;
            indGpModel.fix(end).value = log(invWidthSim(i));
            indGpModel.fix(end+1).index = 4;
            indGpModel.fix(end).value = log(invWidthSim(i));
        elseif strcmp(indGpModel.kernType, 'simwhite')
            indGpModel.fix(end+1).index = 1;
            indGpModel.fix(end).value = log(invWidthSim(i));
            indGpModel.fix(end+1).index = 2;
            indGpModel.fix(end).value = log(2^(1-i));
            indGpModel.fix(end+1).index = 4;
            indGpModel.fix(end).value = log(2^(1-i));
        else
            error('The smart initialisation for this model is not supported');
        end
        indGpModel = modelExpandParam(indGpModel, par);
        indGpModel = modelOptimise(indGpModel, 1, niter2);
        modelDisplay(indGpModel);
        paramsInit(j, :, i) = modelExtractParam(indGpModel);
    end
end

if isfield(options, 'typeLf')
    indGpOptions.kernType = 'simwhite';
    paramsInit2 = zeros(size(y{1}, 2), size(modelExtractParam(indGpModel), 2), numSimLf);
%     namesInit2 = cell(size(y{1}, 2), size(modelExtractParam(indGpModel), 2), numSimLf);
    for i = 1:numSimwhiteLf
        for j = 1:size(y{1}, 2)
            % Creation of the model
            XInd = X{1}(:, j);
            yInd = y{1}(:, j);
            indGpModel = multigpCreate(1, 1, XInd, yInd, indGpOptions);
            % Modification and fixing of some parameters
            indGpModel.X_u = linspace(count(1), count(end), indGpOptions.numActive)';
            indGpModel.kern.comp{1}.comp{1}.isNormalised = isNormalised;
            indGpModel.kern.comp{1}.comp{2}.isNormalised = isNormalised;
            indGpModel.kern.comp{1}.comp{2}.isStationary = isStationary;
            par = modelExtractParam(indGpModel);
            paramInd = paramNameRegularExpressionLookup(indGpModel, ['inverse width']);
            par(paramInd) = log(invWidthSimwhite(i));
            indGpModel.fix(end+1).index = 1;
            indGpModel.fix(end).value = log(invWidthSimwhite(i));
            indGpModel = modelExpandParam(indGpModel, par);
            % Training the model, displaying and saving the parameters
            indGpModel = modelOptimise(indGpModel, 1, niter2);
            modelDisplay(indGpModel);
            paramsInit2(j, :, i) = modelExtractParam(indGpModel);
        end
    end
end

return;
