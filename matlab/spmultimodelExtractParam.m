function  [param, names] = spmultimodelExtractParam(model)

% SPMULTIMODELEXTRACTPARAM Extract the parameters of a SPMULTIMODEL struct.
% FORMAT
% DESC extracts the model parameters from a structure containing
% the information about a multi-output Gaussian process contained inside a
% spmultimodel structure.
% ARG model : the model structure containing the information about
% the model.
% RETURN params : a vector of parameters from the model.
%
% DESC does the same as above, but also returns parameter names.
% ARG model : the model structure containing the information about
% the model.
% RETURN params : a vector of parameters from the model.
% RETURN names : cell array of parameter names.
%
% SEEALSO : spmultimodelCreate, spmultimodelExpandParam, modelExtractParam
%
% COPYRIGHT : Mauricio A Alvarez, 2009

% MULTIGP

if model.varS
    if nargout > 1
        switch model.isSpeed
            case 1
                [paramKern, namesKern] = kernExtractParam(model.kern);
            case 2
                [paramKern, namesKern] = spmultimodelKernExtractParamVar(model.kern, 'component');
            case 3
                [paramKern, namesKern] = spmultimodelKernExtractParam(model.kern, 'component');
            otherwise
                error('Unknown SPEED option')
        end
    else
        switch model.isSpeed
            case 1
                paramKern = kernExtractParam(model.kern);
            case 2
                paramKern = spmultimodelKernExtractParamVar(model.kern, 'component');
            case 3
                paramKern = spmultimodelKernExtractParam(model.kern, 'component');
            otherwise
                error('Unknown SPEED option')
        end
    end
else
    % Extract params from the latent kernel
    if nargout > 1
        [paramKern, namesKern] = spmultimodelKernExtractParam(model.kern, 'component');
    else
        paramKern  = spmultimodelKernExtractParam(model.kern);
    end
end

if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
    if nargout>1
        [meanFuncParams, meanFuncParamNames] = meanExtractParam(model.meanFunction);
        for i = 1:length(meanFuncParamNames)
            meanFuncParamNames{i} = ['mean Func ' meanFuncParamNames{i}];
        end
    else
        meanFuncParams = meanExtractParam(model.meanFunction);
    end
else
    meanFuncParamNames = {};
    meanFuncParams =[];
end

if isfield(model, 'gamma') && ~isempty(model.gamma)
    fhandle = str2func([model.gammaTransform 'Transform']);
    gammaParams = fhandle(model.gamma, 'xtoa');    
    if nargout>1
        gammaParamNames = cell(model.nlf,1);
        for i = 1:length(gammaParams)
            gammaParamNames{i} = ['Gamma ' num2str(i)];
        end
    end
else
    gammaParamNames = {};
    gammaParams =[];
end

if isfield(model, 'beta') && ~isempty(model.beta)
    fhandle = str2func([model.betaTransform 'Transform']);
    betaParams = fhandle(model.beta, 'xtoa');    
    if nargout>1
        betaParamNames = cell(model.nout,1);
        for i = 1:length(betaParams)
            betaParamNames{i} = ['Beta ' num2str(i)];
        end
    end
else
    betaParamNames = {};
    betaParams =[];
end

if isfield(model, 'gammas') && ~isempty(model.gammas)
    fhandle = str2func([model.gammasTransform 'Transform']);
    gammasParams = fhandle(model.gammas, 'xtoa');    
    if nargout>1
        gammasParamNames = cell(model.nout*model.nlf,1);
        for i = 1:length(gammasParamNames)
            gammasParamNames{i} = ['Gammas ' num2str(i)];
        end
    end
else
    gammasParamNames = {};
    gammasParams =[];
end

param = [paramKern meanFuncParams gammaParams betaParams gammasParams];

% Fix the value of the parameters

if isfield(model, 'fix')
    for i = 1:length(model.fix)
        param(model.fix(i).index) = model.fix(i).value;
    end
end

if nargout > 1
    names = {namesKern{:}, meanFuncParamNames{:}, gammaParamNames{:}, ...
        betaParamNames{:}, gammasParamNames{:}};
end


function [params, names] = spmultimodelKernExtractParamVar(kern, type) 

if nargin<2
    type = [];
end

params = zeros(1, kern.nParams);
if nargout>1
    names = cell(1, kern.nParams);
end
startVal = 1;
endVal = 0;

for i =1:length(kern.comp)
    endVal = endVal + kern.comp{i}.nParams;
    if nargout > 1
        [params(startVal:endVal), names(startVal:endVal)] = ...
            cmpndKernExtractParamRed(kern.comp{i});
        for j=1:kern.comp{i}.nParams
            names{(startVal-1)+j} = [type ' ' num2str(i) ' ' names{(startVal-1)+j}];
        end        
    else
        params(startVal:endVal) = cmpndKernExtractParamRed(kern.comp{i});
    end
    startVal = endVal + 1;
end

function  [params, namesTemp] = cmpndKernExtractParamRed(kern)


params = zeros(1, kern.nParams);
if nargout > 1
  namesTemp = cell(1, kern.nParams);
end
startVal = 1;
endVal = 0;
storedTypes = cell(0);
for i = 1:length(kern.comp)
  endVal = endVal + kern.comp{i}.nParams;
  if nargout > 1
    [params(1, startVal:endVal), namesTemp(startVal:endVal)] = kernExtractParam(kern.comp{i});
    instNum = sum(strcmp(kern.comp{i}.type, storedTypes)) + 1;
    for j = startVal:endVal
      namesTemp{1, j} = [kern.comp{i}.type ' ' num2str(instNum) ' ' ...
                         namesTemp{1, j}];
    end
    storedTypes{end+1} = kern.comp{i}.type;
  else
    params(1, startVal:endVal) = kernExtractParam(kern.comp{i});
  end
    startVal = endVal + 1;
end


function [params, names] = spmultimodelKernExtractParam(kern, type) 

if nargin<2
    type = [];
end

params = zeros(1, kern.nParams);
if nargout>1
    names = cell(1, kern.nParams);
end
startVal = 1;
endVal = 0;

for i =1:length(kern.comp)
    endVal = endVal + kern.comp{i}.nParams;
    if nargout > 1
        [params(startVal:endVal), names(startVal:endVal)] = ...
            kernExtractParam(kern.comp{i});
        for j=1:kern.comp{i}.nParams
            names{(startVal-1)+j} = [type ' ' num2str(i) ' ' names{(startVal-1)+j}];
        end        
    else
        params(startVal:endVal) = kernExtractParam(kern.comp{i});
    end
    startVal = endVal + 1;
end

