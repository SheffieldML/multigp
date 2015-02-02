function model = spmultimodelExpandParam(model, params)

% SPMULTIMODELEXPANDPARAM Expand the parameters into a SPMULTIMODEL struct.
% FORMAT
% DESC expands the model parameters from a structure containing
% the information about a sparse multi model multi-output Gaussian process.
% ARG model : the model structure containing the information about
% the model.
% ARG params : a vector of parameters from the model.
% RETURN model : the model structure containing the information about
% the model updated with the new parameter vector.
%
% SEEALSO : spmultimodelCreate, modelExpandParam, modelExtractParam
%
% COPYRIGHT : Mauricio A. Alvarez, 2009

% MULTIGP

paramPart = real(params);

if isfield(model, 'fix')
    for i = 1:length(model.fix)
       paramPart(model.fix(i).index) = model.fix(i).value;
    end
end

startVal = 1;
endVal = model.kern.nParams;
kernParams = paramPart(startVal:endVal);

if length(kernParams) ~= model.kern.nParams
    error('kern parameter vector is incorrect length');
end

if model.varS
    switch model.isSpeed
        case 1
            model.kern = kernExpandParam(model.kern, kernParams);
        case 2
            model.kern = spmultimodelKernExpandParamVar(model.kern, kernParams);
        case 3
            model.kern = spmultimodelKernExpandParam(model.kern, kernParams);
        otherwise
            error('Unknown SPEED option')
    end
else
    model.kern = spmultimodelKernExpandParam(model.kern, kernParams);
end

% Check if there is a mean function.
if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
    startVal = endVal + 1;
    endVal = endVal + model.meanFunction.nParams;
    model.meanFunction = meanExpandParam(model.meanFunction, ...
        paramPart(startVal:endVal));
end

% Check if there is a gamma parameter.
if isfield(model, 'gamma') && ~isempty(model.gamma)
    startVal = endVal + 1;
    endVal = endVal + model.nlf;
    fhandle = str2func([model.gammaTransform 'Transform']);
    model.gamma = fhandle(paramPart(startVal:endVal), 'atox');
end

% Check if there is a beta parameter.
if isfield(model, 'beta') && ~isempty(model.beta)
    startVal = endVal + 1;
    endVal = endVal + model.nout;
    fhandle = str2func([model.betaTransform 'Transform']);
    model.beta = fhandle(paramPart(startVal:endVal), 'atox');
end

% Check if there is a gammas parameter.
if isfield(model, 'gammas') && ~isempty(model.gammas)
    startVal = endVal + 1;
    endVal = endVal + model.nout*model.nlf;
    fhandle = str2func([model.gammasTransform 'Transform']);
    model.gammas = fhandle(paramPart(startVal:endVal), 'atox');
end


% Substract the mean

model.m = spmultimodelComputeM(model);
model = spmultimodelKernCompute(model);
initt = cputime;
if model.isSpeed == 1 && isfield(model, 'subSpeed') && model.subSpeed == 2
    model = spmultimodelUpdateAD2(model);
else
    model = spmultimodelUpdateAD(model);
end
endtt = cputime - initt;
a = 1;

function kern = spmultimodelKernExpandParamVar(kern, params) 

startVal = 1;
endVal = 0;

for i =1:length(kern.comp)
    endVal = endVal + kern.comp{i}.nParams;
    kern.comp{i} = cmpndKernExpandParamRed(kern.comp{i},params(startVal:endVal));
    startVal = endVal + 1;
end

function kern = cmpndKernExpandParamRed(kern, params)

startVal = 1;
endVal = 0;
for i = 1:length(kern.comp)
  endVal = endVal + kern.comp{i}.nParams;
  kern.comp{i} = kernExpandParam(kern.comp{i}, params(1, startVal:endVal));
  startVal = endVal + 1;
end


function kern = spmultimodelKernExpandParam(kern, params) 

startVal = 1;
endVal = 0;

for i =1:length(kern.comp)
    endVal = endVal + kern.comp{i}.nParams;
    kern.comp{i} = kernExpandParam(kern.comp{i},params(startVal:endVal));
    startVal = endVal + 1;
end

function m = spmultimodelComputeM(model)

if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
    mu = meanCompute(model.meanFunction, mat2cell(ones(model.nout,1), ones(model.nout,1), 1), 0);
else
    mu = zeros(model.nout,1);
end
m = cell(model.nout,1);
for j=1:model.d
    m{j} = model.y{j} - mu(j);
    if model.bias(j)~=0
        m{j} = m{j} - model.bias(j);
    end
    if model.scale(j)~=1
        m{j} = m{j}/model.scale(j);
    end
end





