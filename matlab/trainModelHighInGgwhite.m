function  [logModel, model] =  trainModelHighInGgwhite(X, y, numActive, nout, ... 
    iters, noisePerOutput, inverseWidth, sensitivity,...
    nVIKs, paramInit, methodInit, methodInitMGP)

% TRAINMODELHIGHINGGWHITE Trains the DTC VAR approx for GGWHITE kernel
% FORMAT
% DESC Trains a multigp model with DTC VAR approx with fixed parameters
% ARG X : set of training inputs
% ARG y : set of training observations
% ARG numActive : number of pseudo points used for the approx.
% ARG nout : number of outputs
% ARG iters : number of iterations for the optimization
% ARG noisePerOutput : noise value for the parameters of the Output kernels
% ARG inverseWidth : inverse width for the parameters of the Output kernels
% ARG sensitivity : sensitivity value for the parameters of the Output
% kernels
% ARG nVIKs : number of inducing kernels used for the approx (Default = numActive)
% ARG paramInit : initial value for the inverse width of the VIKs
% ARG methodInit : method to initialize the positions of the pseudo points.
% ARG methodInitMGP : method to initialize the multi GP (Default = optimize
% all inducing and pseudo inputs).
% RETURN logModel : the bound given by the approx.
% RETURN model : trained model for the particular setup
%
% SEE ALSO : multigpCreate, muligpOptimise
%
% COPYRIGHT : Mauricio Alvarez, 2009

% MULTIGP 

if nargin <12
    if nargin < 11
        if nargin < 10
            if nargin < 9
                nVIKs = [];
            end
            paramInit = 1e2;
        end
        methodInit = 'randomComplete';
    end
    methodInitMGP = false;   
end

% if ~isempty(nVIKs)
%     helperCreateNames(size(X{1},2), numActive);
% end

%inverseWidth = [10 5 5 100 30 200 400 20 1 40]; % First option
%inverseWidth = [4   5.5 5 7 3   2 4 4.5 3.5 6.5]; % Second option
%sensitivity =  [3.5 5   2 2 3.5 4 4 8   5   5];
invNoisePerOutput = 1./noisePerOutput;

options = multigpOptions('dtc');
options.kernType = 'simwhite';
options.optimiser = 'scg';
options.nlf = 1;
options.initialInducingPositionMethod = methodInit;
options.numActive = numActive;
options.beta = 1e-4*ones(1, nout);
options.fixInducing = methodInitMGP;
options.tieOptions.selectMethod = 'nofree';
options.includeNoise = 1;
%options.isArd = 1;
options.nVIKs = nVIKs;

q = size(X{1},2);
d = size(y, 2);

% Creates the model
model = multigpCreate(q, d, X, y, options);
% Fix parameters
model = fixParamsMultigp(model, options, inverseWidth, sensitivity, ...
    invNoisePerOutput);
model = initInducingMultigp(model, options, inverseWidth, paramInit);
displayGradchek = 1;
% Trains the model
if methodInitMGP
   itersInit = 5000;
 else
   itersInit = iters;
end
model = multigpOptimise(model, displayGradchek, itersInit);

if methodInitMGP
    modelInit = model;
    options.fixInducing = false;
    model = multigpCreate(q, d, X, y, options);
    model.X = modelInit.X;
    % Fix parameters
    model = fixParamsMultigp(model, options, inverseWidth, sensitivity, ...
			     noisePerOutput);
    paramsInit = modelExtractParam(modelInit);    
    paramInitVIK = exp(paramsInit(1:options.nVIKs));
    model = initInducingMultigp(model, options, inverseWidth, paramInitVIK);
    params = modelExtractParam(model);    
    paramInd = paramNameRegularExpressionLookup...
      (modelInit, '.* white 1 variance');
    paramInd2 = paramNameRegularExpressionLookup...
      (model, '.* white 1 variance');
    params(paramInd2) = paramsInit(paramInd);
    model = modelExpandParam(model, params);    
    % Trains the model 
    displayGradchek = 1;   
    model = multigpOptimise(model, displayGradchek, iters);
end
logModel = multigpLogLikelihood(model); 
