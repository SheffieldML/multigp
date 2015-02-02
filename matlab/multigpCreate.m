function  model = multigpCreate(q, d, X, y, options)

% MULTIGPCREATE creates a multi output GP based on a convolution.
% Creates a multiple output Gaussian process using the idea of
% the convolution process. The multigp model could be either FULL in the
% sense that the full covariance of the model is employed or be SPARSE in
% the sense that a low rank approximation is employed. The outputs of the
% model are generated according to the convolution operation
%
%       f_q(x) = \int_{-\infty}^{\infty}k_q(x-z)u(z)dz
%
% where f_q(x) is an output function, k_q(x) is a smoothing kernel and u(z)
% is a latent function, which is considered as a Gaussian process. In
% principle, more latent functions could be employed as well as the
% addition of an independent process is desirable. In that case, the output
% function is now given as
%
%      y_q(x) = f_q(x) + w_q(x)
%
% where w_q(x) corresponds to the independent process.
%
% The code also includes the case when k_q(x-z)=a_q\delta(x-z), which
% corresponds to the linear model of corregionalization (LMC).
%
% FORMAT
% DESC returns a structure for the multiple output Gaussian process model.
% RETURN model : the structure for the multigp model
% ARG q : input dimension size.
% ARG d : output dimension size.
% ARG Y : set of training observations
% ARG X : set of training inputs
% ARG options : contains the options for the MULTIGP model, which includes
% if the model is approximated or full, the number of latent functions, the
% number of output functions.
%
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008
%
% MODIFICATIONS : David Luengo, 2009

% MULTIGP

if iscell(X)
    if size(X{end}, 2) ~= q
        error(['Input matrix X does not have dimension ' num2str(q)]);
    end
else
    if size(X, 2) ~= q
        error(['Input matrix X does not have dimension ' num2str(q)]);
    end
end
if iscell(y)
    % ensure it is a row vector of cells.
    y = y(:)';
    if size(y, 2) ~= d
        error(['Target cell array Y does not have dimension ' num2str(d)]);
    end
    for i = 1:size(y, 2)
        if(size(y{i}, 2)>1)
            error('Each element of the cell array should be a column vector.')
        end
    end
else
    if size(y, 2)~=d
        error(['Target matrix Y does not have dimension ' num2str(d)]);
    end
end


model.type = 'multigp';
model.q = q;

% Ask if this is part of a multimodel structure
if isfield(options, 'isSubmodel') && ~isempty(options.isSubmodel)
    model.isSubmodel = options.isSubmodel;
end

% Short hand for the kernel type
model.kernType = options.kernType;

% Number of latent functions
model.nlf = options.nlf;
if isfield(options, 'typeLf') && ~isempty(options.typeLf)
    model.typeLf = options.typeLf;
end

% Number of output functions
model.d = d;

switch options.approx
    case 'ftc'
        model.nout = d - options.nlf;
    case {'dtc','fitc','pitc', 'dtcvar'}
        model.nout = d;
end

model.approx = options.approx;

% Set up default scale and bias for outputs
if isfield(options, 'scale') && ~isempty(options.scale)
    model.scale = options.scale;
else
    model.scale = ones(1, model.d);
end
if isfield(options, 'bias') && ~isempty(options.bias)
    model.bias = options.bias;
else
    model.bias = zeros(1, model.d);
end

% Initialization of the inputs for the model


switch model.approx
    case 'ftc'
        model.X = X;
        model.y = [];
        for i = 1:length(y)
            model.y = [model.y; y{i}];
        end
        model.N = size(model.y,1);        
    case {'dtc','fitc','pitc', 'dtcvar'}
        if numel(options.numActive) ~= options.nlf
            if numel(options.numActive) == 1,
                options.numActive = options.numActive*ones(1, options.nlf);
            else
                options.numActive = options.numActive(1)*ones(1, options.nlf);
                warning(['The number of latent functions does not match'...
                    ' the number of elements provided in options.numActive']);
            end
        elseif isfield(options, 'tieInducing') && ~isempty(options.tieInducing) && ...
                    options.tieInducing
                options.numActive = options.numActive(1)*ones(1, options.nlf);
                warning(['Since options.tieInducing is true '...
                    ' only the first element in options.numActive is choosen']);        
        end
        if isfield(options, 'tieInducing') && ~isempty(options.tieInducing) && ...
                options.tieInducing
            effPS = 1;
        else
            effPS = options.nlf;
        end
        for i = 1: effPS,
            numActive = options.numActive(i);
            posX = zeros(numActive, q);
            switch options.initialInducingPositionMethod
                case 'espaced'
                    % Especially useful for 1D case. It allows to move the
                    % pseudo-inputs a further to the original input range
                    for j = 1:q,
                        factor = 0.1;
                        med = (max(X{1}(:,j)) - min(X{1}(:,j)))/2;
                        posX(:,j) = linspace(min(X{1}(:,j)) - factor*med, ...
                            max(X{1}(:,j)) + factor*med, numActive)';
                    end
                case 'espacedInRange'
                    % Especially useful for 1D case. It restricts the
                    % initial position of pseudo-inputs to be within the input range
                    for j = 1:q,
                        posX(:,j) = linspace(min(X{1}(:,j)), max(X{1}(:,j)), numActive)';
                    end
                    
                case 'nonrandomDataIsotopic'
                    % The initial positions of pseudo-inputs are taken from
                    % the data. 
                    if size(X{1},1) >= numActive,
                        totX = cell2mat(X(1));
                    else
                        totX = cell2mat(X');
                    end                    
                    posX = totX(1:numActive,:);
                    
                case 'randomDataIsotopic'
                    % The initial positions of pseudo-inputs are taken from
                    % the data. In the isotopic case, since all inputs all
                    % equal for each output, we use only X(1) to select the
                    % positions, as long as the number of numActive is less
                    % than the size of inputs.
                    if size(X{1},1) >= numActive,
                        totX = cell2mat(X(1));
                    else
                        totX = cell2mat(X');
                    end
                    index = randperm(size(totX,1));
                    posX = totX(index(1:numActive),:);
                case 'randomDataHeterotopic'
                    totX = cell2mat(X');
                    index = randperm(size(totX,1));
                    posX = totX(index(1:numActive),:);
                case 'randomComplete'
                    posX = 0.5*rand(numActive, q);
                case 'fixIndices'
                    posX = X{1}(options.fixIndices,:);
                case 'kmeansIsotopic'
                    posX = kmeanlbg(X{1},numActive);
                case 'kmeansHeterotopic'
                    totX = cell2mat(X');
                    posX = kmeanlbg(totX,numActive);
                otherwise
                    error('This is not valid initialization method for the input variables');
            end
            model.X{i,1} = posX;
        end
        if isfield(options, 'tieInducing') && ~isempty(options.tieInducing) && ...
                options.tieInducing
            for i=2:model.nlf
               model.X{i, 1} = model.X{1}; 
            end
        end
        model.y = [];
        for i = 1:length(y)
            model.y = [model.y; y{i}];
            model.X{i + options.nlf,1} = X{i};
        end
        model.N = size(model.y,1);
        X = model.X;
    otherwise
        error('Unknown model approximation')
end

if strcmp(model.kernType, 'lmc')
    % In the LMC isotopic, only one of the input data is needed.
    if isfield(options, 'isIsotopic') && ~isempty(options.isIsotopic) ...
            && options.isIsotopic
        model.X = X(1:model.nlf+1);
        model.isIsotopic = true;
    else
        model.isIsotopic = false;
    end
    if isfield(options, 'rankCorregMatrix') && ~isempty(options.rankCorregMatrix) ...
            && options.rankCorregMatrix
        model.rankCorregMatrix = options.rankCorregMatrix;
    else
        model.rankCorregMatrix = 1;
    end    
end

model.includeInd = options.includeInd;
model.tieIndices = options.tieOptions.tieIndices;
model.includeNoise = options.includeNoise;
kernType = cell(1,sum(options.nlf) + options.includeNoise + ...
    options.includeInd);
cont = 0;

if isfield(options, 'isSpeedUp') && ~isempty(options.isSpeedUp)
    model.isSpeedUp = options.isSpeedUp;
    switch options.isSpeedUp
        case 1
            % This trick only works for outputs that have the same input space
            if strcmp(model.approx, 'ftc')
                upperBoundKernel = model.nlf+2;
            else
                upperBoundKernel = 2;
            end
            for i = 1:options.nlf
                cont = cont + 1;
                kernType{cont} = multigpKernComposer(options.kernType, upperBoundKernel, model.nlf, model.approx, i, options);
            end
            % To include noise
            if model.includeNoise
                cont = cont + 1;
                if strcmp(model.kernType, 'lmc')
                    kernType{cont} = multigpKernComposer('whiteblock', upperBoundKernel, model.nlf, model.approx, [], options);
                else
                    kernType{cont} = multigpKernComposer('white', upperBoundKernel, model.nlf, model.approx);
                end
            end
            kernTemplate = kernCreate(X(1:model.nlf+2), {'cmpnd', kernType{:}});
            % Make a corrections
            kernTemplateEach = cell(1, model.nlf+1);
            crossFunction = cell(1, model.nlf+1);
            for i=1:model.nlf+1
                kernTemplateEach{i} = kernTemplate.comp{i}.comp{end};
                crossFunction{i} = kernTemplate.comp{i}.block{end}.cross{end};
            end
            kern = kernTemplate;
            kern.nParams = 0;
            for i =1:model.nlf+1
                for j=3:model.nout
                    kern.comp{i}.comp{model.nlf+j} = kernTemplateEach{i};
                    kern.comp{i}.nParams = kern.comp{i}.nParams + kernTemplateEach{i}.nParams;
                    kern.comp{i}.numBlocks = kern.comp{i}.numBlocks + 1;
                    kern.comp{i}.block{model.nlf+j}.cross = [kern.comp{i}.block{model.nlf+j-1}.cross {crossFunction{i}}];
                    kern.comp{i}.block{model.nlf+j}.transpose = [kern.comp{i}.block{model.nlf+j-1}.transpose 0];
                    kern.comp{i}.inputDimension(model.nlf+j) = kernTemplate.comp{i}.inputDimension(end);
                    kern.comp{i}.diagBlockDim{model.nlf+j} = kernTemplate.comp{i}.diagBlockDim{end};
                end
                kern.comp{i}.paramGroups = speye(kern.comp{i}.nParams);
                kern.nParams = kern.nParams + kern.comp{i}.nParams;
            end
            kern.inputDimension = kern.comp{1}.inputDimension;
            kern.numBlocks = kern.comp{1}.numBlocks;
            kern.paramGroups = speye(kern.nParams);
            model.kern = kern;
        case 2
            % Creates a global kernel. This is a more restrictive option
            % in the multigp toolbox, but useful for some very
            % repetitive cases. WARINING : so far, it only works for cases
            % in which there is no need to optimize the inducing variables.
            % It is certainly less flexible than the standard multigp
            % version, but much faster.
            kern.type = [options.kernType 'global'];
            kern.inputDimension = q;
            if isfield(options, 'kern') && ~isempty(options.kern)
                kern.options = options.kern;             
            end
            kern.options.nlf = options.nlf;
            kern.options.nout = model.nout;
            kern.options.approx = model.approx;                     
            kern = kernParamInit(kern);
            kernType = multigpKernComposer(options.kernType, 2, 1, 'ftc', 1, options);            
            kern.template.latent = kernCreate(X{model.nlf+1}, kernType{2});
            kern.template.output = kernCreate(X{model.nlf+1}, kernType{3});
            if isfield(options, 'kern') && ~isempty(options.kern)
                typeKernLat = kernType{2}{3};
                typeKernOut = kernType{3}{3};
            else
                typeKernLat = kernType{2};
                typeKernOut = kernType{3};
            end
            kern.funcNames.computeLat = str2func([typeKernLat  'KernCompute']);
            kern.funcNames.computeDiagLat = str2func([typeKernLat  'KernDiagCompute']);
            kern.funcNames.computeOut = str2func([typeKernOut  'KernCompute']);
            kern.funcNames.computeDiagOut = str2func([typeKernOut  'KernDiagCompute']);
            kern.funcNames.computeOutCrossOut = str2func([typeKernOut 'X' typeKernOut 'KernCompute']);
            kern.funcNames.computeOutCrossLat = str2func([typeKernOut 'X' typeKernLat 'KernCompute']);
            kern.funcNames.gradientLat = str2func([typeKernLat 'KernGradient']);
            if isfield(options, 'useKernDiagGradient') && options.useKernDiagGradient ...
                   && (strcmp(options.approx, 'fitc') || strcmp(options.approx, 'dtcvar'))
                kern.funcNames.gradientOut = str2func([typeKernOut 'KernDiagGradient']);
            else                
                kern.funcNames.gradientOut = str2func([typeKernOut 'KernGradient']);
            end
            kern.funcNames.gradientOutCrossOut = str2func([typeKernOut 'X' typeKernOut 'KernGradient']);
            kern.funcNames.gradientOutCrossLat = str2func([typeKernOut 'X' typeKernLat 'KernGradient']);
            kern.funcNames.extractLat = str2func([typeKernLat 'KernExtractParam']);
            kern.funcNames.extractOut = str2func([typeKernOut 'KernExtractParam']);
            kern.funcNames.transferParamLat = str2func([typeKernLat 'KernParamTransfer']);
            kern.funcNames.transferParamOut = str2func([typeKernOut 'KernParamTransfer']);
            kern.funcNames.transferGradLat = str2func([typeKernLat 'KernGradTransfer']);
            kern.funcNames.transferGradOut = str2func([typeKernOut 'KernGradTransfer']);
            kern.funcNames.gradInit = str2func([kern.type 'KernGradInit']);
            kern.funcNames.gradCat = str2func([kern.type 'KernGradCat']);
            model.kern = kern;
            model.kernType = kern.type;
            model.kern.paramGroups = speye(model.kern.nParams);            
        otherwise
            error('Unknown SPEED UP option')
    end
else
    % This checks if there are different clases of latente forces
    if isfield(options, 'typeLf') && ~isempty(options.typeLf)
        if options.nlf ~= sum(options.typeLf)
            error('The number of latent functions per type is different to the total number of latent functions')
        end
        for i = 1:options.typeLf(1)
            cont = cont + 1;
            kernType{cont} = multigpKernComposer(options.kernType, model.d, model.nlf, model.approx, i);
        end
        for i = 1:options.typeLf(2)
            cont = cont + 1;
            kernType{cont} = multigpKernComposer([options.kernType 'white'], model.d, model.nlf, model.approx, i+options.typeLf(1));
        end
    else
        for i = 1:options.nlf
            cont = cont + 1;
            kernType{cont} = multigpKernComposer(options.kernType, model.d, model.nlf, model.approx, i, options);
        end
    end
    % To include independent kernel
    if model.includeInd
        cont = cont + 1;
        kernType{cont} = multigpKernComposer('rbf', model.d, model.nlf, model.approx);
    end
    % To include noise
    if model.includeNoise
        cont = cont + 1;
        if strcmp(model.kernType, 'lmc')
            kernType{cont} = multigpKernComposer('whiteblock', model.d, model.nlf, model.approx, [], options);
        else
            kernType{cont} = multigpKernComposer('white', model.d, model.nlf, model.approx);
        end
    end
    model.kern = kernCreate(X, {'cmpnd', kernType{:}});    
end


if ~isfield(options, 'isSpeedUp') || (isfield(options, 'isSpeedUp') && options.isSpeedUp ~= 2)    
    % These are additional options for some kernels that must be handled outside
    % kernCreate. For example, for the GG kernel, the user can choose ARD or
    % not.
    fhandle = [model.kernType 'MultigpKernOptions'];
    if exist(fhandle, 'file')
        fhandle = str2func(fhandle);
        model = fhandle(model, options);
    end
end

if isfield(options, 'optimiser') && ~isempty(options.optimiser)
    model.optimiser = options.optimiser;
end

% NEIL: Learning of scales hasn't been included, although it should be.
model.learnScales = options.learnScales;
model.scaleTransform = optimiDefaultConstraint('positive');

switch model.approx
    case {'dtc','fitc','pitc', 'dtcvar'}
        model = spmultigpCreate( model, options);
end

model.nParams = 0;
% Extract top level parameters from kernels.
model = multigpUpdateTopLevelParams(model);

% Set up a mean function if one is given.
if isfield(options, 'meanFunction') && ~isempty(options.meanFunction)
    if isstruct(options.meanFunction)
        model.meanFunction = options.meanFunction;
    else
        if ~isempty(options.meanFunction)
            model.meanFunction = meanCreate(q, model.nout, X, y, options.meanFunctionOptions);
        end
    end
    model.nParams = model.nParams + model.meanFunction.nParams;
end

%if ~isfield(options, 'isSpeedUp') || (isfield(options, 'isSpeedUp') && options.isSpeedUp ~= 2)
% Tie options according to the particular kernel employed

fhandle = [model.kernType 'MultigpTieParam'];
if exist(fhandle, 'file')
    fhandle = str2func(fhandle);
    tieInd = fhandle(model, options);
else
    if ~isfield(options, 'isSpeedUp') || (isfield(options, 'isSpeedUp') && options.isSpeedUp ~= 2)
        error('Function for tying parameters for this kernel type not implemented yet')
    else
        warning('multigpCreate:NoTyingFunctionGlobalKernel','Function for tying parameters for this kernel type not implemented yet')
    end
end

model.nParams = model.nParams + model.kern.nParams;

% Creates the noise model (Default model: one beta per output)
% switch model.approx
%     case 'ftc'
%         if isfield(options, 'noiseOpt') && ~isempty(options.noiseOpt) && options.noiseOpt
%             model.noiseOpt = options.noiseOpt;            
%             if options.includeNoise
%                 warning('multigpCreate:additionalNoiseParameter', 'Additional noise parameter added')                
%             end
%             model.variance = log(exp(-2));
%             model.varianceTransform = optimiDefaultConstraint('positive');
%             model.nParams = model.nParams + 1;
%             % Looks for nRepeats in options
%             if isfield(options, 'nRepeats') && ~isempty(options.nRepeats)
%                 model.nRepeats = options.nRepeats;
%             else
%                 error('Must provide the number of repeats')
%             end
%         else
%             warning('multigpCreate:missingNoiseModel','Noise model has not been included')
%         end
%     case {'dtc','fitc','pitc', 'dtcvar'}
if isfield(options, 'gamma') && ~isempty(options.gamma)
    if size(options.gamma,2) == model.nlf
        model.gamma = options.gamma;
    else
        model.gamma = options.gamma*ones(1,model.nlf);
    end
    model.gammaTransform =  optimiDefaultConstraint('positive');
    model.nParams = model.nParams + model.nlf;
end
if isfield(options, 'beta') && ~isempty(options.beta)
    if size(options.beta,2) == model.nout
        model.beta = options.beta;
    else
        model.beta = options.beta*ones(1,model.nout);
    end
    model.betaTransform =  optimiDefaultConstraint('positive');
    model.nParams = model.nParams + model.nout;
    if isfield(options, 'noiseOpt') && ~isempty(options.noiseOpt)
        model.noiseOpt = options.noiseOpt;
        switch options.noiseOpt
            case 0
                % Tie the noise parameters
                index = paramNameRegularExpressionLookup(model, 'Beta .*');
                tieInd{end+1} = index;
            case 1
                index = paramNameRegularExpressionLookup(model, 'Beta .*');
                tieInd{end+1} = index;
                % Looks for nRepeats in options
                if isfield(options, 'nRepeats') && ~isempty(options.nRepeats)
                    model.nRepeats = options.nRepeats;
                else
                    error('Must provide the number of repeats')
                end
            case 2
                % This case assumes a noise value for each entry.
                model.nParams = model.nParams - model.nout;
                model.beta = cell(model.nout,1);
                betaParams = 0;
                for k=1:model.nout
                    betaParams = betaParams + size(model.X{model.nlf+k},1);
                    model.beta{k} = options.beta(1)*ones(size(model.X{model.nlf+k},1),1);
                end
                model.nParams = model.nParams + betaParams;
            case 3
                % This case assumes a value for the variance is
                % provided
                model.nParams = model.nParams - model.nout;
                if isfield(options, 'yvar') && ~isempty(options.yvar)
                    model.beta = options.yvar;
                else
                    error('For option.noiseOpt=3, values for variances must be provided')
                end
                model = rmfield(model, 'betaTransform');
        end
    end
end
%[model, tieInd] =  noiseModelCreate(model, options, tieInd);
if ~strcmp(model.approx, 'ftc')
    if ~options.fixInducing
        model.nParams = model.nParams + sum(model.k)*model.q;
    end
end
%end

% Fix parameters options according to the particular kernel employed
fhandle = [model.kernType 'MultigpFixParam'];
if exist(fhandle, 'file')
    fhandle = str2func(fhandle);
    model = fhandle(model, options);
% else
%     warning('multigp:FixedParam','Function for fixing parameters for this kernel type not implemented yet')
end
if ~isfield(options, 'isSpeedUp') || (isfield(options, 'isSpeedUp') && options.isSpeedUp ~= 2)
    if isfield(options, 'fixInducing') && ~options.fixInducing && ...
            isfield(options, 'tieInducing') && ~isempty(options.tieInducing) && ...
            options.tieInducing && options.nlf > 1
        tieInd(end+1:end+(model.k(1)*model.q)) = spmultigpTiePseudoInputs(model);
    end
    model = modelTieParam(model, tieInd);
else
    if exist('tieInd', 'var') && ~isempty(tieInd)
        model = spmultimodelTieParam(model, tieInd);
    end
end

if isfield(options, 'useKernDiagGradient') && options.useKernDiagGradient
    model.useKernDiagGradient = options.useKernDiagGradient;    
end

params = modelExtractParam(model);
model = modelExpandParam(model, params);

model.alpha = [];
end
