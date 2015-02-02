function model = spmultimodelCreate( q, d, X, y, options )

% SPMULTIMODELCREATE Convolved multigp model with connectivity matrix.
% Creates a multigp model with given connectivity matrix. This function is
% based on MULTIGPCREATE. It appeals to a more flexible structure which
% saves memory compared to MULTIGP and allows the specification in a direct
% way of a connectivity matrix. This connectivity matrix must be supplied
% in the structure OPTIONS with the name connect. It is a matrix with the
% number of rows equal to the number of outputs and the number of columns
% eqaul to the number of latent functions. Each entry is either a zero or a
% one indicating the relation between that particular output and that
% particular latent force. The kernel in the model is a cell where each
% enty is at the same time a multiKern structure, of the form {{rbf},
% {sim}..{sim}} where the 'rbf' part refers to the particular latent
% function and the 'sim' part to the number of outputs asscociated to that
% particular latent force according to the connectivity matrix. Furthermore
% the important information about the connetivity matrix is stored in the
% cells 'indLatGivenOut' and 'indOutGivenOut'. In 'indLatGivenOut', ech
% entry provides the latent functions associated to the given output and in
% 'indOutGivenLat' each entry provides the outputs associated to a
% particular latent function.
%
% FORMAT
% DESC returns a structure for the spmultimodel output Gaussian process model.
% RETURN model : the structure for the spmultimodel multigp model
% ARG q : input dimension size.
% ARG d : output dimension size.
% ARG Y : set of training observations
% ARG X : set of training inputs
% ARG options : contains the options for the model, which includes
% if the model is approximated or full, the number of latent functions, the
% number of output functions, the connectivity matrix.
%
% SEE ALSO : multigpCreate.m
%
% COPYRIGHT : Mauricio A. Alvarez 2009, 2010

% MULTIGP



model.type = 'spmultimodel';

if isfield(options, 'varS') && ~isempty(options.varS)
    if ~strcmp(options.approx, 'dtcvar')
        error('Learning the sensitivities is only possible for dtcvar approximation.')
    end
    model.varS = options.varS;
else
    model.varS = false;
end

if isfield(options, 'connect') && ~isempty(options.connect)
    model.connect = options.connect;
elseif ~model.varS
    error('For this type of model you need to provide a connectivity matrix. Otherwise, use the multigp model')
end

model.kernType = options.kernType;
model.nout = d;
model.d = d;
model.q = q;
model.nlf =options.nlf;
model.approx = options.approx;
model.optimiser = options.optimiser;

if isfield(options, 'scale') && ~isempty(options.scale)
    model.scale = options.scale;
else
    model.scale = ones(1, model.nout);
end
if isfield(options, 'bias') && ~isempty(options.bias)
    model.bias = options.bias;
else
    model.bias = zeros(1, model.nout);
end

switch model.approx
    case 'ftc'
    % Not Implemented yet        
    case {'fitc','pitc', 'dtcvar'}
        if numel(options.numActive) ~= options.nlf
            if numel(options.numActive) == 1,
                options.numActive = options.numActive*ones(1, options.nlf);
            else
                options.numActive = options.numActive(1)*ones(1, options.nlf);
                warning('spmultimodelCreate:noMatchingLF',['The number of latent functions does not match'...
                    ' the number of elements provided in options.numActive']);
            end
        end
        for i = 1: options.nlf,
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
            model.latX{i,1} = posX;
        end        
        for i = 1:length(y)            
            model.outX{i,1} = X{i};
            model.y{i} = y{i};
        end  
    otherwise
        error('Unknown model approximation')
end

model.tieIndices = options.tieOptions.tieIndices;
model.includeNoise = options.includeNoise;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model.isSpeed = options.isSpeed;
if model.isSpeed == 1
    model.subSpeed = options.subSpeed;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the dimension of each output in model, which will be useful for
% computing sizes and things like that[param, names] = spmultimodelExtractParam(model)


for k =1:model.nout,
   model.sizeX(k) = size(X{k},1);
end


if ~model.varS
    % Summarizes statistics contained in model.connect. When it is variational
    % these statistics are  meaningless.
    for k = 1:model.nout,
        model.indLatGivenOut{k,1} = find(model.connect(k, :) == 1);
        model.numLatGivenOut(k,1) = length(model.indLatGivenOut{k,1});
    end

    for k = 1:model.nlf,
        model.indOutGivenLat{k,1} = find(model.connect(:, k) == 1);
        model.numOutGivenLat(k,1) = length(model.indOutGivenLat{k,1});
    end
    model.latConnect    = sparse(model.connect'*model.connect);
    model.latConnectInv = inv(model.latConnect);
end
% MAURICIO : I fixed the number of points to be equal for all latent forces
% but be aware that might not be the case
model.k = options.numActive(1);
%model.k = options.numActive;

numParams = 0;
if ~model.varS
    % Creates a kernel that keeps the latent structures. If there's connect
    % information, it uses the stattistics to refer to the appropriate inputs.
    for i=1:model.nlf,
        fprintf('Creating component : %d\n',i)
        whichOutput = model.indOutGivenLat{i};
        localX = cell(1+model.numOutGivenLat(i),1);
        localX{1,1} = model.latX{i};
        for j = 1:model.numOutGivenLat(i),
            localX{1+j,1} = model.outX{whichOutput(j)};
        end
        if isfield(options, 'kernOptions') && ~isempty(options.kernOptions)
            kernType = multigpKernComposer(options.kernType, model.numOutGivenLat(i), 1, model.approx, 1, options.kernOptions);
        else
            kernType = multigpKernComposer(options.kernType, model.numOutGivenLat(i), 1, model.approx, 1);
        end
        model.kern.comp{i,1} = kernCreate(localX, kernType);
        numParams = numParams + model.kern.comp{i}.nParams;
    end
else
    model.template = zeros(sum(model.sizeX), model.nout);
    localX = cell(1+model.nout,1);
    localX{1,1} = model.latX{1};
    startVal = 1;
    endVal = 0;
    for j = 1:model.nout,
        endVal = endVal + model.sizeX(j);
        model.template(startVal:endVal, j) = ones(model.sizeX(k),1);
        localX{1+j,1} = model.outX{j};
        startVal = endVal + 1;
    end
    switch model.isSpeed
        case 1
            kern.type = [options.kernType 'global'];          
            kern.options.nlf = options.nlf;
            kern.options.nout = model.nout;
            kern.options.approx = model.approx;
            kern = kernParamInit(kern);          
            kernType = multigpKernComposer(options.kernType, 1, 1, model.approx, 1, options.kernOptions);
            kern.template.latent = kernCreate(localX{1}, kernType{2});
            kern.template.output = kernCreate(localX{2}, kernType{3});
            kern.funcNames.computeLat = str2func([kernType{2}{3}  'KernCompute']);
            kern.funcNames.computeOut = str2func([kernType{3}{3}  'KernDiagCompute']);
            kern.funcNames.computeCross = str2func([kernType{3}{3} 'X' kernType{2}{3} 'KernCompute']);
            kern.funcNames.gradientLat = str2func([kernType{2}{3} 'KernGradient']);
            kern.funcNames.gradientOut = str2func([kernType{3}{3} 'KernGradient']);
            kern.funcNames.gradientCross = str2func([kernType{3}{3} 'X' kernType{2}{3} 'KernGradient']);
            kern.funcNames.extractLat = str2func([kernType{2}{3} 'KernExtractParam']);
            kern.funcNames.extractOut = str2func([kernType{3}{3} 'KernExtractParam']);
            model.kern = kern;
            model.kernType = kern.type;
            model.kern.paramGroups = speye(model.kern.nParams);
            numParams = model.kern.nParams;
        case 2
            kernType = multigpKernComposer(options.kernType, 1, 1, model.approx, 1, options.kernOptions);
            model.kernFuncNames.computeLat = str2func([kernType{2}{3}  'KernCompute']);
            model.kernFuncNames.computeOut = str2func([kernType{3}{3}  'KernDiagCompute']);
            model.kernFuncNames.computeCross = str2func([kernType{3}{3} 'X' kernType{2}{3} 'KernCompute']);
            model.kernFuncNames.gradientLat = str2func([kernType{2}{3} 'KernGradient']);
            model.kernFuncNames.gradientOut = str2func([kernType{3}{3} 'KernGradient']);
            model.kernFuncNames.gradientCross = str2func([kernType{3}{3} 'X' kernType{2}{3} 'KernGradient']);
            model.kernFuncNames.extractLat = str2func([kernType{2}{3} 'KernExtractParam']);
            model.kernFuncNames.extractOut = str2func([kernType{3}{3} 'KernExtractParam']);
            kernTemplateForce = kernCreate(localX{1}, kernType{2});
            kernTemplateOutput = kernCreate(localX{2}, kernType{3});
            model.kernFuncNames.transformLat.type = str2func([kernTemplateForce.transforms.type 'Transform']);
            model.kernFuncNames.transformLat.index = kernTemplateForce.transforms.index;
            model.kernFuncNames.transformOut.type = str2func([kernTemplateOutput.transforms.type 'Transform']);
            model.kernFuncNames.transformOut.index = kernTemplateOutput.transforms.index;
            fprintf('Creating component : 1\n')
            model.kern.comp{1,1}.comp{1,1} = kernTemplateForce;
            numParams = numParams + model.kern.comp{1}.comp{1}.nParams;
            for i=2:model.nout+1
                model.kern.comp{1}.comp{i} = kernTemplateOutput;
                numParams = numParams + model.kern.comp{1}.comp{i}.nParams;
            end
            model.kern.comp{1}.nParams = numParams;
            for i=2:model.nlf,
                fprintf('Creating component : %d\n',i)
                model.kern.comp{i,1} = model.kern.comp{1,1};
                numParams = numParams + model.kern.comp{i}.nParams;
            end
        case 3
            kernType = multigpKernComposer(options.kernType, model.nout, 1, model.approx, 1, options.kernOptions);
            fprintf('Creating component : 1\n')
            model.kern.comp{1,1} = kernCreate(localX, kernType);
            numParams = numParams + model.kern.comp{1}.nParams;
            for i=2:model.nlf,
                fprintf('Creating component : %d\n',i)
                model.kern.comp{i,1} = model.kern.comp{1,1};
                numParams = numParams + model.kern.comp{i}.nParams;
            end            
        otherwise
            error('Unknown SPEED option')        
    end
end

model.kern.nParams = numParams;

% To include independent kernel if model.includeInd
model.includeInd = options.includeInd;
if options.includeInd
    for i=1:model.nout,
        model.kernInd.comp{i} = kernCreate(model.indX{i}, 'rbf');
        numParams = numParams + model.kernInd.comp{i}.nParams;
    end
    model.kernInd.nParams = numParams - model.kern.nParams;
end

% Count number of parameters
model.nParams = numParams;

% Change options of the kernel
%fhandle = [model.kernType 'MultimodelKernOptions'];
%if exist(fhandle, 'file')
%    fhandle = str2func(fhandle);
%    model = fhandle(model, options);
%end


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

% Create noise model
switch model.approx
    case {'fitc','pitc', 'dtcvar'}
        % In this structure is easier to put noise in the latent functions
        % at the top level
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
        end
end

% Initializes the q(S) distribution in case of variational learning of
% sensitivities
if model.varS
    % Parameters of the prior
    model.gammas = ones(1, model.nout*model.nlf);
    model.gammasTransform = optimiDefaultConstraint('positive');
    model.nParams = model.nParams + model.nout*model.nlf;
    % Initialization of the qS distribution
    model = spmultimodelVarSInit(model);
end


% Tie options according to the particular kernel employed
fhandle = [model.kernType 'MultimodelTieParam'];
if exist(fhandle, 'file')
    fhandle = str2func(fhandle);
    tieInd = fhandle(model);
else
    warning('spmultimodelCreate:notyingOptions',...
        'Function for tying parameters for this kernel type not implemented yet')
end

if model.isSpeed ~= 1
    % Fix parameters options according to the particular kernel employed
    fhandle = [model.kernType 'MultimodelFixParam'];
    if exist(fhandle, 'file')
        fhandle = str2func(fhandle);
        model = fhandle(model, options);
    else
        tieInd{1} = 1;
        warning('spmultimodelCreate:nofixingOptions',...
            'Function for fixing parameters for this kernel type not implemented yet')
    end
end

if model.isSpeed ~= 3
    model = spmultimodelTieParam(model, tieInd);
else
    model = modelTieParam(model, tieInd);
end

params = modelExtractParam(model);
model = modelExpandParam(model, params);

