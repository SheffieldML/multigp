function kern = ggglobalKernParamInit(kern)

% GGGLOBALKERNPARAMINIT
%
% COPYRIGHT : Mauricio A. Alvarez, 2010 

% MULTIGP

if isfield(kern, 'options') && isfield(kern.options, 'nout')
    kern.nout = kern.options.nout;    
else
    error('Number of outputs is required for this kernel')
end

if isfield(kern, 'options') && isfield(kern.options, 'nlf')
    kern.nlf = kern.options.nlf;    
else
    error('Number of latent forces is required for this kernel')
end

if isfield(kern, 'options') && isfield(kern.options, 'approx')
    kern.approx = kern.options.approx;    
else
    error('Approximation method is required for this kernel')
end

if ~isfield(kern, 'inputDimension')
    warning('ggglobalKernParamInit:noInputDimension', 'Input dimension has not been provided. Assuming is one.')
    kern.inputDimension = 1;
end

if isfield(kern, 'options') && isfield(kern.options, 'isArd') && kern.options.isArd
    kern.isArd = true;
    if isfield(kern.options, 'tieOutputParams') && kern.options.tieOutputsParams
        kern.precisionU = ones(kern.inputDimension, kern.nlf);
        kern.precisionG = ones(kern.inputDimension, kern.nout);
        kern.nParams = kern.inputDimension*(kern.nlf + kern.nout);
        
    else
        kern.precisionU = ones(kern.inputDimension, kern.nlf);
        kern.precisionG = ones(kern.inputDimension, kern.nout, kern.nlf);
        kern.nParams = kern.inputDimension*(kern.nlf + kern.nout*kern.nlf);                   
    end
    kern.lfParamsTemplate = kern.inputDimension;
    kern.outParamsTemplate = kern.inputDimension;
else
    kern.isArd = false;
    if isfield(kern, 'options') && isfield(kern.options, 'tieOutputParams')
        kern.tieOutputParams = kern.options.tieOutputParams;
        if kern.options.tieOutputParams
            kern.precisionU = ones(1, kern.nlf);
            kern.precisionG = ones(1, kern.nout);
            kern.nParams = kern.nlf + kern.nout;            
        else
            kern.precisionU = ones(1, kern.nlf);
            kern.precisionG = ones(1, kern.nout, kern.nlf);
            kern.nParams = kern.nlf + kern.nout*kern.nlf;
        end        
    else
        kern.tieOutputParams = true;
%         kern.precisionU = rand(1, kern.nlf);
%         kern.precisionG = rand(1, kern.nout);
        kern.precisionU = ones(1, kern.nlf);
        kern.precisionG = ones(1, kern.nout);
        kern.nParams = kern.nlf + kern.nout;        
    end
    kern.lfParamsTemplate = 1;
    kern.outParamsTemplate = 1;
end
kern.sensitivity  = ones(kern.nout, kern.nlf);
%kern.sensitivity  = rand(kern.nout, kern.nlf);

kern.diffParams = 3;

kern.nParams = kern.nParams + kern.nout*kern.nlf;
kern.transforms.index = 1:(kern.nParams-(kern.nout*kern.nlf));
kern.transforms.type = optimiDefaultConstraint('positive');
kern.isStationary = true;
