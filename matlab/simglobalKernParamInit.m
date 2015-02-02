function kern = simglobalKernParamInit(kern)

% SIMGLOBALKERNPARAMINIT
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

kern.inverseWidthVector = ones(1, kern.nlf);
kern.decayVector = ones(1, kern.nout);

if isfield(kern, 'options') && isfield(kern.options, 'isVarS') ...        
        && kern.options.isVarS,
    kern.isVarS = true;
    kern.nParams = kern.nlf + kern.nout;
    nTransform = kern.nlf + kern.nout;    
else
    kern.isVarS = false;
    kern.sensitivity = ones(kern.nout, kern.nlf);        
    if isfield(kern, 'options') && isfield(kern.options, 'isNegativeS') ...
            && kern.options.isNegativeS,
        kern.isNegativeS = true;
        nTransform = kern.nlf + kern.nout;                
    else
        kern.isNegativeS = false;
        nTransform = kern.nlf + kern.nout + kern.nout*kern.nlf;                        
    end
    kern.nParams = kern.nlf + kern.nout + kern.nout*kern.nlf;    
end

kern.transforms.index = 1:nTransform;
kern.transforms.type = optimiDefaultConstraint('positive');

if isfield(kern, 'options') && isfield(kern.options, 'isStationary') ...
        && kern.options.isStationary,
   kern.isStationary = true;
else
   kern.isStationary = false;
end

if isfield(kern, 'options') && isfield(kern.options, 'isNormalised') ...
        && kern.options.isNormalised,
    kern.isNormalised = true;
else
    kern.isNormalised = false;
end

kern.positiveTime = true;
