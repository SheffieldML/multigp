function kern = lfmglobalKernParamInit(kern)

% LFMGLOBALKERNPARAMINIT
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

kern.springVector = ones(1, kern.nout);
kern.damperVector = ones(1, kern.nout);
kern.inverseWidthVector = ones(1, kern.nlf);
kern.sensitivity = ones(kern.nout, kern.nlf);        

if isfield(kern, 'options') && isfield(kern.options, 'isMassFixed')    
    if isfield(kern.options.isMassFixed, 'state') && kern.options.isMassFixed.state
        kern.isMassFixed = true;
        if isfield(kern.options.isMassFixed, 'massFixedValue')            
            kern.massFixedVal = kern.options.isMassFixed.massFixedValue;
        else
            % The Default value assumed for the masses would be one
            kern.massFixedVal = 1;            
        end
        nTransform = kern.nlf + 2*kern.nout;            
    else
        kern.isMassFixed = false;
        kern.massVector = ones(1, kern.nout);
        nTransform = kern.nlf + 3*kern.nout;
    end
else
    kern.isMassFixed = false;
    kern.massVector = ones(1, kern.nout);
    nTransform = kern.nlf + 3*kern.nout;
end
kern.nParams = nTransform + kern.nout*kern.nlf;
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
