function kern = simglobalKernExpandParam(kern, params)

% SIMGLOBALKERNEXPANDPARAM
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% MULTIGP

kern.decayVector = params(1:kern.nout);
kern.inverseWidthVector = params(kern.nout+1:kern.nout+kern.nlf);
if ~kern.isVarS
    kern.sensitivity = reshape(params(kern.nlf+kern.nout+1:end), kern.nout, kern.nlf);
end    
