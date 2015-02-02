function kern = simglobalKernGradInit(kern)

% SIMGLOBALKERNGRADINIT
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% MULTIGP

kern.grad.inverseWidthVector = zeros(size(kern.inverseWidthVector));
kern.grad.decayVector = zeros(size(kern.decayVector));

if ~kern.isVarS
    kern.grad.sensitivity = zeros(size(kern.sensitivity));
end
