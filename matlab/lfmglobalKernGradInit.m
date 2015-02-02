function kern = lfmglobalKernGradInit(kern)

% LFMGLOBALKERNGRADINIT
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% MULTIGP

if ~kern.isMassFixed
    kern.grad.massVector = zeros(size(kern.massVector));
end
kern.grad.springVector = zeros(size(kern.springVector));
kern.grad.damperVector = zeros(size(kern.damperVector));
kern.grad.inverseWidthVector = zeros(size(kern.inverseWidthVector));
kern.grad.sensitivity = zeros(size(kern.sensitivity));

