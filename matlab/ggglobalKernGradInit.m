function kern = ggglobalKernGradInit(kern)

% GGGLOBALKERNGRADINIT
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% MULTIGP

kern.grad.precisionU = zeros(size(kern.precisionU));
kern.grad.precisionG = zeros(size(kern.precisionG));
kern.grad.sensitivity = zeros(size(kern.sensitivity));
