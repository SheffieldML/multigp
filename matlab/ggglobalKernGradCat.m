function g = ggglobalKernGradCat(kern)

% GGGGLOBALKERNGRADCAT
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% MULTIGP

g = [kern.grad.precisionU(:)' kern.grad.precisionG(:)' kern.grad.sensitivity(:)'];

