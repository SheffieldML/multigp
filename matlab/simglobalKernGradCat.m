function  g = simglobalKernGradCat(kern) 

% SIMGLOBALKERNGRADCAT 
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% MULTIGP

g = [kern.grad.decayVector kern.grad.inverseWidthVector];

if ~kern.isVarS
    g = [g  kern.grad.sensitivity(:)'];
end
