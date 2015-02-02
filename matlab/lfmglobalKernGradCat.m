function  g = lfmglobalKernGradCat(kern) 

% LFMGLOBALKERNGRADCAT 
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% MULTIGP

if kern.isMassFixed
    g = [kern.grad.springVector kern.grad.damperVector ...
        kern.grad.inverseWidthVector];
else
    g = [kern.grad.massVector kern.grad.springVector kern.grad.damperVector ...
        kern.grad.inverseWidthVector];
end
g = [g  kern.grad.sensitivity(:)'];

