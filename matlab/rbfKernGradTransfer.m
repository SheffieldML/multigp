function   kern = rbfKernGradTransfer(kern, kernLat, localGrad, whichLat)

% RBFKERNGRADTRANSFER
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% MULTIGP

params = kern.funcNames.extractLat(kernLat);

gradInvWith = localGrad(1)*params(1);

kern.grad.inverseWidthVector(whichLat) = kern.grad.inverseWidthVector(whichLat) ...
    + gradInvWith';





