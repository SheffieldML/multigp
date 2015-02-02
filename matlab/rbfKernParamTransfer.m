function rbfKern = rbfKernParamTransfer(kern, rbfKern, whichLatent)

% RBFKERNPARAMTRANSFER
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% MULTIGP

rbfKern.inverseWidth = kern.inverseWidthVector(whichLatent);

