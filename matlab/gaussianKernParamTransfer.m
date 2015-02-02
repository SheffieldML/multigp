function gaussianKern = gaussianKernParamTransfer(kern, gaussianKern, whichLatent)

% GUASSIANKERNPARAMTRANSFER
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% MULTIGP

gaussianKern.precisionU = kern.precisionU(:, whichLatent);

