function ggKern = ggKernParamTransfer(kern, ggKern, whichOutput, whichLatent)

% GGKERNPARAMTRANSFER
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% MULTIGP

ggKern.precisionU = kern.precisionU(:,whichLatent);
if kern.tieOutputParams
    ggKern.precisionG = kern.precisionG(:,whichOutput);
else
    ggKern.precisionG = kern.precisionG(:,whichOutput, whichLatent);
end
ggKern.sensitivity = kern.sensitivity(whichOutput, whichLatent);

