function lfmKern = lfmKernParamTransfer(kern, lfmKern, whichOutput, whichLatent)

% LFMKERNPARAMTRANSFER
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% MULTIGP
if kern.isMassFixed
    lfmKern.mass = kern.massFixedVal;
else
    lfmKern.mass = kern.massVector(whichOutput);
end
lfmKern.spring = kern.springVector(whichOutput);
lfmKern.damper = kern.damperVector(whichOutput);
lfmKern.inverseWidth = kern.inverseWidthVector(whichLatent);
lfmKern.sensitivity = kern.sensitivity(whichOutput, whichLatent);

