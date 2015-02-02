function   kern = lfmKernGradTransfer(kern, kernOut, localGrad, whichOut, whichLat)

% LFMKERNGRADTRANSFER
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% MULTIGP

params = kern.funcNames.extractOut(kernOut);

if ~kern.isMassFixed
    gradMass = localGrad(1)*params(1);
    kern.grad.massVector(whichOut) = kern.grad.massVector(whichOut) + gradMass;
end

gradSpring = localGrad(2)*params(2);
gradDamper = localGrad(3)*params(3);
gradInvWidth = localGrad(4)*params(4);
kern.grad.springVector(whichOut) = kern.grad.springVector(whichOut) + gradSpring;
kern.grad.damperVector(whichOut) = kern.grad.damperVector(whichOut) + gradDamper;
kern.grad.inverseWidthVector(whichLat) = kern.grad.inverseWidthVector(whichLat) + gradInvWidth;

gradSensitivity = localGrad(5);
kern.grad.sensitivity(whichOut, whichLat) = kern.grad.sensitivity(whichOut, whichLat) + gradSensitivity;

    







