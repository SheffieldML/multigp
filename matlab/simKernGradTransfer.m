function   kern = simKernGradTransfer(kern, kernOut, localGrad, whichOut, whichLat)

% SIMKERNGRADTRANSFER
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% MULTIGP

params = kern.funcNames.extractOut(kernOut);

gradDecay = localGrad(1)*params(1);
gradInvWidth = localGrad(2)*params(2);

kern.grad.decayVector(whichOut) = kern.grad.decayVector(whichOut) + gradDecay;
kern.grad.inverseWidthVector(whichLat) = kern.grad.inverseWidthVector(whichLat) + gradInvWidth;

if ~kern.isVarS
    if kern.isNegativeS 
        gradSensitivity = localGrad(3);
    else
        gradSensitivity = localGrad(3)*params(3);
    end
    kern.grad.sensitivity(whichOut, whichLat) = kern.grad.sensitivity(whichOut, whichLat) + gradSensitivity;
end
    







