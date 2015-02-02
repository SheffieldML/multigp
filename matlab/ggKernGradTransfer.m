function   kern = ggKernGradTransfer(kern, kernOut, localGrad, whichOut, whichLat)

% GGKERNGRADTRANSFER
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% MULTIGP

params = kern.funcNames.extractOut(kernOut);
gradInvWithLatent = localGrad(1:kern.lfParamsTemplate).*params(1:kern.lfParamsTemplate);
gradInvWidthOutput = localGrad(kern.lfParamsTemplate+1:kern.lfParamsTemplate+kern.outParamsTemplate)...
    .*params(kern.lfParamsTemplate+1:kern.lfParamsTemplate+kern.outParamsTemplate);

gradSensitivity = localGrad(end);

kern.grad.precisionU(:, whichLat) = kern.grad.precisionU(:, whichLat) + gradInvWithLatent';
if kern.tieOutputParams
    kern.grad.precisionG(:, whichOut) = kern.grad.precisionG(:, whichOut) + gradInvWidthOutput;
else
    kern.grad.precisionG(:, whichOut, whichLat) = kern.grad.precisionG(:, whichOut, whichLat) + gradInvWidthOutput;
end

kern.grad.sensitivity(whichOut, whichLat) = kern.grad.sensitivity(whichOut, whichLat) + gradSensitivity;








