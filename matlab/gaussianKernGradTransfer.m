function kern = gaussianKernGradTransfer(kern, kernLat, localGrad, whichLat)

% GUASSIANKERNGRADTRANSFER
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% MULTIGP

params = kern.funcNames.extractLat(kernLat);

gradInvWithLatent = localGrad(1:kern.lfParamsTemplate).*params(1:kern.lfParamsTemplate);

kern.grad.precisionU(:, whichLat) = kern.grad.precisionU(:, whichLat) + gradInvWithLatent';
