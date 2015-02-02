function kern = lfmglobalKernExpandParam(kern, params)

% LFMGLOBALKERNEXPANDPARAM
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% MULTIGP

if kern.isMassFixed
    kern.springVector = params(1:kern.nout);
    kern.damperVector = params(kern.nout+1:2*kern.nout);
    kern.inverseWidthVector = params((2*kern.nout)+1:(2*kern.nout)+kern.nlf);
    kern.sensitivity = reshape(params(kern.nlf+2*kern.nout+1:end), kern.nout, kern.nlf);
    massVector = kern.massFixedVal*ones(1, kern.nout);
else
    kern.massVector = params(1:kern.nout);
    kern.springVector = params(kern.nout+1:2*kern.nout);
    kern.damperVector = params((2*kern.nout)+1:3*kern.nout);
    kern.inverseWidthVector = params((3*kern.nout)+1:(3*kern.nout)+kern.nlf);
    kern.sensitivity = reshape(params(kern.nlf+3*kern.nout+1:end), kern.nout, kern.nlf);
    massVector = kern.massVector;
end

kern.alphaVector = kern.damperVector./(2*massVector);
kern.omegaVector = sqrt(kern.springVector./massVector-kern.alphaVector.^2);
kern.gammaVector = kern.alphaVector + 1i*kern.omegaVector;

kern.zetaVector = kern.damperVector./(2*sqrt(massVector.*kern.springVector));
kern.omega_0Vector = sqrt(kern.springVector./massVector);
