function g = simglobalKernGradient(kern, outX, latX, dLdKyy, dLdKuy, dLdKuu)

% SIMGLOBALKERNGRADIENT
%
% COPYRIGTH : Mauricio A. Alvarez, 2010

% MULTIGP

kernLat = kern.template.latent;
kernOut = kern.template.output;

gradInverseWidthVector = zeros(1, kern.nlf);
gradDecayVector = zeros(1, kern.nout);

% Compute derivatives of parameters appearing in Kuu
for k = 1:kern.nlf
    % First we need to expand the parameters in the vector to the local
    % kernel
    kernLat.inverseWidth = kern.inverseWidthVector(k);
    % gradLocal(1) is inverse width in rbfKernGradient
    gradLocal = real(kern.funcNames.gradientLat(kernLat, latX{k}, dLdKuu{k}));
    % Before adding the value of the gradient, we must perform a 
    % transformation of the parameters. We implicitly assume this is the 
    % correct transformation, but must be careful for other choices of kernels
    gradInverseWidthVector(k) = gradInverseWidthVector(k) ...
        + gradLocal(1)*kernLat.inverseWidth; 
end
for i = 1:kern.nout,
    % Expand the parameter decay
    kernOut.decay = kern.decayVector(i);
    for j = 1: kern.nlf
        % Expand the parameter inverseWidth
        kernOut.inverseWidth = kern.inverseWidthVector(j);
        kernLat.inverseWidth = kern.inverseWidthVector(j);
        % Gradient of Kyy
        % gradLocal(1) is decay and gradLocal(2) is inversewidth according to simKernGradient        
        gradLocal = real(kern.funcNames.gradientOut(kernOut, outX{i}, dLdKyy{i,j}));
        gradDecayVector(i) = gradDecayVector(i) + gradLocal(1)*kernOut.decay; 
        gradInverseWidthVector(j) = gradInverseWidthVector(j) + gradLocal(2)*kernLat.inverseWidth;
        % Gradient of Kyu
        % gradLocal1(1) is decay and gradLocal1(2) is inversewidth according
        % to simXrbfKernGradient
        gradLocal1 =  real(kern.funcNames.gradientCross(kernOut, kernLat, ...
            outX{i}, latX{j}, (dLdKuy{j,i})'));
        gradDecayVector(i) = gradDecayVector(i) + gradLocal1(1)*kernOut.decay;
        gradInverseWidthVector(j) = gradInverseWidthVector(j) + gradLocal1(2)*kernLat.inverseWidth;
    end
end

g = [gradInverseWidthVector gradDecayVector];