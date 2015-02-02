function [Kyy, Kyu, Kuu] = simglobalKernComputeTest(kern, latX, outXs, latXs, gamma)

% SIMGLOBALKERNCOMPUTETEST
%
% COPYRIGTH : Mauricio A. Alvarez, 2010

% MULTIGP

if nargin < 5
    gamma = [];
end

Kuu = cell(kern.nlf,1);
Kyu = cell(kern.nout, kern.nlf);
Kyy = cell(kern.nout, kern.nlf);


kernLat = kern.template.latent;
kernOut = kern.template.output;

% Compute Kuu
for k = 1:kern.nlf
    % First we need to expand the parameters in the vector to the local
    % kernel
    kernLat.inverseWidth = kern.inverseWidthVector(k);
    Kuu{k} = real(kern.funcNames.computeLat(kernLat, latXs, latX{k}));
    if ~isempty(gamma)
        Kuu{k} = Kuu{k} + gamma(k)*eye(size(Kuu{k}));
    end
end
for i = 1:kern.nout,
    % Expand the parameter decay
    kernOut.decay = kern.decayVector(i);
    for j = 1: kern.nlf
        % Expand the parameter inverseWidth
        kernOut.inverseWidth = kern.inverseWidthVector(j);
        kernLat.inverseWidth = kern.inverseWidthVector(j);        
        % Compute Kff        
        Kyy{i,j} = real(kern.funcNames.computeOut(kernOut, outXs{i}));
        % Compute Kfu, which corresponds to K_{\hat{f}}u, really.
        Kyu{i,j} = real(kern.funcNames.computeCross(kernOut, kernLat, outXs{i}, latX{j}));
    end
end
