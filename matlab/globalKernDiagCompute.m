function [kyy, kuu, kyyInd] = globalKernDiagCompute(kern, outX, latX, gamma)

% GLOBALKERNDIAGCOMPUTE
%
% COPYRIGHT : Mauricio A. Alvarez

% MULTIGP

if nargin < 4
    gamma = [];
    if nargin <3
        latX = [];
    end
end


kernOut = kern.template.output;

switch kern.approx
    case 'ftc'
        dimX = zeros(1, kern.nout);
        for i=1:kern.nout
            dimX(i) = size(outX{i},1);
        end        
        kyy = zeros(sum(dimX), 1);
        for i=1:kern.nlf
            startOne = 1;
            endOne = 0;
            for j=1:kern.nout
                endOne = endOne + dimX(j);
                kernOut = kern.funcNames.transferParamOut(kern, kernOut, j, i);
                localK = kern.funcNames.computeDiagOut(kernOut, outX{j});
                kyy(startOne:endOne) = kyy(startOne:endOne) + localK;                
                startOne = endOne + 1;
            end
        end        
    case {'dtc', 'fitc', 'pitc', 'dtcvar'}        
        kuu = cell(kern.nlf,1);        
        if ~strcmp(kern.approx, 'dtc')
            kyy = cell(kern.nout,1);
            kyyInd = cell(kern.nout, kern.nlf);
        end
        kernLat = kern.template.latent;
        if nargout > 1
            % Compute Kuu
            for i = 1:kern.nlf
                % First we need to expand the parameters in the vector to the local
                % kernel
                kernLat = kern.funcNames.transferParamLat(kern, kernLat, i);
                kuu{i} = real(kern.funcNames.computeDiagLat(kernLat, latX{i}));
                if ~isempty(gamma)
                    kuu{i} = kuu{i} + gamma(i);
                end
            end
        end
        for i = 1:kern.nout,
            kyy{i} = zeros(size(outX{i},1),1);
            for j = 1: kern.nlf
                kernLat = kern.funcNames.transferParamLat(kern, kernLat, j);
                kernOut = kern.funcNames.transferParamOut(kern, kernOut, i, j);
                kyyInd{i,j} =  kern.funcNames.computeDiagOut(kernOut, outX{i});
                kyy{i} = kyy{i} + kyyInd{i,j};
            end
        end
end
