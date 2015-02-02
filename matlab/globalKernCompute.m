function [Kyy, Kyu, Kuu, KyyInd] = globalKernCompute(kern, outX, latX, gamma)

% GLOBALKERNCOMPUTE
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
        if ~iscell(outX{1}) || (iscell(outX{1}) && length(outX) == 1)
            if (iscell(outX{1}) && length(outX) == 1)
               outX = outX{1}; 
            end
            dimX = zeros(1, kern.nout);
            for i=1:kern.nout
                dimX(i) = size(outX{i},1);
            end
            kernOut2 = kern.template.output;
            Kyy = zeros(sum(dimX), sum(dimX));
            for i=1:kern.nlf
                startOne = 1;
                endOne = 0;
                for j=1:kern.nout
                    endOne = endOne + dimX(j);
                    kernOut = kern.funcNames.transferParamOut(kern, kernOut, j, i);
                    localK = real(kern.funcNames.computeOut(kernOut, outX{j}));
                    Kyy(startOne:endOne, startOne:endOne) = Kyy(startOne:endOne, startOne:endOne) ...
                        + localK;
                    startTwo = 1;
                    endTwo = 0;
                    for k=1:j-1
                        endTwo = endTwo + dimX(k);
                        kernOut2 = kern.funcNames.transferParamOut(kern, kernOut2, k, i);
                        localK = real(kern.funcNames.computeOutCrossOut(kernOut, kernOut2, outX{j}, outX{k}));
                        Kyy(startOne:endOne, startTwo:endTwo) = Kyy(startOne:endOne, startTwo:endTwo) ...
                            + localK;
                        Kyy(startTwo:endTwo, startOne:endOne) = (Kyy(startOne:endOne, startTwo:endTwo)');
                        startTwo = endTwo + 1;
                    end
                    startOne = endOne + 1;
                end
            end
        else
            outX1 = outX{1};
            outX2 = outX{2};
            dimX1 = zeros(1, kern.nout);
            dimX2 = zeros(1, kern.nout);
            for i=1:kern.nout
                dimX1(i) = size(outX1{i},1);
                dimX2(i) = size(outX2{i},1);
            end
            kernOut2 = kern.template.output;
            Kyy = zeros(sum(dimX1), sum(dimX2));
            for i=1:kern.nlf
                startOne = 1;
                endOne = 0;
                startThree = 1;
                endThree = 0;
                for j=1:kern.nout
                    endOne = endOne + dimX1(j);
                    endThree = endThree + dimX2(j);
                    kernOut = kern.funcNames.transferParamOut(kern, kernOut, j, i);                    
                    localK = real(kern.funcNames.computeOut(kernOut, outX1{j}, outX2{j}));
                    Kyy(startOne:endOne, startThree:endThree) = Kyy(startOne:endOne, startThree:endThree) ...
                        + localK;
                    startTwo = 1;
                    endTwo = 0;
                    startFour = 1;
                    endFour = 0;
                    for k=1:j-1
                        endTwo = endTwo + dimX1(k);
                        endFour = endFour + dimX2(k);                        
                        kernOut2 = kern.funcNames.transferParamOut(kern, kernOut2, k, i);
                        localK = real(kern.funcNames.computeOutCrossOut(kernOut, kernOut2, outX1{j}, outX2{k}));
                        Kyy(startOne:endOne, startFour:endFour) = Kyy(startOne:endOne, startFour:endFour) ...
                            + localK;
                        localK = real(kern.funcNames.computeOutCrossOut(kernOut2, kernOut, outX1{k}, outX2{j}));
                        Kyy(startTwo:endTwo, startThree:endThree) = Kyy(startTwo:endTwo, startThree:endThree) ...
                            + localK;
%                         Kyy(startTwo:endTwo, startThree:endThree) = Kyy(startTwo:endTwo, startThree:endThree) ...
%                             + localK';
                        startTwo = endTwo + 1;
                        startFour = endFour + 1;
                    end
                    startOne = endOne + 1;
                    startThree = endThree + 1;
                end
            end
        end
    case {'dtc', 'fitc', 'pitc', 'dtcvar'}       
        Kuu = cell(kern.nlf,1);
        Kyu = cell(kern.nout, kern.nlf);
        if ~strcmp(kern.approx, 'dtc')
            Kyy = cell(kern.nout,1);
            KyyInd = cell(kern.nout, kern.nlf);
        end
        kernLat = kern.template.latent;
        % Compute Kuu
        if ~iscell(latX{1}) || (iscell(latX{1}) && length(latX) == 1)
            for i = 1:kern.nlf
                % First we need to expand the parameters in the vector to the local
                % kernel
                kernLat = kern.funcNames.transferParamLat(kern, kernLat, i);
                Kuu{i} = real(kern.funcNames.computeLat(kernLat, latX{i}));
                if ~isempty(gamma)
                    Kuu{i} = Kuu{i} + gamma(i)*eye(size(Kuu{i}));
                end
            end
        else
            latX1 = latX{1};
            latX2 = latX{2};
            for i = 1:kern.nlf
                % First we need to expand the parameters in the vector to the local
                % kernel
                kernLat = kern.funcNames.transferParamLat(kern, kernLat, i);
                Kuu{i} = real(kern.funcNames.computeLat(kernLat, latX1{i}, latX2{i}));                
            end
            latX = latX{2};
        end        
        for i = 1:kern.nout,
            switch kern.approx
                case {'fitc', 'dtcvar'}
                    Kyy{i} = zeros(size(outX{i},1),1);
                case 'pitc'
                    Kyy{i} = zeros(size(outX{i},1));
            end
            for j = 1: kern.nlf
                kernLat = kern.funcNames.transferParamLat(kern, kernLat, j);
                kernOut = kern.funcNames.transferParamOut(kern, kernOut, i, j);
                switch kern.approx
                    case {'fitc', 'dtcvar'}
                        KyyInd{i,j} =  real(kern.funcNames.computeDiagOut(kernOut, outX{i}));
                        Kyy{i} = Kyy{i} + KyyInd{i,j};
                    case 'pitc'
                        KyyInd{i,j} = real(kern.funcNames.computeOut(kernOut, outX{i}));
                        Kyy{i} = Kyy{i} + KyyInd{i,j};
                end
                Kyu{i,j} = real(kern.funcNames.computeOutCrossLat(kernOut, kernLat, outX{i}, latX{j}));
            end
        end
        
        if strcmp(kern.approx, 'dtc')
            Kyy = [];            
        end        
end
