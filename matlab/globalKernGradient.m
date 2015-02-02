function g = globalKernGradient(kern, outX, dLdKyy, latX, dLdKuy, dLdKuu)

% GLOBALKERNGRADIENT
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% MULTIGP


dimX = zeros(1, kern.nout);
for i=1:kern.nout
    dimX(i) = size(outX{i},1);
end

kernOut = kern.template.output;
kern = kern.funcNames.gradInit(kern);

switch kern.approx
    case 'ftc'
        kernOut2 = kern.template.output;
        for i=1:kern.nlf
            startOne = 1;
            endOne = 0;
            for j=1:kern.nout
                endOne = endOne + dimX(j);
                kernOut = kern.funcNames.transferParamOut(kern, kernOut, j, i);
                g1 = kern.funcNames.gradientOut(kernOut, outX{j}, dLdKyy(startOne:endOne, startOne:endOne));
                kern = kern.funcNames.transferGradOut(kern, kernOut, g1, j, i);
                startTwo = 1;
                endTwo = 0;
                for k=1:j-1
                    endTwo = endTwo + dimX(k);
                    kernOut2 = kern.funcNames.transferParamOut(kern, kernOut2, k, i);
                    [g1, g2] = kern.funcNames.gradientOutCrossOut(kernOut, kernOut2, outX{j}, ...
                        outX{k}, dLdKyy(startOne:endOne, startTwo:endTwo));
                    kern = kern.funcNames.transferGradOut(kern, kernOut, 2*g1, j, i);
                    kern = kern.funcNames.transferGradOut(kern, kernOut2, 2*g2, k, i);
                    startTwo = endTwo + 1;
                end
                startOne = endOne + 1;
            end
        end
    case {'dtc', 'fitc', 'pitc', 'dtcvar'}
        kernLat = kern.template.latent;
        % Compute derivatives of parameters appearing in Kuu
        for i = 1:kern.nlf
            kernLat = kern.funcNames.transferParamLat(kern, kernLat, i);
            g1 = kern.funcNames.gradientLat(kernLat, latX{i}, dLdKuu{i});
            kern = kern.funcNames.transferGradLat(kern, kernLat, g1, i);
        end
        % Here we have two options, according to the form for dLdKyy. It
        % could have been derivative per force and output or the sum over
        % forces.     
        if strcmp(kern.approx, 'dtc')
            for i = 1:kern.nout,
                for j = 1: kern.nlf
                    kernOut = kern.funcNames.transferParamOut(kern, kernOut, i, j);
                    kernLat = kern.funcNames.transferParamLat(kern, kernLat, j);
                    [g1, g2] = kern.funcNames.gradientOutCrossLat(kernOut, kernLat, outX{i}, ...
                        latX{j}, (dLdKuy{j,i})');
                    kern = kern.funcNames.transferGradOut(kern, kernOut, g1, i, j);
                    kern = kern.funcNames.transferGradLat(kern, kernLat, g2, j);
                end
            end            
        else
            colDer = size(dLdKyy,2);
            if colDer == 1
                for i = 1:kern.nout,
                    for j = 1: kern.nlf
                        kernOut = kern.funcNames.transferParamOut(kern, kernOut, i, j);
                        g1 = kern.funcNames.gradientOut(kernOut, outX{i}, dLdKyy{i});
                        kern = kern.funcNames.transferGradOut(kern, kernOut, g1, i, j);
                        kernLat = kern.funcNames.transferParamLat(kern, kernLat, j);
                        [g1, g2] = kern.funcNames.gradientOutCrossLat(kernOut, kernLat, outX{i}, ...
                            latX{j}, (dLdKuy{j,i})');
                        kern = kern.funcNames.transferGradOut(kern, kernOut, g1, i, j);
                        kern = kern.funcNames.transferGradLat(kern, kernLat, g2, j);
                    end
                end
            else                
                for i = 1:kern.nout,
                    for j = 1: kern.nlf
                        kernOut = kern.funcNames.transferParamOut(kern, kernOut, i, j);
                        g1 = kern.funcNames.gradientOut(kernOut, outX{i}, dLdKyy{i,j});
                        kern = kern.funcNames.transferGradOut(kern, kernOut, g1, i, j);
                        kernLat = kern.funcNames.transferParamLat(kern, kernLat, j);
                        [g1, g2] = kern.funcNames.gradientOutCrossLat(kernOut, kernLat, outX{i}, ...
                            latX{j}, (dLdKuy{j,i})');
                        kern = kern.funcNames.transferGradOut(kern, kernOut, g1, i, j);
                        kern = kern.funcNames.transferGradLat(kern, kernLat, g2, j);
                    end
                end
            end
        end
end
g = kern.funcNames.gradCat(kern);

