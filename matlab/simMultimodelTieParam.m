function tieInd  = simMultimodelTieParam(model)

% SIMMULTIMODELTIEPARAM Tie parameters for a sparse multimodel with SIM
% FORMAT
% DESC Tie the parameters for a sparse multimodel that uses SIM kernel.
% RETURN tieInd : cell with elements containing the indexes of parameters
% to tie.
% ARG model : model created
%
% COPYRIGHT : Mauricio A. Alvarez, 2009, 2010

% MULTIGP

[params, names] = modelExtractParam(model);

cont = 0;
% Tie initialization for inverse widths
indexInvWidths = findExpression(names, '.* inverse width');
namesInvWidths = names(indexInvWidths);
tieInd = cell(1, model.nlf);

% Tie initialization for decays 
indexDecaysKernel = findExpression(names, 'component .* decay');
indexDecaysMean = findExpression(names, 'mean Func .* decay');
namesDecaysKernel = names(indexDecaysKernel);

if ~model.varS
    % Since for this we have the structure, we carefully select the
    % corresponding indexes for tying
    for i = 1:model.nlf
        index = findExpression(namesInvWidths, ['component ' num2str(i) ' .*']);
        cont = cont + 1;
        tieInd{cont} = indexInvWidths(index);
    end
    for i=1:model.nout
        whichLatent = model.indLatGivenOut{i};
        indDecayKernel = zeros(1, model.numLatGivenOut(i));
        for j=1:model.numLatGivenOut(i),
            whichOutput = model.indOutGivenLat{whichLatent(j)};
            indDecayKernel(j) = findExpression(namesDecaysKernel,  ['component ' num2str(whichLatent(j))...
                ' .* '  num2str(find(whichOutput == i))  ' decay']);
        end
        index = [indexDecaysKernel(indDecayKernel) indexDecaysMean(i)];
        cont = cont + 1;
        tieInd{cont} = index;
    end
else
    % Exploit the periodicity of the inverse widths
    %index1 = findExpression(namesInvWidths, 'component .* rbf .* inverse width');
    startIndex = 1; % First RBF component
    endIndex = 0;
    int = model.nout + 1;
    for i = 1:model.nlf
        endIndex = endIndex + int;
        cont = cont + 1;
        tieInd{cont} = indexInvWidths(startIndex:endIndex);
        startIndex = endIndex + 1;
    end
       
    % Exploit the periodicity in which the decay terms appear
    
    for i=1:model.nout
        indDecayKernel = i:model.nout:length(indexDecaysKernel);
        if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
            index = [indexDecaysKernel(indDecayKernel) indexDecaysMean(i)];
        else
            index = indexDecaysKernel(indDecayKernel);
        end
        cont = cont + 1;
        tieInd{cont} = index;
    end

end

function ind = findExpression(names, pattern)
ind = [];
for i = 1:length(names)
    if(regexp(names{i}, pattern))
        ind = [ind i];
    end
end