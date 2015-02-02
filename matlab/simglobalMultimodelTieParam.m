function tieInd  = simglobalMultimodelTieParam(model)

% SIMGLOBALMULTIMODELTIEPARAM Tie parameters for a sparse multimodel 
% FORMAT
% DESC Tie the parameters for a sparse multimodel that uses SIM kernel.
% RETURN tieInd : cell with elements containing the indexes of parameters
% to tie.
% ARG model : model created
%
% COPYRIGHT : Mauricio A. Alvarez, 2009, 2010

% MULTIGP

[params, names] = modelExtractParam(model);

% Tie decays
indexDecaysKernel = findExpression(names, 'kernel .* decay');
indexDecaysMean = findExpression(names, 'mean Func .* decay');
indexDecays = [indexDecaysKernel' indexDecaysMean'];
tieInd = (mat2cell(indexDecays, ones(1, model.nout), 2))';

function ind = findExpression(names, pattern)
ind = [];
for i = 1:length(names)
    if(regexp(names{i}, pattern))
        ind = [ind i];
    end
end