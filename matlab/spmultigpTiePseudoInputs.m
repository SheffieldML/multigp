function tieInd = spmultigpTiePseudoInputs(model)

% SPMULTIGPTIEPSEUDOINPUTS Gives the indeces to tie the pseudo inputs
% DESC Gives the indexes to tie the pseudo points
% ARG model : input sparse model.
% ARG tieInd : a cell where each component refers to the elements of the
% pseudo inputs to be tied.
%
% COPYRIGHT : Mauricio A. Alvarez,  2010

% MULTIGP

[params, names] = modelExtractParam(model);
indexes = zeros(model.k(1)*model.q, model.nlf);

for k=1:model.nlf
    indexes(:,k) = findExpression(names, ['X_u .* Force ' num2str(k) '\.']);   
end

tieInd = mat2cell(indexes, ones(1, model.k(1)*model.q), model.nlf)';

function ind = findExpression(names, pattern)
ind = [];
for i = 1:length(names)
    if(regexp(names{i}, pattern))
        ind = [ind i];
    end
end