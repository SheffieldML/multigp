function model = simMultimodelFixParam(model, options)

% SIMMULTIMODELFIXPARAM Fix parameters for a sparse multi model with SIM
% FORMAT
% DESC Fix the parameters for a sparse multi model that uses SIM kernel.
% RETURN model : model with fixed parameters included
% ARG model    : model before fixing the parameters
% ARG options  : options for fixing parameters
%
% COPYRIGHT : Mauricio A. Alvarez, 2009

% MULTIGP

if isfield(options, 'typeLf') && ~isempty(options.typeLf)
    
else
    [params, names] = modelExtractParam(model);
    index = findExpression(names, '.* rbf .* variance');
    count = 0;
    for k=1:length(index);
        count = count + 1;
        model.fix(count).index = index(k);
        model.fix(count).value = expTransform(1, 'xtoa');
    end
end

function ind = findExpression(names, pattern)
ind = [];
for i = 1:length(names)
    if(regexp(names{i}, pattern))
        ind = [ind i];
    end
end