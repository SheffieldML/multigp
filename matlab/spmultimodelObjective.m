function f = spmultimodelObjective(params, model)

% SPMULTIMODELOBJECTIVE Wrapper function for MODELOPTIMISE objective.

% MULTIGP

model = modelExpandParam(model, params);
f = - modelLogLikelihood(model);
