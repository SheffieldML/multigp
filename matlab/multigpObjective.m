function f = multigpObjective(params, model)

% MULTIGPOBJECTIVE Wrapper function for MULTIGPOPTIMISE objective.

% MULTIGP

model = modelExpandParam(model, params);
f = - multigpLogLikelihood(model);
