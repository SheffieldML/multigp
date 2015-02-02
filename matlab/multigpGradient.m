function g = multigpGradient(params, model)

% MULTIGPGRADIENT Gradient wrapper for a MULTIGP model.

% MULTIGP

model = modelExpandParam(model, params);
g = - modelLogLikeGradients(model);
