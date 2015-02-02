function g = spmultimodelGradient(params, model)

% SPMULTIMODELGRADIENT Gradient wrapper for a SPMULTIMODEL model.

% MULTIGP

model = modelExpandParam(model, params);
g = - modelLogLikeGradients(model);

