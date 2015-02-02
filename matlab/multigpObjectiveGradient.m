function [f, g] = multigpObjectiveGradient(params, model)

% MULTIGPOBJECTIVEGRADIENT Wrapper for MULTIGP objective and gradient.
% FORMAT
% DESC returns the negative log likelihood of a Gaussian process
% model given the model structure and a vector of parameters. This
% allows the use of NETLAB minimisation functions to find the model
% parameters.
% ARG params : the parameters of the model for which the objective
% will be evaluated.
% ARG model : the model structure for which the objective will be
% evaluated.
% RETURN f : the negative log likelihood of the GP model.
% RETURN g : the gradient of the negative log likelihood of the GP
% model with respect to the parameters.
%
% SEEALSO : minimize, gpCreate, gpGradient, gpLogLikelihood, gpOptimise
% 
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% MULTIGP

% Check how the optimiser has given the parameters
if size(params, 1) > size(params, 2)
  % As a column vector ... transpose everything.
  transpose = true;
  model = modelExpandParam(model, params');
else
  transpose = false;
  model = modelExpandParam(model, params);
end

f = - multigpLogLikelihood(model);
if nargout > 1
  g = - modelLogLikeGradients(model);
end
if transpose
  g = g';
end