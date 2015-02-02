function [model, params] = multigpOptimise(model, display, iters)

% MULTIGPOPTIMISE Optimise the inducing variable multigp based kernel.
% FORMAT
% DESC optimises the Gaussian
%	process  model for a given number of iterations.
% RETURN model : the optimised model.
% RETURN params : the optimised parameter vector.
% ARG model : the model to be optimised.
% ARG display : whether or not to display while optimisation proceeds,
%	   set to 2 for the most verbose and 0 for the least verbose.
% ARG iters : number of iterations for the optimisation.
%	
% SEEALSO : scg, conjgrad, multigpOptimiseCreate,
% multigpOptimiseGradient, multigpOptimiseObjective
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006
%
% MODIFIED : Mauricio A. Alvarez, 2008

% MULTIGP

if nargin < 3
  iters = 1000;
  if nargin < 2
    display = 1;
  end
end

params = modelExtractParam(model);

options = optOptions;
if display
    options(1) = 1;
    if length(params) <= 100 && display > 1
        options(9) = 1;
    end
end
options(14) = iters;

if isfield(model, 'optimiser')
  optim = str2func(model.optimiser);
else
  optim = str2func('conjgrad');
end

if strcmp(func2str(optim), 'optimiMinimize')
  % Carl Rasmussen's minimize function 
  params = optim('multigpObjectiveGradient', params, options, model);
else
  % NETLAB style optimization.
  params = optim('multigpObjective', params,  options, ...
                 'multigpGradient', model);
end

%model = multigpExpandParam(model, params);

model = modelExpandParam(model, params);
