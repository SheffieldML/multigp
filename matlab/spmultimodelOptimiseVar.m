function model = spmultimodelOptimiseVar(model, display, iters, niters)

% SPMULTIMODELOPTIMISEVAR Variational learning of sensitivities
% FORMAT
% DESC
% Function that makes the optimization of the hyperparameters of the
% spmultimodel and the variational updates.
% ARG model : the model to be optimised
% ARG display : flag that indicates if the optimization of the bound should
% be shown
% ARG iters  : number of iterations for optimization cycle
% ARG niters : number of iterations for optimization of hyperparameters
% RETURN model : the optimised model
%
% COPYRIGHT : Mauricio Alvarez, 2010
%
% MULTIGP

if nargin<4
    niters = 50;
end

options = optOptions;

j = 1;
finit = modelLogLikelihood(model);
fprintf('Initial value of the BOUND: %f\n', finit)
while j<= iters
    fprintf('\n')
    fprintf('ITERATION NUMBER: %d \n', j)
    % Optimize hyperparameters
    fprintf('OPTIMIZING HYPERPARAMETERS. \n')
    model = modelOptimise(model, [], [], display, niters);
    % Evaluate the moments
    params = modelExtractParam(model);
    model = modelExpandParam(model, params);
    fopt = modelLogLikelihood(model);
    fprintf('Bound AFTER optimizing the hyperparemeters: %f\n', fopt)
    % Update q(u) and q(S) and use those parameters to recompute the
    % likelihood
    fprintf('UPDATING DISTRIBUTIONS. \n')
    if model.isSpeed == 1 && isfield(model, 'subSpeed') && model.subSpeed == 2
        model = spmultimodelUpdateVariational2(model);
    else
        model = spmultimodelUpdateVariational(model);
    end
    params = modelExtractParam(model);
    model = modelExpandParam(model, params);
    fdis = modelLogLikelihood(model);
    fprintf('Bound AFTER updating distributions and optimizing: %f\n', fdis)
    fprintf('Bound INCREASE: %f\n', fdis - fopt)
    if abs(fdis - fopt) < options(3)
        break;
    end
    j = j + 1;
end




