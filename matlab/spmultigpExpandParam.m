function model = spmultigpExpandParam(model, params)

% SPMULTIGPEXPANDPARAM Expand a parameter vector into a SPMULTIGP model.
% FORMAT
% DESC expands the model parameters to a structure containing
% the information about a sparse multi-output Gaussian process.
% ARG model : the sparse model structure containing the information about
% the model.
% ARG params : a vector of parameters from the model.
% RETURN model : the model structure containing the information about
% the sparse model updated with the new parameter vector.
%
% SEEALSO : multigpCreate, spmultigpExtractParam, multigpExtractParam, 
%
%
% COPYRIGHT : Mauricio Alvarez, 2008


% MULTIGP

model.X_u = reshape(params, sum(model.k), model.q);
startVal = 1;
endVal = 0;
for i = 1:model.nlf
    endVal = endVal + model.k(i);
    model.X{i} = model.X_u(startVal:endVal,:);
    startVal = endVal + 1;
end




