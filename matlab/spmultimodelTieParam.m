function model = spmultimodelTieParam(model, paramsList)

% SPMULTIMODELTIEPARAM Tie parameters of a spmultimodel together.
% FORMAT
% DESC groups of parameters of a model to be seen as one parameter
% during optimisation of the model. It is a simplified version of
% modelTieParam.m
% ARG model : the model for which parameters are being tied
% together.
% ARG paramsList : indices of parameteres to group together. The
% indices are provided in a cell array. Each cell in the array
% contains a vector of indices of parameters that should be
% considered as one parameter. Each group of parameters in each
% cell should obviously be mutually exclusive.
% Alternatively each element of the cell array can be a string which
% is interpreted as a regular expression of names of parameters
% (as returned by modelExtractParam) to be tied.
% RETURN model : the model with the parameters grouped together.
%
% SEEALSO : modelTieParam.m
% 
% COPYRIGHT : Neil D. Lawrence, 2003, 2006, 2008
%
% MODIFICATIONS: Antti Honkela, 2009
%
% MODIFICATIONS : Mauricio A. Alvarez, 2010

% MULTIGP

if ~isfield(model, 'paramGroups')
  if isfield(model, 'nParams')
    model.paramGroups = speye(model.nParams);
  elseif isfield(model, 'numParams')
    model.paramGroups = speye(model.numParams);
  else
    error('Model does not list number of parameters.');
  end
end
colToDelete = [];

for i = 1:length(paramsList)
    paramIndices = sort(paramsList{i});
    model.paramGroups(paramIndices, paramIndices(1)) = 1;
    colToDelete = [colToDelete paramIndices(2:end)];
end
  
model.paramGroups(:, colToDelete) = [];
if isfield(model, 'nParams')
  % Update to the new number of parameters.
  model.nParams = size(model.paramGroups, 2);
elseif isfield(model, 'numParams')
  model.numParams = size(model.paramGroups, 2);
end
