function [params, names] =  meanExtractParam(meanFunction)

% MEANEXTRACTPARAM Extract parameters from a MEAN FUNCTION structure.
% FORMAT
% DESC Extract parameters from a mean funtion structure 
% into a vector of parameters for optimisation.
% ARG model : the mean function structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the mean function. 
%
% FORMAT
% DESC Extract parameters and their names from a mean funtion structure 
% ARG model : the mean function structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the mean function. 
% RETURN names : cell array of strings containing parameter names.
%
% SEEALSO  meanCreate, meanExpandParam, kernCreate,
% kernExtractParam
%
% COPYRIGHT : Mauricio Alvarez and Neil D. Lawrence, 2008

% MULTIGP

fhandle = str2func([meanFunction.type 'MeanExtractParam']);

if nargout > 1
  [params, names] = fhandle(meanFunction);
else
  params = fhandle(meanFunction);
end

% Check if parameters are being optimised in a transformed space.
if ~isempty(meanFunction.transforms)
  for i = 1:length(meanFunction.transforms)
    index = meanFunction.transforms(i).index;
    fhandle = str2func([meanFunction.transforms(i).type 'Transform']);
    params(index) = fhandle(params(index), 'xtoa');
  end
end

params = params*meanFunction.paramGroups;
