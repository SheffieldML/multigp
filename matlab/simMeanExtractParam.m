function [params, names] = simMeanExtractParam(meanFunction)

% SIMMEANEXTRACTPARAM Extract parameters from the SIM MEAN FUNCTION structure.
% FORMAT
% DESC Extract parameters from the mean funtion structure of the sim model
% into a vector of parameters for optimisation.
% ARG meanFunction : the mean function structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. 
%
% FORMAT
% DESC Extract parameters and their names from mean funtion structure of
% the sim model
% ARG meanFunction : the mean function structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. 
% RETURN names : cell array of strings containing parameter names.
%
% SEEALSO  simMeanCreate, simMeanExpandParam, simKernCreate,
% simkernExtractParam
%
% COPYRIGHT : Mauricio Alvarez and Neil D. Lawrence, 2008

% MULTIGP

params = [meanFunction.basal' meanFunction.decay'];
if nargout > 1
    names = cell(1, 2*meanFunction.nParams/2);
    for i=1:meanFunction.nParams/2        
        names{i} = ['sim ' num2str(i) ' basal'];
    end    
    for i=meanFunction.nParams/2+1:2*meanFunction.nParams/2        
        names{i} = ['sim ' num2str(i-meanFunction.nParams/2) ' decay'];
    end    
end
