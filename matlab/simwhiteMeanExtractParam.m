function [params, names] = simwhiteMeanExtractParam(meanFunction)

% SIMWHITEMEANEXTRACTPARAM Extract parameters from the SIM-WHITE mean function structure.
% FORMAT
% DESC Extract parameters from the mean funtion structure of the SIM-WHITE
% model into a vector of parameters for optimisation.
% ARG meanFunction : the mean function structure containing the parameters
% to be extracted.
% RETURN param : vector of parameters extracted from the kernel. 
%
% DESC Extract parameters and their names from mean funtion structure of
% the SIM-WHITE model
% ARG meanFunction : the mean function structure containing the parameters
% to be extracted.
% RETURN param : vector of parameters extracted from the kernel. 
% RETURN names : cell array of strings containing parameter names.
%
% SEEALSO  simwhiteMeanCreate, simwhiteMeanExpandParam, simwhiteKernCreate,
% simwhiteKernExtractParam
%
% COPYRIGHT : Mauricio Alvarez and Neil D. Lawrence, 2008
%
% MODIFICATIONS : David Luengo, 2009

% MULTIGP

params = [meanFunction.basal' meanFunction.decay'];
if nargout > 1
    names = cell(1, 2*meanFunction.nParams/2);
    for i=1:meanFunction.nParams/2        
        names{i} = ['simwhite ' num2str(i) ' basal'];
    end    
    for i=meanFunction.nParams/2+1:2*meanFunction.nParams/2        
        names{i} = ['simwhite ' num2str(i-meanFunction.nParams/2) ' decay'];
    end    
end
