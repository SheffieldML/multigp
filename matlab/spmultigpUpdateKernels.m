function model = spmultigpUpdateKernels(model)

% SPMULTIGPUPDATEKERNELS 
% FORMAT
% DESC Update the kernels that are needed for the sparse multigp
% ARG model : the model structure containing the model parameters
% RETURN model : the model structure with updated kernels.
%
% COPYRIGHT : Mauricio A. Alvarez, 2008, 2009

% MULTIGP

switch model.approx
    case {'dtc','fitc', 'pitc','dtcvar'}
        model = sparseKernCompute(model);        
    otherwise        
end

if isfield(model, 'beta') && ~isempty(model.beta)
    model = spmultigpUpdateAD(model);
end