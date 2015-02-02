function alpha = multigpComputeAlpha(model, m)

% MULTIGPCOMPUTEALPHA Update the vector `alpha' for computing posterior mean quickly.
% FORMAT
% DESC updates the vectors that are known as `alpha' in the support
% vector machine, in other words invK*y, where y is the target values.
% ARG model : the model for which the alphas are going to be
% updated.
% ARG m : the values of m for which the updates will be made.
% RETURN model : the model with the updated alphas.
%
% SEEALSO : multigpCreate, multigpUpdateAD, multigpUpdateKernels
%
% COPYRIGHT : Neil D. Lawrence, 2008

% MULTIGP

if nargin < 2
  m = model.m;
end

switch model.approx
 case 'ftc'
  alpha = model.invK*m;
 case {'dtc','fitc','pitc', 'dtcvar'}
     if isfield(model, 'beta') && ~isempty(model.beta)
         alpha = model.AinvKuyDinvy;
     else
         alpha =[];
     end
 otherwise
  error('Alpha update not yet implemented for sparse kernels');
end