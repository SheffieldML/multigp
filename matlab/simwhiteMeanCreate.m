function meanFunction = simwhiteMeanCreate(q, d, options)

% SIMWHITEMEANCREATE mean function structure for the SIM-WHITE kernel.
% FORMAT
% DESC creates the mean function for a multi output
% GP model based in the SIM-WHITE kernel (first order differential equation
% with white noise input process). The outputs of the model are generated
% according to
%
%       mean_q = B_q/D_q
%
% where mean_q is an output constant corresponding to the mean of the 
% output function q, B_q is basal transcription and D_q is the decay
% constant.
% RETURN model : the structure for the multigp model
% ARG q : input dimension size.
% ARG d : output dimension size.
% ARG options : contains the options for the MEAN of the MULTIGP model.
%
% SEE ALSO: simwhiteKernParamInit, simwhiteKernCompute
%
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008 
%
% MODIFICATIONS : David Luengo, 2009

% MULTIGP

if q > 1
  error('SIM-WHITE MEAN FUNCTION only valid for one-D input.')
end

meanFunction.type = 'simwhite';
meanFunction.basal  = ones(d,1);
meanFunction.decay = ones(d,1);
meanFunction.transforms.index = d+1:2*d;
meanFunction.transforms.type = optimiDefaultConstraint('positive');
% Only the parameters of basal rates are counted. The springs are already
% counted in the kernel
meanFunction.nParams = 2*d;
