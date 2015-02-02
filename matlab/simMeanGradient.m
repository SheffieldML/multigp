function g = simMeanGradient(meanFunction, varargin)

% SIMMEANGRADIENT Gradient of the parameters of the mean function in the
% multigp model with SIM kernel
% FORMAT
% DESC gives the gradient of the objective function for the parameters of
% the mean function in the multigp model with LFM kernel (second order
% differential equation).
% ARG meanFunction : mean function structure to optimise.
% ARG P1, P2, P3 ... : optional additional arguments.
% RETURN g : the gradient of the error function to be minimised.
% 
% SEEALSO : simMeanCreate, simMeanOut
%
% COPYRIGHT : Mauricio Alvarez and Neil D. Lawrence, 2008

% MULTIGP

gmu = varargin{1}';
gB = gmu./meanFunction.decay;
gD = -gmu.*meanFunction.basal./(meanFunction.decay.*meanFunction.decay);
g = [gB' gD'];
