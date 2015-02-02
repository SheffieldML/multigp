function g = simwhiteMeanGradient(meanFunction, varargin)

% SIMWHITEMEANGRADIENT Gradient of the parameters of the mean function in
% the multigp model with SIM-WHITE kernel.
% FORMAT
% DESC gives the gradient of the objective function for the parameters of
% the mean function in the multigp model with SIM-WHITE kernel (second order
% differential equation with white noise process input).
% ARG meanFunction : mean function structure to optimise.
% ARG P1, P2, P3 ... : optional additional arguments.
% RETURN g : the gradient of the error function to be minimised.
% 
% SEEALSO : simwhiteMeanCreate, simwhiteMeanOut
%
% COPYRIGHT : Mauricio Alvarez and Neil D. Lawrence, 2008
%
% MODIFICATIONS : David Luengo, 2009

% MULTIGP

gmu = varargin{1}';
gB = gmu./meanFunction.decay;
gD = -gmu.*meanFunction.basal./(meanFunction.decay.*meanFunction.decay);
g = [gB' gD'];
