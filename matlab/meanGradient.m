function g = meanGradient(meanFunction, varargin)

% MEANGRADIENT Gradient of the parameters of the mean function in the
% multigp model
% FORMAT
% DESC gives the gradient of the objective function for the parameters of
% the mean function in the multigp model 
% ARG meanFunction : mean function structure to optimise.
% ARG P1, P2, P3 ... : optional additional arguments.
% RETURN g : the gradient of the error function to be minimised.
% 
% SEEALSO : meanCreate, meanCompute
%
% COPYRIGHT : Mauricio Alvarez and Neil D. Lawrence, 2008

% MULTIGP

fhandle = str2func([meanFunction.type 'MeanGradient']);
g = fhandle( meanFunction, varargin{:});

factors = meanFactors(meanFunction, 'gradfact');
g(factors.index) = g(factors.index).*factors.val;

