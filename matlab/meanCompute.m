function m = meanCompute(meanFunction, X, varargin)

% MEANCOMPUTE Give the output of the lfm mean function model for given X.
% FORMAT
% DESC gives the output of a mean function model for a given input X.
% ARG meanFunction : structure specifying the model.
% ARG X : input location(s) for which output is to be computed.
% RETURN Y : output location(s) corresponding to given input
% locations.
%
% SEEALSO : meanCreate
%
% COPYRIGHT : Mauricio Alvarez and Neil D. Lawrence, 2008

% MULTIGP

fhandle = str2func([meanFunction.type 'MeanCompute']);
m = fhandle(meanFunction, X, varargin{:});
