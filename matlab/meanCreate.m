function  model = meanCreate(q, d, X, y, options)

% MEANCREATE creates the mean function for a multi output GP
%
% FORMAT
% DESC returns a structure for the mean function for the multiple output 
% Gaussian process model 
% RETURN model : the structure for the mean function of the multigp model
% ARG q : input dimension size.
% ARG d : output dimension size.
% ARG X : set of training inputs
% ARG y : set of training observations
% ARG options : contains the options for the MEAN of the MULTIGP model.
%
% SEE ALSO: meanCompute
%
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008 

% MULTIGP

fhandle = str2func([options.type  'MeanCreate' ]);
model = fhandle(q , d, options);
model.paramGroups = speye(model.nParams);