function model = spmultigpCreate(model, options)

% SPMULTIGPCREATE
% DESC incorporates into de model options retaled with the sparse methods
% RETURN model : the structure for the sparse multigp model
% ARG model : input sparse model.
% ARG options : contains the options for the sparse multigp model

% COPYRIGHT : Mauricio Alvarez  2008

% MODIFICATIONS : Mauricio Alvarez 2009

% MULTIGP

switch options.approx
    case {'dtc','fitc', 'pitc', 'dtcvar'}
        % Sub-sample inducing variables.
        model.k = options.numActive;
        model.fixInducing = options.fixInducing;
        model.X_u = model.X{1};
        for k=2:options.nlf
            model.X_u = [model.X_u; model.X{k}];
        end        
end
if isfield(options, 'tieInducing') && ~isempty(options.tieInducing) && ...
        options.tieInducing
    effPS = 1;
else
    effPS = sum(model.k);
end
if effPS>model.N
    error('Number of active points cannot be greater than number of data.')
end
