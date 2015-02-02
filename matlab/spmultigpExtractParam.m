function [params, names] = spmultigpExtractParam(model, paramPart, names)

% SPMULTIGPEXTRACTPARAM Extract a parameter vector from a sparse MULTIGP model.
% FORMAT
% DESC extracts the model parameters from a structure containing
% the information about a sparse multi-output Gaussian process.
% ARG model : the model structure containing the information about
% the model.
% ARG paramPart: the parameters corresponding to the sparse part of the 
% multigp model
% RETURN params : a vector of parameters from the model.
%
% DESC extracts the model parameters from a structure containing
% the information about a sparse multi-output Gaussian process.
% ARG model : the model structure containing the information about
% the model.
% ARG paramPart: the parameters corresponding to the sparse part of the 
% multigp model
% ARG names : cell array correspondig to the names of the basic multigp
% model structure.
% RETURN params : a vector of parameters from the model.
% RETURN names : cell array of parameter names.
%
% SEEALSO : multigpCreate, multigpExpandParam, modelExtractParam
%
% COPYRIGHT : Mauricio A Alvarez, 2008
%
% MODIFICATIONS : Mauricio A Alvarez, 2009

% MULTIGP

switch model.approx
    case {'dtc','fitc', 'pitc', 'dtcvar'}
        if model.fixInducing
            params = paramPart;
            model.X_u = model.X{1};
            for k=2:model.nlf
                model.X_u = [model.X_u; model.X{k}];
            end
        else
            model.X_u = model.X{1};
            for k=2:model.nlf
                model.X_u = [model.X_u; model.X{k}];
            end
            params =  [paramPart model.X_u(:)'];            
            if nargout>1
                X_uNames = cell(size(model.X_u));
                if exist('Xunames.txt', 'file')
                    fidNames = fopen('Xunames.txt','r');
                    for i = 1:size(model.X_u, 1)
                        for j = 1:size(model.X_u, 2)
                            X_uNames{i, j} = fgetl(fidNames);
                        end
                    end
                    fclose(fidNames);
                else 
                    row = 0;
                    for nlf =1:model.nlf
                        for i = 1:size(model.X{nlf}, 1)
                            row = row + 1;
                            for j = 1:size(model.X{nlf}, 2)
                                X_uNames{row, j} = ['X_u (' num2str(i) ', ' num2str(j) ')' ' Force ' num2str(nlf) '.'];
                            end
                        end
                    end
                end
                names = {names{:}, X_uNames{:}};
            end
        end
    otherwise
        %
end