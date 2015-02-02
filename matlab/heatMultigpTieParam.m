function tieInd = heatMultigpTieParam(model, options)

% HEATMULTIGPTIEPARAM Tie the parameters for a multigp model with HEAT kernel
% FORMAT
% DESC Tie the parameters for a multigp model that uses HEAT kernel.
% RETURN tieInd : cell with elements containing the indexes of parameters
% to tie.
% ARG model : model created
% ARG options : options for tying parameters
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% MULTIGP



tieInd = cell(1,2*(model.nlf+model.nout));
cont = 0;
for i = 1:options.nlf
    cont = cont + 1;
    tieInd{cont} = paramNameRegularExpressionLookup(model, ['multi ' num2str(i) ...
        ' .* inverse width time']);
    cont = cont + 1;
    tieInd{cont} = paramNameRegularExpressionLookup(model, ['multi ' num2str(i) ...
        ' .* inverse width space\.']);
    if isfield(options, 'kern')
        if length(options.kern) > 1
            if isfield(options.kern(i), 'includeIC') && options.kern(i).includeIC
                cont = cont + 1;
                tieInd{cont} = paramNameRegularExpressionLookup(model, ['multi ' num2str(i) ...
                    ' .* inverse width space IC\.']);
            end
        else
            if isfield(options.kern, 'includeIC') && options.kern.includeIC
                cont = cont + 1;
                tieInd{cont} = paramNameRegularExpressionLookup(model, ['multi ' num2str(i) ...
                    ' .* inverse width space IC\.']);
            end
        end
    end
end
if model.nlf > 1
    posInit   =  paramNameRegularExpressionLookup(model, '.* 1 decay');
    posNext =  paramNameRegularExpressionLookup(model, '.* 2 decay');
    increment = posNext - posInit;
    indexes = posInit - increment;
    for i = 1:model.nout
        indexes = indexes + increment;
        cont = cont+1;
        tieInd{cont} = indexes;  
    end
    posInit   =  paramNameRegularExpressionLookup(model, '.* 1 diffusion rate');
    posNext =  paramNameRegularExpressionLookup(model, '.* 2 diffusion rate');
    increment = posNext - posInit;
    indexes = posInit - increment;
    for i = 1:model.nout
        indexes = indexes + increment;
        cont = cont+1;
        tieInd{cont} = indexes;  
    end
else
    tieInd = tieInd(1:cont);
end
