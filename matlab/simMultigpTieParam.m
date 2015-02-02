function tieInd = simMultigpTieParam(model, options)

% SIMMULTIGPTIEPARAM Tie the parameters for a multigp model with SIM kernel
% FORMAT
% DESC Tie the parameters for a multigp model that uses SIM kernel.
% RETURN tieInd : cell with elements containing the indexes of parameters
% to tie.
% ARG model : model created
% ARG options : options for tying parameters
%
% COPYRIGHT : Mauricio A. Alvarez, David Luengo 2009

% MULTIGP


if isfield(options, 'typeLf') && ~isempty(options.typeLf)
    for i = 1:options.typeLf(1)
        tieInd{i} = paramNameRegularExpressionLookup(model, ['multi ' num2str(i) ...
            ' .* inverse width']);
    end
    for i = 1+options.typeLf(1):sum(options.typeLf)
        tieInd{end+1} = paramNameRegularExpressionLookup(model, ['multi ' num2str(i) ...
            ' .* variance']);
    end
    if strcmp(options.tieOptions.selectMethod, 'typeLf')
        % Tying separately the decays for SIM and SIM-WHITE kernels
        % (i.e. for each output there is a separate decay for the
        % SIM kernels and for the SIM-WHITE kernels)
        for i = 1:model.nout
            tieInd{end+1} = paramNameRegularExpressionLookup(model, ...
                ['sim ' num2str(i)  ' decay']);
            tieInd{end+1} = paramNameRegularExpressionLookup(model, ...
                ['simwhite ' num2str(i)  ' decay']);
        end
    else
        for i = 1:model.nout
            tieInd{end+1} = paramNameRegularExpressionLookup(model, ['.* ' num2str(i)  ' decay']);
        end
    end
else
    tieInd = cell(1,model.nlf+model.nout);
    cont = 0;
    for i = 1:options.nlf
        cont = cont + 1;
        tieInd{cont} = paramNameRegularExpressionLookup(model, ['multi ' num2str(i) ...
            ' .* inverse width']);   
    end
    if ~(isfield(options, 'isSubmodel') && ~isempty(options.isSubmodel))
        if model.nout > 1
            posInit   =  paramNameRegularExpressionLookup(model, '.* 1 decay');
            posNext =  paramNameRegularExpressionLookup(model, '.* 2 decay');
            increment = posNext - posInit;
            indexes = posInit - increment;
            for i = 1:model.nout
                indexes = indexes + increment;
                cont = cont+1;
                tieInd{cont} = indexes;
                % tieInd{cont} = paramNameRegularExpressionLookup(model, ['.* ' num2str(i)  ' decay']);
            end
        else
            cont = cont + 1;
            tieInd{cont} = paramNameRegularExpressionLookup(model, '.* 1 decay');
        end
    else
        tieInd = tieInd(1:model.nlf);
    end
end