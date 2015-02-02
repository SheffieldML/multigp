function tieInd = ggwhiteMultigpTieParam(model, options)

% GGWHITEMULTIGPTIEPARAM Tie parameters for a multigp with GGWHITE kernel
% FORMAT
% DESC Tie the parameters for a multigp model that uses GG kernel.
% RETURN tieInd : cell with elements containing the indexes of parameters
% to tie.
% ARG model   : model created
% ARG options : options for tying parameters
%
% COPYRIGHT : Mauricio A. Alvarez, 2009

% MULTIGP


if isfield(options, 'typeLf') && ~isempty(options.typeLf)
else
    if isfield(options, 'isArd') && ~isempty(options.isArd)
        if options.isArd
            tieInd = tieIfArd(model, options);
        else
            tieInd = tieIfNoArd(model, options);
        end
    else
        tieInd = tieIfNoArd(model);
    end
end

function tieInd = tieIfArd(model, options)

cont = 0;
tieInd = cell(1, model.nlf);
for i = 1:model.nlf    
    cont = cont + 1;
    tieInd{cont} = paramNameRegularExpressionLookup(model, ['multi ' num2str(i) ...
        ' .* variance latent']);
end

if ~strcmp(options.tieOptions.selectMethod, 'free')
    if model.nlf > 1
        index1 = paramNameRegularExpressionLookup(model, ['multi 1 ggwhite .*' ...
            ' inverse width output 1\.']);
        index2 = paramNameRegularExpressionLookup(model, ['. ggwhite 1' ...
            ' inverse width output 1\.']);
        increm = diff(index2);
        for i = 1:model.nout
            posInd = index1(i);
            for j = 1:model.q,
                tposInd = zeros(1, model.nlf);
                tposInd(1) = posInd;
                for k =2:model.nlf,
                    tposInd(k) = tposInd(k-1) + increm(k-1);
                end
                cont = cont + 1;
                tieInd{cont} = tposInd;                
                posInd = index1(i) + j;
%                 cont = cont + 1;
%                 tieInd{cont} = paramNameRegularExpressionLookup(model, ['. ggwhite ' num2str(i)...
%                     ' inverse width output ' num2str(j) '\.']);
            end
        end        
    end
end

function tieInd = tieIfNoArd(model)

cont = 0;
tieInd = cell(1, model.nlf);
for i = 1:model.nlf
    cont = cont + 1;
    tieInd{cont} = paramNameRegularExpressionLookup(model, ['multi ' num2str(i) ...
        ' .* variance latent']);
end
if model.nlf > 1
    for i = 1:model.nout
        cont = cont + 1;
        tieInd{cont} = paramNameRegularExpressionLookup(model, ['. ggwhite ' num2str(i)...
            ' inverse width output 1.']);
    end
end
