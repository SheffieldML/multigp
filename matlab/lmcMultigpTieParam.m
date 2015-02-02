function tieInd = lmcMultigpTieParam(model, options)

% LMCMULTIGPTIEPARAM Tie the parameters for a multigp model with LMC kernel
% FORMAT
% DESC Tie the parameters for a multigp model that uses LMC kernel.
% RETURN tieInd : cell with elements containing the indexes of parameters
% to tie.
% ARG model : model created
% ARG options : options for tying parameters
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% MULTIGP


if isfield(options.kern, 'isArd') && ~isempty(options.kern.isArd) && ...
        options.kern.isArd
    tieInd = tieIfArd(model);
else
    tieInd = tieIfNoArd(model);
end

function tieInd = tieIfArd(model)

cont = 0;
tieInd = cell(1, model.nlf*(model.q + 1));
for i = 1:model.nlf
    index1 = paramNameRegularExpressionLookup(model, ['multi ' num2str(i) ...
            ' .* inverse width latent 1\.']);
    %increm = diff(index1);     
    posLat =  index1(1) - 1;
    for j = 1:model.q
        posInit = index1(2) -1 + j;
        %posEnd = posInit + increm*(model.nout-1);
        cont = cont + 1;
        %tieInd{cont} = [(posLat+j) posInit:increm(2):posEnd];
        tieInd{cont} = [(posLat+j) posInit];
    end
    cont = cont + 1;
    tieInd{cont} = paramNameRegularExpressionLookup(model, ['multi ' num2str(i) ...
        ' .* variance latent']);
end

function tieInd = tieIfNoArd(model)

cont = 0;
tieInd = cell(1, 2*model.nlf);
for i = 1:model.nlf
    cont = cont + 1;
    tieInd{cont} = paramNameRegularExpressionLookup(model, ['multi ' num2str(i) ...
        ' .* inverse width latent 1.']);
    cont = cont + 1;
    tieInd{cont} = paramNameRegularExpressionLookup(model, ['multi ' num2str(i) ...
        ' .* variance latent']);
end

% cont = cont + 1;
% tieInd{cont} = paramNameRegularExpressionLookup(model, 'multi .* inverse width latent 1.');
% cont = cont + 1;
% tieInd{cont} = paramNameRegularExpressionLookup(model, 'multi .* variance latent');
% 
% indexT = zeros(model.nlf, model.k(1));
% for i=1:model.nlf
%     indexT(i,:) = paramNameRegularExpressionLookup(model, ['X_u .* Force ' num2str(i) '\.']);
% end
% 
% for i=1:model.k(1)
%     cont = cont + 1;
%     tieInd{cont} = indexT(:,i)';
% end
% 

    
