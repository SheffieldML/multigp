function tieInd = ggMultigpTieParam(model, options)

% GGMULTIGPTIEPARAM Tie the parameters for a multigp model with GG kernel
% FORMAT
% DESC Tie the parameters for a multigp model that uses GG kernel.
% RETURN tieInd : cell with elements containing the indexes of parameters
% to tie.
% ARG model : model created
% ARG options : options for tying parameters
%
% COPYRIGHT : Mauricio A. Alvarez, 2009
%
% MODIFICATIONS : Mauricio A. Alvarez, 2010

% MULTIGP

if isfield(options, 'icmStyle') && ~isempty(options.icmStyle) && ...
        options.icmStyle
    tieInd = tieIfICM(model, options);
else
    if isfield(options, 'typeLf') && ~isempty(options.typeLf)
    else
        if isfield(options, 'isArd') && ~isempty(options.isArd)
            if options.isArd
                tieInd = tieIfArd(model, options);
            else
                tieInd = tieIfNoArd(model, options);
            end
        else
            tieInd = tieIfNoArd(model, options);
        end
    end
end

function tieInd = tieIfICM(model, options)

if isfield(options, 'isArd') && ~isempty(options.isArd) && options.isArd
    cont = 0;
    tieInd = cell(1, model.q + 1);
    for i=1:model.q
        index1 = paramNameRegularExpressionLookup(model, ['multi .* inverse width .* ' num2str(i) '\.']);
        cont = cont + 1;
        tieInd{i} = index1;        
    end    
    cont = cont + 1;
    tieInd{cont} = paramNameRegularExpressionLookup(model, 'multi  .* variance latent');
else
    tieInd{1} = paramNameRegularExpressionLookup(model, 'multi .* inverse width  .* 1\.');
    tieInd{2} = paramNameRegularExpressionLookup(model, 'multi  .* variance latent');
end


function tieInd = tieIfArd(model, options)

cont = 0;
tieInd = cell(1, model.nlf*(model.q + 1));
for i = 1:model.nlf
    index1 = paramNameRegularExpressionLookup(model, ['multi ' num2str(i) ...
            ' .* inverse width latent 1\.']);
    increm = diff(index1);     
    %posInit = index1(2) -1;
    posLat =  index1(1) - 1;
    for j = 1:model.q
        posInit = index1(2) -1 + j;%posInit + j;
        posEnd = posInit + increm(2)*(model.nout-1);
        cont = cont + 1;
        tieInd{cont} = [(posLat+j) posInit:increm(2):posEnd];
%         tieInd{cont} = paramNameRegularExpressionLookup(model, ['multi ' num2str(i) ...
%             ' .* inverse width latent ' num2str(j) '\.']);    
    end
    cont = cont + 1;
    tieInd{cont} = paramNameRegularExpressionLookup(model, ['multi ' num2str(i) ...
        ' .* variance latent']);
end

if ~strcmp(options.tieOptions.selectMethod, 'free')
    if model.nlf > 1
        index1 = paramNameRegularExpressionLookup(model, ['multi 1 gg .*' ...
            ' inverse width output 1\.']);
        index2 = paramNameRegularExpressionLookup(model, ['. gg 1' ...
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
            end
        end        
    end
end

function tieInd = tieIfNoArd(model, options)

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
if ~strcmp(options.tieOptions.selectMethod, 'free')    
    if model.nlf > 1
        index1 = paramNameRegularExpressionLookup(model, ['multi 1 gg .*' ...
            ' inverse width output 1\.']);
        index2 = paramNameRegularExpressionLookup(model, ['. gg 1' ...
            ' inverse width output 1\.']);
        increm = diff(index2);
        for i = 1:model.nout
            posInd = index1(i);
            tposInd = zeros(1, model.nlf);
            tposInd(1) = posInd;
            for k =2:model.nlf,
                tposInd(k) = tposInd(k-1) + increm(k-1);
            end
            cont = cont + 1;
            tieInd{cont} = tposInd;
        end
%         for i = 1:model.nout
%             cont = cont + 1;
%             tieInd{cont} = paramNameRegularExpressionLookup(model, ['. gg ' num2str(i)...
%                 ' inverse width output 1.']);
%         end
    end
end




