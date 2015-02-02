function model = fixParamsMultigp(kernelType, model, options, inverseWidth, ...
    sensitivity, noisePerOutput, inverseWidthLatent)

% FIXPARAMSMULTIGP Fix the parameteres of a DTC multigp model
% FORMAT
% DESC Fix the parameters of a multigp model with DTC approximation and
% ggwhite and gaussianwhite kernels.
% ARG kernelType : type of kernel employed
% ARG model : created model
% ARG options : options for the created model
% ARG inverseWidth : inverse widths corresponding to the smoothing 
% kernels of the outputs.
% ARG sensitivity : sensitivity values for the outputs
% ARG noisePerOutput : noise added per output
% RETURN model : model with fixed parameters included
%
% COPYRIGHT : Mauricio A Alvarez, 2009

% MULTIGP 

% Fix inverse width parameters

if nargin < 7
    inverseWidthLatent = [];
end

modelAux = model;
modelAux.paramGroups = speye(size(model.paramGroups,1));
count = length(model.fix);
for i =1:model.nout
    if strcmp(options.tieOptions.selectMethod,'free')
        paramInd = paramNameRegularExpressionLookup(modelAux, ['. ' kernelType ' ' num2str(i)...
            ' inverse width output .*']);
        for k = 1:length(paramInd)
            count = count + 1;
            model.fix(count).index = paramInd(k);
            model.fix(count).value = expTransform(inverseWidth(i), 'xtoa');
        end
    else
        paramInd = paramNameRegularExpressionLookup(modelAux, ['. ' kernelType ' ' num2str(i)...
            ' inverse width output']);
        for k = 1:length(paramInd)
            count = count + 1;
            model.fix(count).index = paramInd(k);
            model.fix(count).value = expTransform(inverseWidth(i), 'xtoa');
        end
    end
end
% Fix sensitivities
for i =1:model.nout
    paramInd = paramNameRegularExpressionLookup(modelAux, ['. ' kernelType ' ' num2str(i)...
                                ' sensitivity']);    
    for k = 1:length(paramInd)
        count = count + 1;
        model.fix(count).index = paramInd(k);
        model.fix(count).value = sensitivity(i);
    end                     
end
% Fix variance parameter of the noise
% for i =1:model.nout
%     paramInd = paramNameRegularExpressionLookup(modelAux, ['Beta ' num2str(i)]);
%     count = count + 1;
%     model.fix(count).index = paramInd(1);
%     model.fix(count).value = expTransform(noisePerOutput(i), 'xtoa');
% end
for i =1:model.nout
    paramInd = paramNameRegularExpressionLookup(modelAux, ['. white ' num2str(i+1)...
                                ' variance']);
    for k = 1:length(paramInd)
        count = count + 1;
        model.fix(count).index = paramInd(k);
        model.fix(count).value = expTransform(1e-5, 'xtoa');
    end
end

if ~isempty(inverseWidthLatent)
    paramInd = paramNameRegularExpressionLookup(modelAux, '.* inverse width latent .*');
    for k = 1:length(paramInd)
        count = count + 1;
        model.fix(count).index = paramInd(k);
        model.fix(count).value = expTransform(inverseWidthLatent, 'xtoa');
    end
end

