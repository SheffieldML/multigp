function model = initInducingMultigp(kernelType, model, options, inverseWidth, paramInit)

% INITINDUCINGMULTIGP Initiliaze the inducing parameteres of a DTC VAR multigp
% FORMAT
% DESC Initialize the inducing parameters of a DTC VAR multigp model.
% ARG kernelType : type of kernel
% ARG model : created model
% ARG options : options for the created model
% ARG inverseWidth : inverse widths corresponding to the smoothing
% kernels of the outputs.
% ARG paramInit : initial value for the inducing kernel parameters
% RETURN model : model with fixed parameters included
%
% COPYRIGHT : Mauricio A Alvarez, 2009

% MULTIGP


switch kernelType

    case 'ggwhite'

        if isscalar(paramInit)
            paramInit = paramInit*ones(1,options.nVIKs);
        end
        fixingLengthScale = 0;
        if fixingLengthScale
            for i=1:model.nlf
                paramInd = paramNameRegularExpressionLookup(modelAux, ...
                    ['multi ' num2str(i) ' gaussianwhite .* inverse width latent .*']);
                for k = 1:length(paramInd)
                    count = count + 1;
                    model.fix(count).index = paramInd(k);
                    model.fix(count).value = expTransform(inverseWidth(i), 'xtoa');
                end
            end
        else
            if strcmp(options.tieOptions.selectMethod,'free')
                params =  modelExtractParam(model);
                paramInd = paramNameRegularExpressionLookup...
                    (model, '. gaussianwhite .* inverse width latent .*');
                params(paramInd) = log(mean(inverseWidth(1:nout))); % The value of the last simulation
                model = modelExpandParam(model, params);
            else
                params =  modelExtractParam(model);
                if isempty(options.nVIKs)
                    paramInd = paramNameRegularExpressionLookup...
                        (model, '. gaussianwhite .* VIK .* inverse width latent');
                    params(paramInd) = log(rand(1, options.numActive));
                else
                    startVal = 1;
                    interVal = options.numActive/options.nVIKs;
                    if options.nlf > 1
                        for k=1:options.nVIKs,
                            paramInd = paramNameRegularExpressionLookup(model,...
                                ['. gaussianwhite .* VIK ' num2str(startVal) ...
                                'inverse width latent']);
                            params(paramInd) = log(paramInit(k));
                            startVal = startVal + interVal;
                        end
                    else
                        if interVal ~=1
                            for k=1:options.nVIKs,
                                params(k) = log(paramInit(k));
                            end
                        else
                            params(1:options.nVIKs) = log(paramInit);
                        end
                    end
                end
                model = modelExpandParam(model, params);
            end
        end
    case 'gg'

        fixingLengthScale = 0;
        if fixingLengthScale
            for i=1:model.nlf
                paramInd = paramNameRegularExpressionLookup(modelAux, ...
                    ['multi ' num2str(i) ' gaussian .* inverse width latent .*']);
                for k = 1:length(paramInd)
                    count = count + 1;
                    model.fix(count).index = paramInd(k);
                    model.fix(count).value = expTransform(inverseWidth(i), 'xtoa');
                end
            end
        else
            if strcmp(options.tieOptions.selectMethod,'free')
                params =  modelExtractParam(model);
                paramInd = paramNameRegularExpressionLookup...
                    (model, '. gaussian .* inverse width latent .*');
                params(paramInd) = log(mean(inverseWidth(1:nout))); % The value of the last simulation
                model = modelExpandParam(model, params);
            else
                params =  modelExtractParam(model);
                paramInd = paramNameRegularExpressionLookup...
                    (model, '. gaussian .* inverse width latent');
                params(paramInd) = log(paramInit);
                model = modelExpandParam(model, params);
            end
        end
end



