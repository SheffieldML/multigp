function model = ggwhiteMultigpFixParam(model, options)

% GGWHITEMULTIGPFIXPARAM Fix parameters for a multigp with GGWHITE kernel
% FORMAT
% DESC Fix the parameters for a multigp model that uses GGWHITE kernel.
% RETURN model : model with fixed parameters included
% ARG model    : model before fixing the parameters
% ARG options  : options for fixing parameters
%
% COPYRIGHT : Mauricio A. Alvarez, 2009

% MULTIGP

count = 0;
% if model.nlf>1
%     base = 0.5;
% else
base = 1;
% end
for i = 1:model.nlf
    index = paramNameRegularExpressionLookup(model, ...
        ['multi [' num2str(i) '] .* variance latent']);
    for k=1:length(index);
        count = count + 1;
        model.fix(count).index = index(k);
        model.fix(count).value = expTransform(base^(i-1), 'xtoa');
    end
end
if ~strcmp(model.approx, 'ftc')
    % In the approximations fitc, pitc and dtc, this function is
    % accomplished by the parameter beta
    % If there is noise then is at the end and it's not necessary to llok
    % for them using paramNameRegularExpressionLookUp
    nParamsKern = 0;
    nToLook = options.nlf + options.includeInd;
    for k =1:nToLook,
        nParamsKern = nParamsKern + model.kern.comp{k}.nParams;        
    end
    index = (nParamsKern + options.nlf + 1):...
        (nParamsKern + options.nlf +model.nout);     
    for k=1:length(index);
        count = count + 1;
        model.fix(count).index = index(k);
        model.fix(count).value = expTransform(1e-9, 'xtoa');
    end    
%     for k = 1:model.nout
%         index = paramNameRegularExpressionLookup(model, ...
%             ['multi ' num2str(length(model.kern.comp)) ' white ' ...
%             num2str(model.nlf+k) ' .*']);
%         count = count + 1;
%         model.fix(count).index = index;
%         model.fix(count).value = expTransform(1e-9, 'xtoa');
%     end
end