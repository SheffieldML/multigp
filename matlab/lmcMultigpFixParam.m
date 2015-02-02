function model = lmcMultigpFixParam(model, options)

% LMCMULTIGPFIXPARAM Fix the parameters for a multigp model with LMC kernel
% FORMAT
% DESC Fix the parameters for a multigp model that uses LMC kernel.
% RETURN model : model with fixed parameters included
% ARG model    : model before fixing the parameters
% ARG options  : options for fixing parameters
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% MULTIGP


index = paramNameRegularExpressionLookup(model, 'multi .* variance latent');
count = 0;
for k=1:length(index);
    count = count + 1;
    model.fix(count).index = index(k);
    model.fix(count).value = expTransform(1, 'xtoa');
end

% index = paramNameRegularExpressionLookup(model, 'Beta .*');
% for k=1:length(index);
%     count = count + 1;
%     model.fix(count).index = index(k);
%     model.fix(count).value = expTransform(1e-2, 'xtoa');
% end



% if ~strcmp(model.approx, 'ftc')
%     % In the approximations fitc, pitc and dtc, this function is
%     % accomplished by the parameter beta
%     % If there is noise then is at the end and it's not necessary to look
%     % for them using paramNameRegularExpressionLookUp
%     nParamsKern = 0;
%     nToLook = options.nlf + options.includeInd;
%     for k =1:nToLook,
%         nParamsKern = nParamsKern + model.kern.comp{k}.nParams;        
%     end
%     index = (nParamsKern + options.nlf + 1):...
%         (nParamsKern + options.nlf +model.nout);     
%     for k=1:length(index);
%         count = count + 1;
%         model.fix(count).index = index(k);
%         model.fix(count).value = expTransform(1e-9, 'xtoa');
%     end    
% end