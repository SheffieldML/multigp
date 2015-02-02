function model = multigpUpdateTopLevelParams(model)

% MULTIGPUPDATETOPLEVELPARAMS Update parameters at top level from kernel.
% FORMAT
% DESC updates parameters which are embedded in the kernel but also
% associated with the top level (such as baseline forces, decays in lfm &
% sim model).
% ARG model : model to be updated.
% RETURN model : updated model.
%
% COPYRIGHT : Mauricio Alvarez, 2008
%
% COPYRIGHT : Neil D. Lawrence, 2008
%
% SEEALSO : multigpExpandParam, multigpCreate

% MULTIGP


if strcmp(model.kernType,'lfm') || strcmp(model.kernType,'sim') ...
        || strcmp(model.kernType,'lfmwhite') || strcmp(model.kernType,'simwhite') ...
        || strcmp(model.kernType,'simnorm') 
    switch model.approx
        case 'ftc'
            if (model.includeNoise) || (model.nlf>1)
                cont = 0;
                for j=1:model.d-model.nlf
                    switch model.kernType
                        case {'lfm', 'lfmwhite'}
                            model.mu_D(j) = model.kern.comp{1}.comp{j+model.nlf}.spring;
                            for k =1:model.nlf,
                                cont = cont + 1;
                                model.S(cont) = model.kern.comp{k}.comp{j+model.nlf}.sensitivity;
                            end
                        case {'sim', 'simwhite', 'simnorm'}
                            model.mu_D(j) = model.kern.comp{1}.comp{j+model.nlf}.decay;
                            for k = 1:model.nlf
                                cont = cont + 1;
                                model.S(cont) = sqrt(model.kern.comp{k}.comp{j+model.nlf}.variance);
                            end
                        otherwise
                            %
                    end
                end
            else
                for j=1:model.d-model.nlf,
                    switch model.kernType
                        case {'lfm', 'lfmwhite'}
                            model.mu_D(j) = model.kern.comp{j+model.nlf}.spring;
                            model.S(j) = model.kern.comp{j+model.nlf}.sensitivity;
                        case {'sim', 'simwhite', 'simnorm'}
                            model.mu_D(j) = model.kern.comp{j+model.nlf}.decay;
                            model.S(j) = sqrt(model.kern.comp{j+model.nlf}.variance);
                        otherwise
                            %
                    end
                end
            end
        case {'dtc','fitc','pitc', 'dtcvar'}
            if model.includeNoise
                cont = 0;
                for j=1:model.nout,
                    switch model.kernType
                        case {'lfm', 'lfmwhite'}
                            model.mu_D(j) = model.kern.comp{1}.comp{j+model.nlf}.spring;
                            for k =1:model.nlf,
                                cont = cont + 1;
                                model.S(cont) = model.kern.comp{k}.comp{j+model.nlf}.sensitivity;
                            end
                        case {'sim', 'simwhite', 'simnorm'}
                            model.mu_D(j) = model.kern.comp{1}.comp{j+model.nlf}.decay;
                            if isfield(model, 'typeLf') && ~isempty(model.typeLf)
                                for k = 1:model.typeLf(1)
                                    cont = cont + 1;
                                    model.S(cont) = sqrt(model.kern.comp{k}.comp{j+model.nlf}.variance);
                                end
                                for k = 1+model.typeLf(1):sum(model.typeLf)
                                    cont = cont + 1;
                                    model.S(cont) = sqrt(model.kern.comp{k}.comp{j+model.nlf}.sensitivity);
                                end
                            else
                                for k = 1:model.nlf
                                    cont = cont + 1;
                                    model.S(cont) = sqrt(model.kern.comp{k}.comp{j+model.nlf}.variance);
                                end
                            end
                        otherwise
                            %
                    end
                end
            else
            end
    end
end

