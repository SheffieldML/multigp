function [gKernParam, gX_u] = sparseKernGradient(model,dLdKyy, dLdKuy, dLdKuu)

% SPARSEKERNGRADIENT Computes the gradients of the parameters for sparse multigp kernel.
%
% FORMAT
% DESC computes the gradients of the parameters for sparse multigp kernel.
%  When using more than one latent function, the partial
% derivative is organized in a suitable way so that the kernGradient
% routine can be used.
% RETURN gKernParam :  derivatives wrt to the parameters of the kernels of latent,
%    output and independent functions
% RETURN gX_u : derivatives wrt to the inducing points of the latent functions
% ARG model : the model for which the gradients are to be computed.
% ARG dldKyy : the gradient of the likelihood with respect to the block or
% 	   diagonal terms.
% ARG dLdKuu :  the gradient of the likelihood with respect to the
%	   elements of K_uu.
% ARG dLdKuy : the gradient of the likelihood with respect to the
%	   elements of K_uf.
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% MULTIGP

if nargout>1
   gX_u = cell(model.nlf,1);
   for k=1:model.nlf
       gX_u{k} =  zeros(size(model.X{k}));
   end
end

if strcmp(model.kernType, 'lmc')
    A = cell(model.nlf,1);
    B = cell(model.nlf,1);
    for r =1:model.nlf
        A{r} = model.kern.comp{r}.comp{model.nlf+1}.A;
        B{r} = model.kern.comp{r}.comp{model.nlf+1}.B;
    end     
    gKernParam = zeros(1, size(model.kern.paramGroups, 1));
    startVal3 = 1;
    endVal3 = 0;
    for c = 1:length(model.kern.comp)
        endVal3 = endVal3 + model.kern.comp{c}.nParams;
        gMulti = zeros(1, size(model.kern.comp{c}.paramGroups, 1));
        startVal = 1;
        endVal = 0;
        for i = 1:model.kern.comp{c}.numBlocks
            endVal = endVal + model.kern.comp{c}.comp{i}.nParams;
            if i==c
                for r=1:model.rankCorregMatrix
                    linI = (c-1)*model.rankCorregMatrix + r;
                    covPartial = dLdKuu{linI,1};
                    % In case the derivatives wrt to X are necessary
                    if nargout>1
                        gXu = real(multiKernGradientBlockX(model.kern.comp{c}, model.X{c}, covPartial, i, i));
                        gX_u{c} = gX_u{c} + gXu;
                    end
                    gMulti(1, startVal:endVal) = gMulti(1, startVal:endVal) + ... 
                        real(multiKernGradientBlock(model.kern.comp{c}, model.X{i}, covPartial, i, i));
                end
            elseif i>model.nlf
                gBK = zeros(1,model.kern.comp{c}.comp{c}.nParams);
                gA = zeros(model.nout, model.rankCorregMatrix);
                for j=1:model.nout
                    for r=1:model.rankCorregMatrix
                        linI = (c-1)*model.rankCorregMatrix + r;
                        covPartial = dLdKuy{linI,j}';
                        gBK = gBK + A{c}(j,r)*real(multiKernGradientBlock(model.kern.comp{c}, ...
                            model.X{j+model.nlf}, model.X{c}, covPartial, c, c));
                        gA(j,r) = sum(sum(dLdKuy{linI, j}.*(model.Kyu{j,c}')));
                        if nargout>1
                            gXu = real(multiKernGradientBlockX(model.kern.comp{c},...
                               model.X{c},  model.X{j+model.nlf}, dLdKuy{linI,j}, c, c));
                            gX_u{c} = gX_u{c} + A{c}(j,r)*gXu;
                        end
                    end                    
                    if ~strcmp(model.approx, 'dtc')
                        gBK  = gBK + B{c}(j,j)*real(multiKernGradientBlock(model.kern.comp{c}, ...
                            model.X{j+model.nlf}, dLdKyy{j}, c,c));
                        if strcmp(model.approx, 'pitc')
                            temp = sum(sum(dLdKyy{j}.*model.Kyy{j}));
                        else
                            temp = (sum(diag(dLdKyy{j}).*model.Kyy{j}));
                        end
                        gA(j,:) = gA(j,:) + 2*A{c}(j,:)*temp;
                    end
               end
               gMulti(1, startVal:endVal) = [gBK (gA(:))'];
            end
            startVal = endVal + 1;
         end
        gKernParam(1, startVal3:endVal3) = gMulti;
        startVal3 = endVal3 + 1;
    end
else
    gKernParam = zeros(1, size(model.kern.paramGroups, 1));
    startVal3 = 1;
    endVal3 = 0;
    for c = 1:length(model.kern.comp)
        endVal3 = endVal3 + model.kern.comp{c}.nParams;
        gMulti = zeros(1, size(model.kern.comp{c}.paramGroups, 1));
        startVal = 1;
        endVal = 0;
        for i = 1:model.kern.comp{c}.numBlocks
            endVal = endVal + model.kern.comp{c}.comp{i}.nParams;
            if i<=model.nlf
                covPartial = dLdKuu{i,1};
                % In case the derivatives wrt to X are necessary
                if nargout>1 && i==c
                    gXu = real(multiKernGradientBlockX(model.kern.comp{c}, model.X{c}, covPartial, i, i));
                    gX_u{c} = gX_u{c} + gXu;
                end
            else
                covPartial = dLdKyy{i - model.nlf,1};
            end
            if ~isempty(covPartial)
                if isfield(model, 'useKernDiagGradient') && model.useKernDiagGradient && i>model.nlf
                    gMulti(1, startVal:endVal) = real(multiKernDiagGradientBlock(model.kern.comp{c}, model.X{i}, ...
                        covPartial, i));
                else                    
                    gMulti(1, startVal:endVal) = real(multiKernGradientBlock(model.kern.comp{c}, model.X{i}, ...
                        covPartial, i, i));
                end
            end
            startVal2 = 1;
            endVal2 = 0;
            if i > model.nlf && c <= model.nlf
                endVal2 = endVal2 + model.kern.comp{c}.comp{c}.nParams;
                [g1, g2] = multiKernGradientBlock(model.kern.comp{c}, model.X{i}, ...
                    model.X{c}, dLdKuy{c,i - model.nlf}', i, c);
                gMulti(1, startVal:endVal) = gMulti(1, startVal:endVal) + real(g1);
                gMulti(1, startVal2:endVal2) = gMulti(1, startVal2:endVal2) + real(g2);
                % In case the derivatives wrt to X are necessary
                if nargout>1
                    gXu = real(multiKernGradientBlockX(model.kern.comp{c}, model.X{i}, ...
                        model.X{c}, dLdKuy{c,i - model.nlf}, i, c));
                    gX_u{c} = gX_u{c} + gXu;
                end
            end
            startVal = endVal + 1;
        end
        gKernParam(1, startVal3:endVal3) = gMulti;
        startVal3 = endVal3 + 1;
    end
end

if nargout>1
   gX_u = cell2mat(gX_u);
end