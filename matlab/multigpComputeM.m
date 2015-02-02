function m = multigpComputeM(model)

% MULTIGPCOMPUTEM Compute the matrix m given the model.
% FORMAT
% DESC computes the matrix m (the scaled, bias and mean function
% removed matrix of the targets), given the model.
% ARG model : the model for which the values are to be computed.
% ARG m : the scaled, bias and mean function removed values.
%
% SEEALSO : multigpCreate, multigpComputeAlpha, multigpUpdateAD
%
% COPYRIGHT : Neil D. Lawrence, 2008
%
% COPYRIGHT : Mauricio Alvarez, 2008
%
% MODIFICATIONS : David Luengo, 2009

% MULTIGP

switch model.approx
    case 'ftc'
       if sum(cellfun('length',model.X)) == model.N,
            base = 0;
        else
            base = model.nlf;
        end
        if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
            med = meanCompute(model.meanFunction, model.X, model.nlf);
            Llf = sum(cellfun('length',model.X(1:base)));
            m = model.y - med(Llf+1:end);
        else
            m = model.y;
        end
        startVal=1;
        endVal=0;
        for j=1:model.d - base
            if strcmp(model.kernType, 'lmc') && model.isIsotopic
                endVal = endVal + size(model.X{base+1}, 1);                
            else
                endVal = endVal + size(model.X{j+base}, 1);
            end
            m(startVal:endVal, 1) = m(startVal:endVal, 1) - model.bias(j+base);
            if model.scale(j)~=1
                m(startVal:endVal, 1) = m(startVal:endVal, 1)/model.scale(j+base);
            end
            startVal = endVal+1;
        end
    case {'dtc','fitc','pitc', 'dtcvar'}
        if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)   
            mu = meanCompute(model.meanFunction, mat2cell(ones(model.nout,1), ones(model.nout,1), 1), 0);
        else
            mu = zeros(model.nout,1);
        end
        m = cell(model.nout,1);
        startVal = 1;
        endVal = 0;
        for j=1:model.d            
            endVal = endVal + size(model.X{j+model.nlf}, 1);            
            m{j} = model.y(startVal:endVal) - mu(j);
            m{j} = m{j} - model.bias(j);
            if model.scale(j)~=1
                m{j} = m{j}/model.scale(j);
            end
            startVal = endVal+1;
        end
end


