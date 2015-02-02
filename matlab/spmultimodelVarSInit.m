function model = spmultimodelVarSInit(model)

% SIMMULTIMODELVARSINIT Initialize variational distribution for S
% FORMAT
% DESC Initialize the variational distribution of sensitivities for a 
% sparse multi model.
% RETURN model : model with initialized distribution.
% ARG model    : model before initializing distribution.
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% MULTIGP


for k=1:model.nout;   
    if isfield(model, 'connect') && ~isempty(model.connect)
        model.qs.mean(k,:) = model.connect(k,:);
    else
        model.qs.mean(k,:) = ones(model.nlf,1); 
    end
    %A = rand(model.nlf);
    %B = A*A';
    %model.qs.Sigma(:,:,k) = B - diag(diag(B)) + eye(model.nlf);
    model.qs.Sigma(:,:,k) = 1e-2*eye(model.nlf);
    %model.qs.Sigma(:,:,k) = zeros(model.nlf);
end

