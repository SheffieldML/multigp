function llApprox = highInputsGgwhite(X, y, nInducing, nOuts, ...
    iters, noisePerOutput, inverseWidth, sensitivity,...
    paramInit, methodInit, methodInitMGP, folds)

% HIGHINPUTDIMDTC Trains the DTC VAR approx for GGWHITE kernel
% FORMAT
% DESC Trains a multigp model with DTC VAR approx with fixed parameters
% ARG X : set of training inputs
% ARG y : set of training observations
% ARG nInducing : vector of pseudo points used for the approx.
% ARG nout : number of outputs
% ARG iters : number of iterations for the optimization
% ARG noisePerOutput : noise value for the parameters of the Output kernels
% ARG inverseWidth : inverse width for the parameters of the Output kernels
% ARG sensitivity : sensitivity value for the parameters of the Output kernels
% ARG paramInit : initial value for the inverse width of the VIKs
% ARG methodInit : method to initialize the positions of the pseudo points.
% ARG methodInitMGP : method to initialize the multi GP (Default = optimize
% all inducing and pseudo inputs).
% ARG folds : number of folds to repeat the experiment
% RETURN logApprox : a matrix containing the value of the bound. The
% rows indicate the fold, the columns, the index over inducing values and
% the thid dimension indicates the VIKs.
% RETURN model : trained model for the particular setup
%
% SEE ALSO : multigpCreate, muligpOptimise
%
% COPYRIGHT : Mauricio Alvarez, 2009

% MULTIGP 



llApprox = zeros(folds,length(nInducing),2);

for k=1:length(nInducing)
    nVIKs = [1 nInducing(k)];
    for j=1:length(nVIKs),
        for n =1:folds
            rand('twister',10^(n));
            randn('state',10^(n));
            llApprox(n,k,j) =  trainModelHighInGgwhite(X, y, nInducing(k), ...
                nOuts, iters, noisePerOutput,inverseWidth, sensitivity,...
                nVIKs(j), paramInit, methodInit, methodInitMGP);
        end
    end
end

