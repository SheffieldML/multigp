function [ll, X, y, noisePerOutput, kern] = generateDataset(kernName, ...
    inputDim, nout, nDataPoints, lengthScale, sensitivity, lengthScaleLatent, nlf)

% GENERATEDATASET Generates samples from a multigp model 
% FORMAT 
% DESC Generates data, sampling from a multigp model with specified kernel
% ARG kernName : name of the kernel for the multigp construction
% ARG inputDim : dimension of the input space
% ARG nout : number of outputs to be generated
% ARG nDataPoints : number of data points per output
% ARG inverseWidth : inverse width for the parameters of the Output 
% kernels
% ARG sensitivity : sensitivity value for the parameters of the Output
% kernels
% ARG inverseWidthLatent : inverse width of latent function in case gg
% kernel
% RETURN ll : likelihood given by the full GP
% RETURN X : input data sampled from a Gaussian distribution
% RETURN y : output data generated
% RETURN Kout : covariance matrix for the outputs
% RETURN noisePerOutput : noise added to each output (Default 10 percent of
% the variance).
%
% COPYRIGHT : Mauricio Alvarez, 2009

% MULTIGP

inverseWidth = 1./(lengthScale.^2);
inverseWidthLatent = 1./(lengthScaleLatent.^2);

if nargin< 7
    inverseWidthLatent = [];
end

%inverseWidth = [10 5 5 100 30 200 400 20 1 40]; % First option
%sensitivity =  [2.5  2  2.2 1.5 3.5  3   4   2  5  5];
%inverseWidth = [4   5.5 5 7 3   2 4 4.5 3.5 6.5]; % Second option
%sensitivity =  [3.5 5   2 2 3.5 4 4 8   5   5];
%nlf = 1;
d = nlf + nout;
% Sample the inputs unitl the proportion between the mean distance and the
% mean length scale for the outputs be at least 4:1
% X1 = gsamp(zeros(inputDim,1), eye(inputDim)/sqrt(inputDim), nDataPoints);
for i=0.0001:0.0001:1
    X1 = gsamp(zeros(inputDim,1), i*eye(inputDim)/sqrt(inputDim), nDataPoints);
    %X1 = unifrnd(0, i*ones(nDataPoints, inputDim));
    hh = dist2(X1,X1);
    ii = triu(sqrt(hh));
    jj = ii((ii~=0)); 
    %B = mean(jj)^2;
    B = mean(jj);
    % criteria sets the percentage between the distance between points and
    % the lengthscales for the output kernels. We have chosen 4, to
    % indicate that the distance between points must be at least 400%
    % bigger than the mean length scale for the output kernels.
    % Experimental results show that value between 100% and 600% are good
    % enough. Value over 600% might lead to numerical difficulties
    switch kernName
        case 'gg'
            outA1 = mean(2*sqrt(1./inverseWidth) + sqrt(1/inverseWidthLatent));
            outA2 = sqrt(1/inverseWidthLatent);
            if  outA2 > B %&&  outA1 > B,
            %if mean(outA) > B
                %criteria = mean(outA)/B;
                criteria1 = outA1/B;
                criteria2 = outA2/B;
                if  criteria2 > 5 && criteria2 < 15 && criteria1 > 10 %&& criteria1<30 
                    break
                end
            end
            if i == 1
                error('The input space was not properly chosen')
            end
        case 'ggwhite'
            criteria = B/mean(A(1:nout));
            if criteria > 3 && criteria < 4
                break
            end
    end
end
X = cell(1, nout+nlf);
kernType = cell(1,nlf+1);
for j=1:nlf
    kernType{j} = multigpKernComposer(kernName, d, nlf, 'ftc', j);
end
kernType{nlf+1} = multigpKernComposer('white',  d, nlf, 'ftc');
for j=1:nlf
    X{j} = zeros(1, inputDim);
end
for k=1:nout
    X{k+nlf} = X1;
end
kern = kernCreate(X,  {'cmpnd', kernType{:}});
kern.comp{1}.comp{1}.variance = 1;
noiseLatOut = zeros(nlf, nout);
noisePerOutput = zeros(1, nout);
for j=1:nlf
    for k = 1:nout,        
        kern.comp{j}.comp{nlf+k}.precisionG = inverseWidth(k);
        switch kernName
            case 'ggwhite'
                kern.comp{j}.comp{nlf+k}.variance = sensitivity(k);
            case 'gg'
                kern.comp{j}.comp{nlf+k}.sensitivity = sensitivity(k,j);
                kern.comp{j}.comp{nlf+k}.precisionU = inverseWidthLatent(j);
        end
        varOutput = kernDiagCompute(kern.comp{j}.comp{nlf+k},X1);
        noiseLatOut(j,k) = varOutput(1);
    end        
end
for j=1:nlf
    if strcmp(kernName, 'gg')
        kern.comp{j}.comp{j}.precisionU = inverseWidthLatent(j);
    end
end
for j =1:nout,
    noisePerOutput(j) = 1e-3*sum(noiseLatOut(:,j));
    kern.comp{nlf+1}.comp{nlf+j}.variance = noisePerOutput(j);
end
K = kernCompute(kern, X1);
yu = gsamp(zeros(size(K,1),1), K, 1);
u = yu(1:nlf*size(X1,1));
y = yu((nlf*size(X1,1))+1:end);
Kout = K(1+nlf*size(X1,1):end,1+nlf*size(X1,1):end);
[invKout, U, jitter] = pdinv(Kout);
if any(jitter>1e-4)
    fprintf('Warning: Added jitter of %2.4f\n', jitter)
end
logDetKout = logdet(Kout, U);
dim = size(y, 2);
ll = -dim*log(2*pi) -logDetKout - y*invKout*y';
ll = ll*0.5;
U = reshape(u,size(X1,1),nlf);
Y = reshape(y,size(X1,1),nout);
y = cell(1, nout);
for k =1:nout,
    y{k} = Y(:,k);
end










