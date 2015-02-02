function [K, Kbase, n2]  = gaussianaXgaussianKernCompute(kern, x, x2)

% GAUSSIANAXGAUSSIANKERNCOMPUTE Derivative of acceleration of the Gaussian
% kernel 
% FORMAT
% DESC computes the kernel matrix for the kernel formed when taking the
% second derivative of the Gaussian kernel for the first argument. Only 
% works when x and x2 have dimension one.
% RETURN K : the kernel matrix computed at the given points.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG X : the input matrix associated with the rows of the kernel.
% ARG X2 : the input matrix associated with the columns of the kernel.
%
% FORMAT
% DESC computes the kernel matrix for the kernel formed when taking the
% second derivative of the Gaussian kernel for the first argument. Only 
% works when x have dimension one.
% RETURN K : the kernel matrix computed at the given points.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
%	
% SEEALSO : gaussianKernCompute, gaussianvXgaussianKernCompute
% 
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if size(x,2) > 1
    error('The current version of this kernel only works for 1 D inputs')
end

if nargin > 2
    if size(x2,2) > 1
        error('The current version of this kernel only works for 1 D inputs')
    end
end

if nargin < 3
    x2 = x;
end

n2 = dist2(x, x2);
Kbase = exp(-0.5*kern.precisionU*n2);

X = x(:, ones(1, size(x2,1)));
pX2 = x2';
X2 = pX2(ones(size(x,1),1), :);
X_X2 = X - X2;


K = kern.sigma2Latent*kern.precisionU*(kern.precisionU*(X_X2.^2) - 1).*Kbase;    








