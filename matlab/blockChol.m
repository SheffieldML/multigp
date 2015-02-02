function [L, jitter] = blockChol(model)

% BLOCKCHOL obtains a block Cholesky Factorization of a Matriz A according
% to the algorithm 

% MULTIGP
  
%   
%       http://www.netlib.org/utk/papers/factor/node9.html


dim1 = zeros(1,model.kern.numBlocks);
for i = 1:model.kern.numBlocks
    dim1(i) = size(model.X{i}, 1);
end

AA = model.K;
L = zeros(size(AA));
jitter = zeros(model.kern.numBlocks,1);

for i = 1:model.kern.numBlocks,
    A11 = AA(1:dim1(i),1:dim1(i));
    A21 = AA(dim1(i)+1:end,1:dim1(i));
    A22 = AA(dim1(i)+1:end,dim1(i)+1:end);
    [U11, jitter(i)] = jitChol(A11); 
    L11 = U11';
    invL11T = (L11')\eye(size(A11,1));
    L21 = A21*invL11T;
    AA = A22 - L21*L21'; 
    L(sum(dim1(1:i-1))+1:sum(dim1(1:i)),sum(dim1(1:i-1))+1:sum(dim1(1:i))) = L11;
    L(sum(dim1(1:i))+1:end,sum(dim1(1:i-1))+1:sum(dim1(1:i))) = L21;
end