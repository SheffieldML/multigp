function [K, asim_K] = checkKernelSymmetry( K )

% CHECKKERNELSYMMETRY Check the kernel symmetry.

% MULTIGP
  
asim_K = max(max(K-K.'));
if (asim_K ~= 0)
    K = (K + K.')/2;
end;
