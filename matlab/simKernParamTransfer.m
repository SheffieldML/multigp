function simKern = simKernParamTransfer(kern, simKern, whichOutput, whichLatent)

% SIMKERNPARAMTRANSFER
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% MULTIGP

simKern.decay = kern.decayVector(whichOutput);
simKern.inverseWidth = kern.inverseWidthVector(whichLatent);

if ~kern.isVarS
    if kern.isNegativeS        
        simKern.sensitivity = kern.sensitivity(whichOutput, whichLatent);
    else
        simKern.variance = kern.sensitivity(whichOutput, whichLatent);
    end
end
