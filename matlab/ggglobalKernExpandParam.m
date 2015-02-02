function kern = ggglobalKernExpandParam(kern, params)

% GGGLOBALKERNEXPANDPARAM
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% MULTIGP

if kern.isArd
    nParamsLat = kern.inputDimension*kern.nlf;     
    kern.precisionU = reshape(params(1:nParamsLat), kern.inputDimension, kern.nlf);    
    if kern.tieOutputParams
        nParamsOut = kern.inputDimension*kern.out;
        kern.precisionG = reshape(params(nParamsLat+1:nParamsLat+nParamsOut), kern.inputDimension, kern.nout);    
    else
        nParamsOut = kern.inputDimension*kern.out*kern.nlf;
        kern.precisionG = reshape(params(nParamsLat+1:nParamsLat+nParamsOut), kern.inputDimension, kern.nout, kern.nlf);    
    end    
else
    nParamsLat = kern.nlf;
    kern.precisionU = reshape(params(1:nParamsLat), 1, kern.nlf);  
    if kern.tieOutputParams
        nParamsOut = kern.nout;
        kern.precisionG = reshape(params(nParamsLat+1:nParamsLat+nParamsOut), 1, kern.nout);    
    else
        nParamsOut = kern.nout*kern.nlf;
        kern.precisionG = reshape(params(nParamsLat+1:nParamsLat+nParamsOut), 1, kern.nout, kern.nlf);    
    end    
end

kern.sensitivity = reshape(params(nParamsLat+nParamsOut+1:end), kern.nout, kern.nlf);