function [params, names] = simglobalKernExtractParam(kern)

% SIMGLOBALKERNEXTRACTPARAM
%
% COPYRIGHT

% MULTIGP

if kern.isVarS
    params = [kern.decayVector kern.inverseWidthVector];
else
    params = [kern.decayVector kern.inverseWidthVector kern.sensitivity(:)'];
end

if nargout > 1    
    namesDecay = cell(kern.nout,1);
    namesInvWidth = cell(kern.nlf,1);
    for i=1:kern.nout
        namesDecay{i} = ['decay ' num2str(i) '.'];
    end
    for i=1:kern.nlf
        namesInvWidth{i} = ['inverse width ' num2str(i) '.'];        
    end    
    names = [namesDecay(:)' namesInvWidth(:)'];
    if ~kern.isVarS
        namesSensitivity = cell(kern.nout, kern.nlf);        
        for i=1:kern.nout
            output = num2str(i);
            for j=1:kern.nlf
                force = num2str(j);
                namesSensitivity{i,j} = ['sensitivity output ' output ' force ' force '.'];
            end
        end
    end
    names = [names(:)' namesSensitivity(:)'];
end
