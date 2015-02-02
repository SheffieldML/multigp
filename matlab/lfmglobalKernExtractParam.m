function [params, names] = lfmglobalKernExtractParam(kern)

% LFMGLOBALKERNEXTRACTPARAM
%
% COPYRIGHT

% MULTIGP

if kern.isMassFixed
    params = [kern.springVector kern.damperVector ...
        kern.inverseWidthVector kern.sensitivity(:)'];
else
    params = [kern.massVector kern.springVector kern.damperVector ...
        kern.inverseWidthVector kern.sensitivity(:)'];
end

if nargout > 1    
    namesSpring = cell(kern.nout,1);
    namesDamper = cell(kern.nout,1);
    namesInvWidth = cell(kern.nlf,1);
    for i=1:kern.nout        
        namesSpring{i} = ['spring ' num2str(i) '.'];
        namesDamper{i} = ['damper ' num2str(i) '.'];
    end
    for i=1:kern.nlf
        namesInvWidth{i} = ['inverse width ' num2str(i) '.'];        
    end 
    if ~kern.isMassFixed
        namesMass = cell(kern.nout,1);
        for i=1:kern.nout
            namesMass{i}   = ['mass ' num2str(i) '.'];            
        end        
    end
    if kern.isMassFixed
        names = [namesSpring(:)' namesDamper(:)' namesInvWidth(:)'];
    else
        names = [namesMass(:)' namesSpring(:)' namesDamper(:)' namesInvWidth(:)'];
    end
    namesSensitivity = cell(kern.nout, kern.nlf);
    for i=1:kern.nout
        output = num2str(i);
        for j=1:kern.nlf
            force = num2str(j);
            namesSensitivity{i,j} = ['sensitivity output ' output ' force ' force '.'];
        end
    end    
    names = [names(:)' namesSensitivity(:)'];
end
