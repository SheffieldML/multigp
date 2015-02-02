function [params, names] = ggglobalKernExtractParam(kern)

% GGGLOBALKERNEXTRACTPARAM
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% MULTIGP

params = [kern.precisionU(:)' kern.precisionG(:)' kern.sensitivity(:)'];

if nargout > 1
    sensitivityNames = cell(kern.nout, kern.nlf);
    if kern.isArd
        invWidthLatNames = cell(kern.inputDimension, kern.nlf);
        for i=1:kern.nlf
            force = num2str(i);
            for j=1:kern.inputDimension
                numdim = num2str(j);
                invWidthLatNames{i,j} = ['inverse width latent: force ' force ' dim ' numdim  '.'];
            end
        end
        if kern.tieOutputParams
            invWidthOutNames = cell(kern.inputDimension, kern.nout);
            for i=1:kern.nout
                output = num2str(i);
                for j=1:kern.inputDimension
                    numdim = num2str(j);
                    invWidthOutNames{i,j} = ['inverse width output: output ' output ' dim ' numdim  '.'];
                end
            end                        
        else
            invWidthOutNames = cell(kern.inputDimension, kern.nout, kern.nlf);
            for i=1:kern.nout
                output = num2str(i);
                for j=1:kern.nlf
                    force =  num2str(j);
                    for k=1:kern.inputDimension
                        numdim = num2str(k);
                        invWidthOutNames{k,i,j} = ['inverse width output: force ' force ' output ' output ' dim ' numdim  '.'];
                    end
                end              
            end
        end
    else
        invWidthLatNames = cell(1, kern.nlf);
        for i=1:kern.nlf
            force = num2str(i);
            invWidthLatNames{i} = ['inverse width latent: force ' force '.'];
        end
        if kern.tieOutputParams
            invWidthOutNames = cell(1, kern.nout);
            for i=1:kern.nout
                output = num2str(i);
                invWidthOutNames{i} = ['inverse width output: output ' output '.'];            
            end
        else
            invWidthOutNames = cell(1, kern.nout, kern.nlf);
            for i=1:kern.nout
                output = num2str(i);
                for j=1:kern.nlf
                    force =  num2str(j);                                        
                    invWidthOutNames{1,i,j} = ['inverse width output: force ' force ' output ' output '.'];                    
                end       
            end
        end
    end
    for i=1:kern.nout
        output = num2str(i);
        for j=1:kern.nlf
            force = num2str(j);
            sensitivityNames{i,j} = ['sensitivity output ' output ' force ' force '.'];
        end
    end
    names = [invWidthLatNames(:)' invWidthOutNames(:)' sensitivityNames(:)'];
end