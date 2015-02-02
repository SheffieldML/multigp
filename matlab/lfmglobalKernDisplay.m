function lfmglobalKernDisplay(kern, spacing)

% LFMGLOBALKERNDISPLAY Display parameters of the LFMGLOBAL kernel.
% FORMAT
% DESC displays the parameters of the latent force model
% kernel and the kernel type to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : lfmKernDisplay, modelDisplay, kernDisplay
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
for i=1:kern.nlf
fprintf(spacing);
fprintf('LFM inverse width Force number %d: %2.4f (length scale %2.4f)\n', ...
        i, kern.inverseWidthVector(i), 1/sqrt(kern.inverseWidthVector(i)));    
end
for i=1:kern.nout
    fprintf(spacing);
    fprintf('\n')
    fprintf('OUTPUT NUMBER: %d\n', i)
    fprintf(spacing);
    if kern.isMassFixed
        fprintf('LFM mass: %2.4f\n', kern.massFixedVal)
    else
        fprintf('LFM mass: %2.4f\n', kern.massVector(i))
    end
    fprintf(spacing);
    fprintf('LFM spring: %2.4f\n', kern.springVector(i))
    fprintf(spacing);
    fprintf('LFM damper: %2.4f\n', kern.damperVector(i))
    for j=1:kern.nlf
        fprintf(spacing);
        fprintf('LFM sensitivity Force number %d: %2.4f\n', j, kern.sensitivity(i,j))
    end
    fprintf(spacing);
    fprintf('System Characteristics:\n')
    fprintf(spacing);
    fprintf('LFM omega: %2.4f\n', kern.omegaVector(i))
    fprintf(spacing);
    fprintf('LFM alpha: %2.4f\n', kern.alphaVector(i))
    fprintf(spacing);
    fprintf('LFM gamma: %2.4f\n', kern.gammaVector(i))
    fprintf(spacing);
    fprintf('LFM Damping Ratio: %2.4f\n', kern.zetaVector(i))
    fprintf(spacing);
    fprintf('LFM Undamped Natural Frequency: %2.4f\n', kern.omega_0Vector(i))
end
