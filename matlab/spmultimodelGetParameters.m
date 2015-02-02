function [S, D] = spmultimodelGetParameters(model, tf, dirSaveYeastResults, saveFigures)

% SPMULTIMODELGETPARAMETERS Sensitivities and decays in spmultimodel struc.
% FORMAT
% DESC
% Returns the sensitivities and decays of a spmultimodel structure and
% plots the histogram of the sensitivities, the decays and the parameter
% formed as the sensitivity over the decay.
% ARG model : the spmultimodel structure.
% ARG tf : the transcription factor for which asscociated sensitivities and
% decays are needed.
% ARG dirSaveYeastResults : path to the directory where the plots of the 
% histograms will be saved. 
% ARG saveFigures : indicates whether the plots are saved or not.
% RETURN S : values for the sensitivties associated to the 'tf'
% RETURN D : value of the decays of the outputs asscociated to the 'tf'
%
% COPYRIGHT : Mauricio A Alvarez, 2009

% MULTIGP


S = zeros(1, model.numOutGivenLat(tf));
D = zeros(1, model.numOutGivenLat(tf));



for i=1:model.numOutGivenLat(tf)
    S(i) = model.kern.comp{tf}.comp{i+1}.variance;
    D(i) = model.kern.comp{tf}.comp{i+1}.decay;
end

figure
hist(S, 18)
title('Histogram of sensitivities', 'fontsize', 20)
if saveFigures
    print('-dpng', [dirSaveYeastResults '/histS']);
end
figure
hist(D, 18)
title('Histogram of decays', 'fontsize', 20)
if saveFigures
    print('-dpng', [dirSaveYeastResults '/histD']);
end
figure
hist(S./D, 18)
title('Histogram of S/D', 'fontsize', 20)
if saveFigures
    print('-dpng', [dirSaveYeastResults '/histSD']);
end