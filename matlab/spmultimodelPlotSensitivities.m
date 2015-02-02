function spmultimodelPlotSensitivities(dataSetName, figToPlot, approx, tf,  options, saveFigures)

% SPMULTIMODELPLOTSENSITIVTIES Plot the histograms for sensitivities.

% MULTIGP

load(['S2' tf '.mat']);
varS = -1./diag(gS2);
stdS = sqrt(varS)';
SNR = S./stdS;
figure
hist(SNR, 15, 'Color', 'none')
h = findobj(gca,'Type','patch');
set(h,'FaceColor','k','EdgeColor','w')
g = xlabel(options.xlabel, 'fontsize', options.fontsizeLegend);
prop = get(g);
posXlabel = prop.Position;
posXlabel(2) = posXlabel(2) -1;  
xlabel(options.xlabel, 'fontsize', options.fontsizeLegend, 'Position', posXlabel);
ylabel(options.ylabel, 'fontsize', options.fontsizeLegend)
set(gca, 'fontname', 'arial')
set(gca, 'ylim', options.ylim, 'xlim', options.xlim, 'fontsize', options.fontsize)
set(gca, 'position', options.position)
high = get(0, 'screensize');
set(gcf, 'position', high)
set(gcf, 'PaperPositionMode', 'auto');
box on

if saveFigures==1    
    fileName = [dataSetName 'Hist' approx num2str(figToPlot)];    
    print('-depsc', ['./resultsYeast/' fileName]);
    saveas(gcf,['./resultsYeast/' fileName],'fig');  
    print('-dpng', ['./resultsYeast/' fileName])   
    print('-dpdf', ['./resultsYeast/' fileName])   
end