function spmultimodelPlotResults(dataSetName, experimentNo, type, ...
    figToPlot, annotation, TransNames, options, saveFigures)

% SPMULTIMODELPLOTRESULTS Plot the results of prediction with spmultimodel.
% FORMAT
% DESC Plots the prediction for the ouputs or the latent functions. These
% predictions are obtained using the spmultimodel structure.
% ARG dataSetName : name of the dataset that contains the gene expressions
% and the names of the genes and trancription factors.
% ARG experimentNo : number of the experiment for which the model was
% trained.
% ARG type : type is either 'output' or 'latent'.
% ARG figToPlot : indicates which latent function or which output is going
% to be plotted.
% ARG annotation : containes the names of the genes.
% ARG TransNames : containes the names of the transcription factors.
% ARG options : contains the options for plotting.
% ARG saveFigures : indicates if the figures are saved or not.
%
% COPYRIGHT : Mauricio Alvarez, 2009

% MULTIGP


% load data
[xTemp, yTemp, tfNames, connect] = mapLoadData(dataSetName);


y= cell(size(yTemp,1),1);
X = cell(size(yTemp,1),1);
for k = 1:size(yTemp,1)
    y{k} = exp(yTemp(k, :))';    
    X{k} = xTemp;
end

capName = dataSetName;
capName(1) = upper(capName(1));
load(['dem' capName num2str(experimentNo) '.mat'], 'model');

if exist('localPredictionMultigp.mat', 'file')
    load('localPredictionMultigp.mat','mu','varsigma', 'Xtt');
else
    Xtt = linspace(xTemp(1)-1, xTemp(end)+1, 200)';
    [mu, varsigma] = spmultimodelPosteriorMeanVar(model, Xtt, true);
%    save('localPredictionMultigp.mat','mu','varsigma','Xtt');
end


if strcmp(type, 'latent')
    localFigToPlot = figToPlot;
    nameFigure = TransNames(figToPlot);
else
    if strcmp(type, 'output')
        localFigToPlot = figToPlot + model.nlf;
        submodel.approx = model.approx;
        nameFigure = annotation(figToPlot);
    end
end

submodel.d    = model.d;
submodel.nlf  = model.nlf;
submodel.X    = model.outX;

a = simPlotResultsOutput(submodel, X, y, Xtt ,mu, varsigma,localFigToPlot);
xlabel(options.xlabel, 'fontsize', options.fontsizeLegend)
ylabel(options.ylabel, 'fontsize', options.fontsizeLegend)
set(gca, 'ylim', options.ylim, 'xlim', options.xlim, 'fontsize', options.fontsize)
set(gca, 'position', options.position)
high = get(0, 'screensize');
set(gcf, 'position', high)
set(gcf, 'PaperPositionMode', 'auto');
%box on

if saveFigures==1    
    fileName = [dataSetName upper(model.approx)  num2str(figToPlot)];    
    print('-depsc', ['./resultsYeast/' fileName]);
    saveas(gcf,['./resultsYeast/' fileName],'fig');  
    print('-dpng', ['./resultsYeast/' fileName])   
    print('-dpdf', ['./resultsYeast/' fileName])   
end

if strcmp(type, 'latent')
    nameGene = find(strcmp(annotation, nameFigure)==1, 1);
    if ~isempty(nameGene)
        figure
        c =plot(xTemp,y{nameGene},'k');        
        set(c,   'lineWidth', 2);
        xlabel(options.gene.xlabel, 'fontsize', options.fontsizeLegend)
        ylabel(options.gene.ylabel, 'fontsize', options.fontsizeLegend)
        set(gca, 'fontname', 'arial', 'Color', 'none')
        set(gca, 'ylim', options.gene.ylim, 'xlim', options.xlim, 'fontsize', options.fontsize)
        set(gca, 'position', options.position)
        high = get(0, 'screensize');
        set(gcf, 'position', high)
        set(gcf, 'PaperPositionMode', 'auto');         
        box on
        if saveFigures==1
            fileName = strcat(fileName, 'GeneExpression');
            print('-depsc', ['./resultsYeast/' fileName]);
            saveas(gcf, ['./resultsYeast/' fileName],'fig');
            print('-dpng', ['./resultsYeast/' fileName])
            print('-dpdf', ['./resultsYeast/' fileName])  
        end
    end
end
    


