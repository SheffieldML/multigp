function helperCreateNames(numberDims, numberInducing)

% create Names
if exist('gaussianwhiteNames.txt', 'file')
    delete('gaussianwhiteNames.txt')
end
if exist('Xunames.txt', 'file')
    delete('Xunames.txt')
end
fidGaussianwhite = fopen('gaussianwhiteNames.txt','a');
fidXu = fopen('Xunames.txt','a');
for k=1:numberInducing,
    message =['VIK ' num2str(k) ' inverse width latent'];
    fprintf(fidGaussianwhite,'%s\n',message);
end
fclose(fidGaussianwhite);
for i = 1:numberInducing
    for j = 1:numberDims
        message = ['X_u (' num2str(i) ', ' num2str(j) ')'];
        fprintf(fidXu,'%s\n',message);
    end
end
fclose(fidXu);