function data = loadFxData(dataSetName);

% LOADFXDATA Load the text file with the foreign exchange rates data set and store it in a format easily accessible by Matlab.
%
% FORMAT
% DESC Loads the desired data set for the FX data from a .txt file and
% returns all the relevant data stored in a structure.
% RETURN data : Structure with all the relevant FX data.
% ARG dataSetName : The data set to load.
%
% COPYRIGHT : David Luengo, 2009
%
% SEEALSO : demSimDtcFxData, demSimwhiteDtcFxData, fxDataResults

% MULTIGP


% Determining the file name
switch dataSetName
    case 'ShortFxData2009-4outputs'
        fileName = 'fxdata2454834-2454897-4outputs.txt';
    case 'ShortFxData2009-6outputs'
        fileName = 'fxdata2454834-2454897-6outputs.txt';
    case 'ShortFxData2008-2009-6outputs'
        fileName = 'fxdata2454774-2454897-6outputs.txt';
    case 'LongFxData2007-13outputs'
        fileName = 'fxdata2454103-2454466-13outputs.txt';
    case 'LongFxData2008-2009-6outputs'
        fileName = 'fxdata2454468-2454897-6outputs.txt';
    case 'ShortFxData2009-13outputs'
        fileName = 'fxdata2454834-2454896-13outputs.txt';
    case 'LongFxData2007-2008'
        fileName = 'fxdata2454103-2454831-13outputs.txt';
    otherwise
        error('Unknown data set requested.')
end

% Reading and processing auxiliary data (first line)
[fid, message] = fopen(fileName, 'r');
if (fid < 0)
    error(['Error opening ' fileName ' : ' message]);
end

auxData = textscan(fid, '%[^\n]');
auxData = auxData{1}{1};

ind = strfind(auxData, ' ');
numExchangeRates = length(ind)-2;
data.XName = auxData(1:ind(1)-1);
data.dateFormat = auxData(ind(1)+1:ind(2)-1);
data.exchangedCurrency = cell(1, numExchangeRates);
data.baseCurrency = cell(1, numExchangeRates);
for i = 3:length(ind)-1
    exchangeRateName = auxData(ind(i)+1:ind(i+1)-1);
    ind2 = strfind(exchangeRateName, '/');
    data.baseCurrency{i-2} = exchangeRateName(1:ind2-1);
    data.exchangedCurrency{i-2} = exchangeRateName(ind2+1:end);
end
exchangeRateName = auxData(ind(end)+1:end);
ind2 = strfind(exchangeRateName, '/');
data.exchangedCurrency{numExchangeRates} = exchangeRateName(1:ind2-1);
data.baseCurrency{numExchangeRates} = exchangeRateName(ind2+1:end);

% Reading exchange data (second and subsequent lines)
status = fseek(fid, length(auxData), 'bof');
if (status < 0)
    error(['Error trying to read exchange rates from ' fileName ' : ' fileferror(fid)]);
end

exchangeData = textscan(fid, ['%u %s %s' repmat('%f', 1, numExchangeRates)]);

fclose(fid);

data.date = exchangeData{2};
data.weekDay = exchangeData{3};
data.y = cell(1, numExchangeRates);

if exist('numSamplesTrainingSet', 'var')
    data.X = exchangeData{1}(1:numSamplesTrainingSet);
    data.XTest = exchangeData{1}(numSamplesTrainingSet+1:end);
    for i = 1:numExchangeRates
        data.y{i} = exchangeData{i+3}(1:numSamplesTrainingSet);
        data.yTest{i} = exchangeData{i+3}(numSamplesTrainingSet+1:end);
    end
else
    data.X = exchangeData{1};
    data.XTest = exchangeData{1};
    for i = 1:numExchangeRates
        data.y{i} = exchangeData{i+3};
        data.yTest{i} = exchangeData{i+3};
    end    
end

return;