% DEMSENSORSIND1 Independent GPs over each temperature sensor data

% MULTIGP

rand('twister',1e6);
randn('state',1e6);

dataSetName = 'sensorsTemperature';
experimentNo = 1;

[XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData(dataSetName);

q = 1;
d = size(yTemp, 2);
nout = size(yTemp, 2);

optionsG = gpOptions('ftc');
isRbf = 1;
if isRbf
    optionsG.kern = {'rbf', 'white'};
else
    optionsG.kern = {'simwhite', 'white'};
end

optionsG.scale2var1 = 1;
itersSingleGp = 50;
displaySingleGp = 1;

X = cell(size(yTemp, 2),1);
y = cell(size(yTemp, 2),1);

for i = 1:size(yTemp, 2)
  y{i} = yTemp{i};
  X{i} = XTemp{i};
end

% Configuration of parameters
gpModel = cell(nout,1);
for k =1:nout,
    gpModel{k} = gpCreate(q, 1, X{k}, y{k}, optionsG);
    gpModel{k} = gpOptimise(gpModel{k}, 1, itersSingleGp);
end
