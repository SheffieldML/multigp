function  gplfmToy(varargin)


% We create the xk and the fr using only the second order ODEs

seed = 1e5;
randn('seed', seed);
rand('seed', seed);

% Parameters of the simulation

nout = 3; % Number of outputs
nlf = 1; % Number of latent functions

T = 300; % Number of points of evaluation of the output functions
T2 = 300; % Number of the points of evaluation of the latent functions

NpointsTrain = 150;  
NpointsTest = 140;

m(1,1) = 1;
m(2,1) = 1;
m(3,1) = 1;
m(4,1) = 1;
m(5,1) = 1;
m(6,1) = 1;
m(7,1) = 1;
m(8,1) = 1;
m(9,1) = 1;
m(10,1) = 1;
m(11,1) = 1;


m = m(1:nout); 


D(1,1) = 0.5;
D(2,1) = 3;
D(3,1) = 1;
D(4,1) = 1;
D(5,1) = 1;
D(6,1) = 2;
D(7,1) = 2.5;
D(8,1) = 2;
D(9,1) = 2.5;
D(10,1) = 0.5;
D(11,1) = 1.4;

D = D(1:nout);

C(1,1) = 2;
C(2,1) = 1;
C(3,1) = 3;
C(4,1) = 1;
C(5,1) = 3;
C(6,1) = 1;
C(7,1) = 3;
C(8,1) = 1;
C(9,1) = 2;
C(10,1) = 1.5;
C(11,1) = 2.5;

C = C(1:nout);
% m(1,1) = 1;
% m(2,1) = 1;
% m(3,1) = 1;
% 
% D(1,1) = 1;
% D(2,1) = 1;
% D(3,1) = 1;
% 
% C(1,1) = 1;
% C(2,1) = 1.6;
% C(3,1) = 4;
% 
% 
B(1,1) = 0.5;
B(2,1) = -2;
B(3,1) = 1;
B(4,1) = 0;
B(5,1) = 0;
B(6,1) = -0.5;
B(7,1) = 0;
B(8,1) = 0;
B(9,1) = 0;
B(10,1) = 0;
B(11,1) = 0;

B = B(1:nout);

S(1,1) = 1.5;
S(1,2) = 0.5;
S(1,3) = 1;
S(1,4) = 2;
S(1,5) = 0.5;
S(1,6) = 0.2;
S(1,7) = 0.6;
S(1,8) = 1;
S(1,9) = 2;
S(1,10) = 0.5;
S(1,11) = 1;



S(2,1) = 1;%S(2,1) = 0.5;
S(2,2) = 1.2;
S(2,3) = 1.5;
S(2,4) = 2;
S(2,5) = 0.5;
S(2,6) = 0.2;
S(2,7) = 0.6;
S(2,8) = 1;
S(2,9) = 2;
S(2,10) = 0.5;
S(2,11) = 1;



S = S(:,1:nout);

% m = ones(nout,1);
% D = ones(nout,1);
% C = 1 + randn(10,1);
% B = 1 + randn(10,1);
% S = 2+rand(1, nout);


sigma2 = [0.1 0.01];
nsamp = 1;


t = linspace(0,20,T); % time points of the outputs

t2 = linspace(0,20,T2); % time points of the inputs

% If there is an independent kernel


% Parameters of the solution

alpha = C./(2*m);
omega = (sqrt(D./m-alpha.^2));


omega0 = sqrt(D./m);
zita = C./(2*sqrt(D.*m));


%  Sampling for a particular vector force

Kff = cell(nlf,1);
lfm = cell(nlf,1);
for r = 1:nlf,    
    lfm{r}.kern.inputDimension = 1;
    lfm{r}.kern = rbfKernParamInit(lfm{r}.kern);
    lfm{r}.kern.inverseWidth = sigma2(r);
end
for r =1:nlf,    
    Kff{r} = rbfKernCompute(lfm{r}.kern, t2');
    asim_Kff = max(max(Kff{r}-Kff{r}.'));
    if (asim_Kff ~= 0)
        warning(['Kf' num2str(r) 'f' num2str(r) ' contains slight asymmetries (difference between' ...
            ' off-diagonal terms smaller than %0.5g). Taking' ...
            ' (Kff+Kff.'')/2 to make it symmetric.'],asim_Kff);
        Kff{r} = (Kff{r} + Kff{r}.')/2;
    end;
    [Vff,Lff] = eig(Kff{r});
    eig_Kff = diag(Lff);
    min_Kff = min(eig_Kff);
    if (min_Kff<0)
        warning(['Kf' num2str(r) 'f' num2str(r) ' contains eigenvalues smaller than 0. Adding noise' ...
            ' with variance %0.5g to the main diagonal to make it positive' ...
            ' semidefinite'],10*abs(min_Kff));
        Kff{r} = Kff{r} + 10*abs(min_Kff)*eye(size(Kff{r}));
    end;
    Kff{r} = Kff{r} + 1e-6*eye(size(Kff{r}));
    warning(['Noise with variance 1e-6 added to the main diagonal to avoid' ...
        ' numerical problems when inverting ' 'Kf' num2str(r) 'f' num2str(r)]);
end

fr = cell(nlf,1);
for r =1:nlf
    fr{r} = gsamp(zeros(size(Kff{r},1),1),Kff{r},nsamp);
%    fr{r} =  fr{r}-mean(fr{r});
%    fr{r} = 3*(sin(10*(r)*pi*t2/3) + sin(2*(r)*pi*t2/3));
%     fr{r} = [100 zeros(1,T2-1)];
%      fr{r} = 50*randn(1,T2);
end

% Solution for xk(t) using the integral

delta_t = diff(t);      % Sampling intervals
ltf = cell(nout,1);

for tau=2:length(t),
    for k=1:nout,
        % sum of forces
        Ft = zeros(1,1:tau-1);
%         convo = zeros(1,1:tau-1);        
        for r=1:nlf
            Ft = Ft + S(r,k)*fr{r}(1:tau-1);
        end
        convo = delta_t(1:tau-1).*Ft.*exp(alpha(k)*t(1:tau-1))...
            .*sin(omega(k)*(t(tau)-t(1:tau-1)));
        ltf{k}(tau) = simpson(convo);
        clear convo Ft;
    end
end



xk = cell(nout,1);
yk = cell(nout,1);

for k=1:nout,
    xk{k} = (1/(m(k)*omega(k))*exp(-alpha(k)*t).*ltf{k})';
    xk{k} = B(k)/D(k) + (xk{k});% - mean(xk{k}));
    pond = var(xk{k}); 
    yk{k} = xk{k} + 0.1*sqrt(pond)*randn(size(xk{k},1),1);
end

close all
figure
for k =1:nout,
    subplot(nout,1,k)
    plot(t,yk{k},'LineWidth',2);    
end

figure
for r =1:nlf,
    subplot(nlf,1,r)
    plot(t2,fr{r},'LineWidth',2);    
end


randomPointSel = 1;
if randomPointSel
    index = randperm(T);    
else
    increment = floor(200/NpointsTrain);  
    indexP = 1:increment:200;
    indexNon = 1:200;
    indexNon(indexP(1:NpointsTrain)) = [];
    index = [indexP(1:NpointsTrain) indexNon];
end


includeInd = 0;

if includeInd
    signal2noiseRatio = [0.05; 0.01; 0.01];
    powSignal = zeros(nout,1); 
    for k =1:nout,
       powSignal(k) = var(yk{k}(index(1:NpointsTrain)));
    end
    senInd = powSignal.*signal2noiseRatio;
    lengthScaleInd = [1 1.2 0.9];
    inverseWidthInd = 1./(lengthScaleInd.^2);
    Kind = cell(nout,1);
    rbfInd = cell(nout,1);
    xkInd = cell(nout,1);
    for k =1:nout,
        rbfInd{k}.kern.inputDimension = 1;
        rbfInd{k}.kern = rbfKernParamInit(rbfInd{k});
        rbfInd{k}.kern.inverseWidth = inverseWidthInd(k);
        rbfInd{k}.kern.variance = senInd(k);
        Kind{k} = rbfKernCompute(rbfInd{k}.kern, t');
        xkInd{k} = gsamp(zeros(size(Kind{k},1),1),Kind{k},nsamp)';
        xk{k} = xk{k} + xkInd{k};
        yk{k} = xk{k};
    end
    figure
    for k =1:nout,
        subplot(nout,1,k)
        plot(t,xkInd{k},'LineWidth',2);
    end
    figure
    for k =1:nout,
        subplot(nout,1,k)
        plot(t,xk{k},'LineWidth',2);
    end

end


Xtrain = cell(nout,1);
Ytrain = cell(nout,1);
Xtest = cell(nout,1);
Ytest = cell(nout,1);
Ftrain = cell(nout,1);
XtrainFull = cell(nout,1);

if randomPointSel
    for i=1:nout
        Xtrain{i} = t(index(1:NpointsTrain))';
        Ytrain{i} = yk{i}(index(1:NpointsTrain));
        %   Ytrain{i} = Ytrain{i} - Ytrain{i}(1);
        Xtest{i} = t(index(NpointsTrain+1:NpointsTrain+NpointsTest))';
        Ytest{i} = yk{i}(index(NpointsTrain+1:NpointsTrain+NpointsTest));
        %  Ytest{i} = Ytest{i} - Ytest{i}(1);
        Ftrain{i} = xk{i};
        XtrainFull{i} = t';
    end
else
    for i=1:nout
        Xtrain{i} = t(index(1:NpointsTrain))';
        Ytrain{i} = yk{i}(index(1:NpointsTrain));
        %   Ytrain{i} = Ytrain{i} - Ytrain{i}(1);
        Xtest{i} = t(index(NpointsTrain+1:NpointsTrain+NpointsTest))';
        Ytest{i} = yk{i}(index(NpointsTrain+1:NpointsTrain+NpointsTest));
        %  Ytest{i} = Ytest{i} - Ytest{i}(1);
        Ftrain{i} = xk{i};
        XtrainFull{i} = t';
    end
end

Ztrain = cell(nlf,1);
Utrain = cell(nlf,1);

for i=1:nlf
    Ztrain{i} = t2';
    Utrain{i} = fr{i}';
end

save(['./data/datasetODE' num2str(nout) num2str(nlf) '_' num2str(NpointsTrain) '_' num2str(log10(seed)) 'Ind' num2str(includeInd)],'nlf','nout', 'Xtrain','Ytrain', 'Xtest','Ytest',...
    'XtrainFull', 'Ftrain', 'Utrain', 'Ztrain', 'B', 'D', 'S', 'm','C', 'sigma2');



