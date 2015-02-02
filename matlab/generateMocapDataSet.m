function generateMocapDataSet(subject, trainSet, testSet, setName, ...
    selectedChannels, decimation, flags);
%
% Function that reads the appropriate MOCAP data obtained from the CMU data
% base and generates a data set for the experiments.


% Parameters of the simulation

t0 = 1;               % Initial time instant
Nf = 120;               % Number of fps in the original data
if (nargin == 7)
    plotChannels = flags(1);    % Boolean: to plot or not to plot...
    generateDataSet = flags(2); % Boolean: set to true to save a new data set
    saveSkeleton = flags(3);    % Boolean: set to true to save skeleton data
else
    plotChannels = false;    % Boolean: to plot or not to plot...
    generateDataSet = true;  % Boolean: set to true to save a new data set
    saveSkeleton = false;    % Boolean: set to true to save skeleton data
end

nout = length(selectedChannels);
selChanStr = int2str(nout);

if length(decimation)==1
    decStr = int2str(decimation);
    decimation = decimation*ones(1, nout);
else
    decStr = [int2str(min(decimation)) '-' int2str(max(decimation))];    
end

subjectStr = int2str(subject);
if (trainSet < 10)
    trainSetStr = ['0' int2str(trainSet)];
else
    trainSetStr = int2str(trainSet);
end
if (testSet < 10)
    testSetStr = ['0' int2str(testSet)];
else
    testSetStr = int2str(testSet);
end

% Reading the skeleton

skel = acclaimReadSkel(strcat('./', int2str(subject), '.asf'));

if saveSkeleton
    fid = fopen(strcat(int2str(subject),'_skeleton.txt'),'w');
    fid2 = fopen(strcat(int2str(subject),'_channels.txt'),'w');
    cont = 1;
    for i=1:length(skel.tree)
        part_name = skel.tree(i).name;
        fprintf(fid,'%s',part_name);
        for k=1:length(skel.tree(i).channels)
            channel_name = skel.tree(i).channels{k};
            fprintf(fid,' %s',channel_name);
            fprintf(fid2,'Channel %2g : %s - %s',cont,part_name,channel_name);
            fprintf(fid2,'\n');
            cont = cont+1;
        end
        fprintf(fid,'\n');
    end
    st = fclose(fid);
    if st
        warning('Skeleton file not closed properly');
    end
    st = fclose(fid2);
    if st
        warning('Channels file not closed properly');
    end
end

% Loading the movement data

eval(['[ChannelsTrain, skel] = acclaimLoadChannels(''./' ...
    int2str(subject) '_' trainSetStr,  '.amc'', skel);']);
eval(['[ChannelsTest, skel] = acclaimLoadChannels(''./' ...
    int2str(subject) '_' testSetStr,  '.amc'', skel);']);

% Decimation of original channels

[nxTrain, np] = size(ChannelsTrain);
ndTrain = ceil(nxTrain./decimation);
ndTrainMax = max(ndTrain);
TTrain = zeros(ndTrainMax, nout);

[nxTest, np] = size(ChannelsTest);
ndTest = ceil(nxTest./decimation);
ndTestMax = max(ndTest);
TTest = zeros(ndTestMax, nout);

decimatedChannelsTrain = zeros(ndTrainMax, nout);
decimatedChannelsTest = zeros(ndTestMax, nout);

for i=1:nout
    TTrain(1:ndTrain(i), i) = t0 + [0:decimation(i)/Nf:decimation(i)/Nf*(ndTrain(i)-1)];
    TTest(1:ndTest(i), i) = t0 + [0:decimation(i)/Nf:decimation(i)/Nf*(ndTest(i)-1)];
    decimatedChannelsTrain(1:ndTrain(i),i) = decimate(ChannelsTrain(:,selectedChannels(i)), decimation(i));
    decimatedChannelsTest(1:ndTest(i),i) = decimate(ChannelsTest(:,selectedChannels(i)), decimation(i));
end

% Generating the data set

if generateDataSet
    nin = 1; nout = length(selectedChannels);
    Xtrain = cell(nout,1);
    Ytrain = cell(nout,1);
    Xtest = cell(nout,1);
    Ytest = cell(nout,1);
    for i=1:nout
        Xtrain{i} = TTrain(1:ndTrain(i),i);
        Ytrain{i} = decimatedChannelsTrain(1:ndTrain(i),i);
%         Ytrain{i} = Ytrain{i} - Ytrain{i}(1);
        Xtest{i} = TTest(1:ndTest(i),i);
        Ytest{i} = decimatedChannelsTest(1:ndTest(i),i);
%         Ytest{i} = Ytest{i} - Ytest{i}(1);
    end
    eval(['save ./datasetMOCAPSubject' num2str(subject) 'Train' trainSetStr ...
         'Test' testSetStr 'Dec' decStr 'SelChan' selChanStr '-' setName ...
         ' Xtrain Ytrain Xtest Ytest nin nout decimation selectedChannels']);
end

% Plotting the channels

if plotChannels
    figure(1);
    for i=1:nout
        clf;
        plot(TTrain(1:ndTrain(i),i), decimatedChannelsTrain(1:ndTrain(i),i), 'b');
        hold on
        plot(TTest(1:ndTest(i),i), decimatedChannelsTest(1:ndTest(i),i), 'r');
        title(['Output number ' num2str(i)]);
        legend(['Subject ' int2str(subject) ' - Trial ' int2str(trainSet)], ...
            ['Subject ' int2str(subject) ' - Trial ' int2str(testSet)], ...
            'Location', 'Best');
        grid on;
        pause
    end
end

return
