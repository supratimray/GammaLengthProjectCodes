% Full data is provided for these three electrodes
electrodeList = [29 83 86]; %86-low, 83-medium, 29-high; 
subjectName = 'alpa'; expDate = '120316'; protocolName = 'GRF_001'; 

cVal=100; gammaFreqRangeHz=[40 60]; thresholdFraction=0.5;

numElectrodes=length(electrodeList);
diffPower=zeros(1,numElectrodes);

methodNames{1} = 'CGT';
methodNames{2} = 'Feingold';
%methodNames{3} = 'MP';

numMethods = length(methodNames);
colorNames = 'cgr';
symbols = 'o+X';

allBurstsAllElectrodes = cell(1,numMethods);
numBursts=zeros(numMethods,numElectrodes);
medianBurstLength=zeros(numMethods,numElectrodes);
    
for i=1:numElectrodes
    electrodeNum=electrodeList(i);
    clear burstLengths
    [burstLengths{1},burstLengths{2},burstLengths{3},diffPower(i)]=getBurstLengthRealData(subjectName,expDate,protocolName,electrodeNum,thresholdFraction,cVal,gammaFreqRangeHz);
    
    allBurstsSingleElectrode=cell(1,numMethods);
    numTrials=length(burstLengths{1});
    
    for ii=1:numMethods
        for iii=1:numTrials
            allBurstsSingleElectrode{ii}=cat(1,allBurstsSingleElectrode{ii},burstLengths{ii}{iii}(:));
        end
        allBurstsAllElectrodes{ii}=cat(1,allBurstsAllElectrodes{ii},allBurstsSingleElectrode{ii});
        numBursts(ii,i)=length(allBurstsSingleElectrode{ii})/numTrials;
        medianBurstLength(ii,i)=median(allBurstsSingleElectrode{ii});
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
histC = 0.02:0.02:2;
for i=1:numMethods
    subplot(211);
    plot(diffPower,medianBurstLength(i,:),'color',colorNames(i),'marker',symbols(i),'linestyle','none'); hold on;
    axis([0 40 0 0.5]);

    subplot(212);
    histVals = hist(allBurstsAllElectrodes{i},histC);
    nHistVals = histVals/sum(histVals);
    plot(histC,nHistVals,'color',colorNames(i)); hold on;
    disp(['Method: ' methodNames{i} ', median burst length: ' num2str(median(medianBurstLength(i,:))) ', se: ' num2str(getSEMedian(medianBurstLength(i,:)))]);
end