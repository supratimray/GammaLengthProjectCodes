% This part requires the full dataset
% subjectName = 'alpa'; expDate = '120316'; protocolName = 'GRF_001'; plotLocation = [1 2 3];
% subjectName = 'kesari'; expDate = '230716'; protocolName = 'GRF_002'; plotLocation = [4 5 6];
% x=load([subjectName 'MicroelectrodeRFData']);
% electrodeList = x.highRMSElectrodes;

% Data is provided for these three electrodes
subjectName = 'alpa'; expDate = '120316'; protocolName = 'GRF_001'; plotLocation = [1 2 3];
electrodeList = [29 83 86]; %86-low, 83-medium, 29-high; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cVal=100; gammaFreqRangeHz=[40 60]; thresholdFraction=0.5;

numElectrodes=length(electrodeList);
diffPower=zeros(1,numElectrodes);

methodNames{1} = 'MP'; 
methodNames{2} = 'CGT';
methodNames{3} = 'Feingold';

runMPAnalysisFlag=1;

numMethods = length(methodNames);
colorNames{1} = 'k';
colorNames{2} = [0.5 0.5 0.5];
colorNames{3} = [0.8 0.8 0.8];
symbols = 'o+X';

allBurstLengthsAllElectrodes = cell(1,numMethods);
allBurstFreqsAllElectrodes = cell(1,numMethods-1);
numBursts=zeros(numMethods,numElectrodes);
medianBurstLength=zeros(numMethods,numElectrodes);
allModsMP=[];
allNormalizedModsMP=[];

% Filter lengths that are too long
lengthLimit = 2;

for i=1:numElectrodes
    electrodeNum=electrodeList(i);
    disp(i);
    clear burstLengths freqLists
    [burstLengths{2},burstLengths{3},burstLengths{1},diffPower(i),freqLists{2},freqLists{1},modListMP]=getBurstLengthRealData(subjectName,expDate,protocolName,electrodeNum,thresholdFraction,cVal,gammaFreqRangeHz,runMPAnalysisFlag);
    
    allBurstLengthsSingleElectrodes=cell(1,numMethods);
    allBurstFreqsSingleElectrode=cell(1,numMethods);
    numTrials=length(burstLengths{1});
    
    allModsMPSingleElectrode=[];
    for iii=1:numTrials
        allModsMPSingleElectrode=cat(1,allModsMPSingleElectrode,modListMP{iii});
    end
    allModsMP=cat(1,allModsMP,allModsMPSingleElectrode);
    allNormalizedModsMP=cat(1,allNormalizedModsMP,allModsMPSingleElectrode/max(allModsMPSingleElectrode));
    
    for ii=1:numMethods
        for iii=1:numTrials
            allBurstLengthsSingleElectrodes{ii}=cat(1,allBurstLengthsSingleElectrodes{ii},burstLengths{ii}{iii}(:));
            if ii<numMethods
                allBurstFreqsSingleElectrode{ii}=cat(1,allBurstFreqsSingleElectrode{ii},freqLists{ii}{iii}(:));
            end
        end
        allBurstLengthsAllElectrodes{ii}=cat(1,allBurstLengthsAllElectrodes{ii},allBurstLengthsSingleElectrodes{ii});
        badIndices = find(allBurstLengthsSingleElectrodes{ii}>lengthLimit);

        if ii<numMethods
            allBurstFreqsAllElectrodes{ii}=cat(1,allBurstFreqsAllElectrodes{ii},allBurstFreqsSingleElectrode{ii});
        end
        burstLengthsTMP = allBurstLengthsSingleElectrodes{ii};
        burstLengthsTMP(badIndices)=[];
        numBursts(ii,i)=length(burstLengthsTMP)/numTrials;
        medianBurstLength(ii,i)=median(burstLengthsTMP);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
histC = 0.025:0.025:2;
for i=1:numMethods
    subplot(2,3,plotLocation(1));
    plot(diffPower,medianBurstLength(i,:),'color',colorNames{i},'marker',symbols(i),'linestyle','none'); hold on;
    axis([0 40 0 0.5]);

    x=allBurstLengthsAllElectrodes{i};
    badIndices = find(x>lengthLimit);
    
    subplot(2,3,plotLocation(2));
    x2=x; x2(badIndices)=[];
    histVals = hist(x2,histC);
    nHistVals = histVals/sum(histVals);
    plot(histC,nHistVals,'color',colorNames{i}); 
    hold on;
    axis([0 2 0 0.2]);
    
    disp(['Method: ' methodNames{i} ', Long trials: ' num2str(100*length(badIndices)/length(x)) '%, median burst length: ' num2str(median(medianBurstLength(i,:))) ', se: ' num2str(getSEMedian(medianBurstLength(i,:))) ', numBursts: ' num2str(mean(numBursts(i,:))) ', sem: ' num2str(std(numBursts(i,:))/sqrt(numElectrodes))]);
    
    [r,p]=corr(diffPower',medianBurstLength(i,:)','type','Spearman');
    disp(['diffPower & BurstLength: Corr=' num2str(r) ',p=' num2str(p)]);

    if i==1
        subplot(2,3,plotLocation(3));
        f=allBurstFreqsAllElectrodes{i};
        f2=f; f2(badIndices)=[];
        numUniqueFreqs=unique(f2);
        lengthUniqueFreqs = length(numUniqueFreqs);
        medianLength = zeros(1,lengthUniqueFreqs);
        seLength = zeros(1,lengthUniqueFreqs);
        
        for j=1:lengthUniqueFreqs
            pos=find(f2==numUniqueFreqs(j));
            medianLength(j) = median(x2(pos));
            if length(pos)>1
                seLength(j)=getSEMedian(x2(pos));
            else
                seLength(j)=0;
            end
        end
        errorbar(numUniqueFreqs,medianLength,seLength,'color',colorNames{i});
        axis([39 61 0 1.25]);
%         y2=allModsMP;
%         y2(badIndices)=[];
%         [r,p]=corr(x2(:),y2(:),'type','Spearman');
%         disp(['MP BurstLength Vs Amplitude, Corr=' num2str(r) ', p=' num2str(p)]);
    end
end

subplot(2,3,4);
xlabel('Change in power (dB)');
ylabel('Burst duration (s)');

subplot(2,3,5);
xlabel('Burst duration (s)');
ylabel('Fraction of bursts');

subplot(2,3,6);
xlabel('Burst center frequency (Hz)');
ylabel('Burst duration (s)');