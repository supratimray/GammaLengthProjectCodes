function testPerformanceSynthData(thresholdFractionList,electrodeNum,numMeanBursts)

subjectName = 'alpa'; expDate = '120316'; protocolName = 'GRF_001';
gridType = 'Microelectrode'; folderSourceString = ''; cVal=100;

numThresholds = length(thresholdFractionList);

% BurstDataParameters
burstLenList = [0.05 0.1:0.1:1];
cvAmp=0.1;
stimulusPeriodS=[0.5 2];
baselinePeriodS=[-1.5 0];
gammaFreqRangeHz=[40 60];
numBurstLengths = length(burstLenList);
synthColorList = jet(numBurstLengths);
displayFlagBurst=0;

% CGT (Xing et al., 2014)
cgtGaborSDList = [12.5 25]/1000;
numCGTSDList = length(cgtGaborSDList);
displayFlagCGT=0;

% Feingold et al., 2015
filterOrderFeingold=4;
displayFlagFeingold=0;

% MP
maxIteration=50; %#ok<*NASGU>
adaptiveDictionaryParam=0.9;
dictionarySize=2500000;
displayFlagMP=0;
showMPResults=1;

% Initialize
medianBurstLengthCGT = zeros(numBurstLengths,numCGTSDList,numThresholds);
seBurstLengthCGT = zeros(numBurstLengths,numCGTSDList,numThresholds);

medianBurstLengthFeingold = zeros(numBurstLengths,numThresholds);
seBurstLengthFeingold = zeros(numBurstLengths,numThresholds);

medianBurstLengthMP = zeros(numBurstLengths,numThresholds);
seBurstLengthMP = zeros(numBurstLengths,numThresholds);

for i=1:numBurstLengths
    burstLen=burstLenList(i);
    synthColorName = synthColorList(i,:);
    disp(['Burst Length: ' num2str(burstLen)]);
    
    [analogData,timeVals] = generateBurstData(subjectName,expDate,protocolName,gridType,folderSourceString,electrodeNum,cVal,burstLen,cvAmp,displayFlagBurst,synthColorName,stimulusPeriodS,gammaFreqRangeHz,numMeanBursts);
    diffPower = getChangeInPower(analogData,timeVals,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz);
    
    for ii=1:length(thresholdFractionList)
        thresholdFactor = diffPower*thresholdFractionList(ii);
        
        % Estimate burst length using CGT (Xing et al., 2012)
        for j=1:numCGTSDList
            burstLengthCGT = getBurstLengthCGT(analogData,timeVals,thresholdFactor,displayFlagCGT,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,cgtGaborSDList(j));
            [medianBurstLengthCGT(i,j,ii),seBurstLengthCGT(i,j,ii)] = getMedianAndSE(burstLengthCGT);
        end
        
        % Estimate burst length using Feingold et al., 2015
        burstLengthFeingold = getBurstLengthFeingold(analogData,timeVals,thresholdFactor,displayFlagFeingold,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,filterOrderFeingold);
        [medianBurstLengthFeingold(i,ii),seBurstLengthFeingold(i,ii)] = getMedianAndSE(burstLengthFeingold);
        
        % Estimate burst length using Stochastic MP
        if showMPResults
            if ii==1
                [burstLengthMP,~,~,gaborInfo,header] = getBurstLengthMP(analogData,timeVals,sqrt(thresholdFactor),displayFlagMP,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,maxIteration,adaptiveDictionaryParam,dictionarySize); %#ok<*UNRCH>
            else
                burstLengthMP = getBurstLengthMP(analogData,timeVals,sqrt(thresholdFactor),displayFlagMP,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,maxIteration,adaptiveDictionaryParam,dictionarySize,gaborInfo,header);
            end
            [medianBurstLengthMP(i,ii),seBurstLengthMP(i,ii)] = getMedianAndSE(burstLengthMP);
        end
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XingDataX = [0.025 0.05 0.075 0.1 0.125];
XingDataY = [0.115 0.122 0.14 0.155 0.165];

for ii=1:numThresholds
    figure(ii);
    plot(XingDataX,XingDataY,'ks','linewidth',2); hold on;
    
    colorNames = jet(numCGTSDList+3);
    for i=1:numCGTSDList
        errorbar(burstLenList,squeeze(medianBurstLengthCGT(:,i,ii)),squeeze(seBurstLengthCGT(:,i,ii)),'color',colorNames(i,:));
    end
    
    errorbar(burstLenList,medianBurstLengthFeingold(:,ii),seBurstLengthFeingold(:,ii),'color',colorNames(numCGTSDList+1,:));
    
    if showMPResults
        errorbar(burstLenList,medianBurstLengthMP(:,ii),seBurstLengthMP(:,ii),'color',colorNames(numCGTSDList+3,:));
    end
    
    %
    xs=0:0.1:1;
    plot(xs,xs,'k');
    axis([0 1 0 1]);
end
end

function [medianAllBurstLength,seAllBurstLength,medianFirstBurstLength,seFirstBurstLength] = getMedianAndSE(lengthList)

allLengths=[]; firstLengths=[];
for j=1:length(lengthList)
    allLengths=cat(1,allLengths,lengthList{j}(:));
    if ~isempty(lengthList{j})
        firstLengths=cat(1,firstLengths,lengthList{j}(1));
    end
end
if ~isempty(allLengths)
    medianAllBurstLength = median(allLengths);
    seAllBurstLength = getSEMedian(allLengths);
    medianFirstBurstLength = median(firstLengths);
    seFirstBurstLength = getSEMedian(firstLengths);
else
    medianAllBurstLength = 0;
    seAllBurstLength = 0;
    medianFirstBurstLength = 0;
    seFirstBurstLength = 0;
end
end