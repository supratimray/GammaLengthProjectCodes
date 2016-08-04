function testPerformanceSynthData(thresholdFactor)
subjectName = 'alpa'; expDate = '120316'; protocolName = 'GRF_001';
gridType = 'Microelectrode'; folderSourceString = '';
electrodeNum=83;cVal=100;
% thresholdFactor = 10;

% BurstDataParameters
burstLenList = [0.05 0.1 0.2 0.3 0.4 0.5];
cvAmp=0.1;
stimulusPeriodS=[0.5 1.5];
baselinePeriodS=[-1 0];
gammaFreqRangeHz=[40 60];
numBurstLengths = length(burstLenList);
synthColorList = jet(numBurstLengths);
displayFlagBurst=0;
numMeanBursts=1;

% CGT (Xing et al., 2014)
cgtGaborSDList = [6.25 12.5 25]/1000;
numCGTSDList = length(cgtGaborSDList);
displayFlagCGT=0;

% Feingold et al., 2015
filterOrderFeingold=4;
displayFlagFeingold=0;

% MP
maxIteration=50; %#ok<*NASGU>
adaptiveDictionaryParam=0.9;
dictionarySize=[];
displayFlagMP=0;
showMPResults=1;

% Initialize
medianBurstLengthCGT = zeros(numBurstLengths,numCGTSDList);
seBurstLengthCGT = zeros(numBurstLengths,numCGTSDList);

medianBurstLengthFeingold = zeros(1,numBurstLengths);
seBurstLengthFeingold = zeros(1,numBurstLengths);

medianBurstLengthMP = zeros(1,numBurstLengths);
seBurstLengthMP = zeros(1,numBurstLengths);

for i=1:numBurstLengths
    burstLen=burstLenList(i);
    synthColorName = synthColorList(i,:);
    disp(['Burst Length: ' num2str(burstLen)]);
    
    [analogData,timeVals] = generateBurstData(subjectName,expDate,protocolName,gridType,folderSourceString,electrodeNum,cVal,burstLen,cvAmp,displayFlagBurst,synthColorName,stimulusPeriodS,gammaFreqRangeHz,numMeanBursts);
    
    % Estimate burst length using CGT (Xing et al., 2012)
    for j=1:numCGTSDList
        burstLengthCGT = getBurstLengthCGT(analogData,timeVals,thresholdFactor,displayFlagCGT,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,cgtGaborSDList(j));
        [medianBurstLengthCGT(i,j),seBurstLengthCGT(i,j)] = getMedianAndSE(burstLengthCGT);
    end
    
    % Estimate burst length using Feingold et al., 2015
    burstLengthFeingold = getBurstLengthFeingold(analogData,timeVals,thresholdFactor,displayFlagFeingold,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,filterOrderFeingold);
    [medianBurstLengthFeingold(i),seBurstLengthFeingold(i)] = getMedianAndSE(burstLengthFeingold);
    
    % Estimate burst length using Stochastic MP
    if showMPResults
        burstLengthMP = getBurstLengthMP(analogData,timeVals,sqrt(thresholdFactor),displayFlagMP,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,maxIteration,adaptiveDictionaryParam,dictionarySize); %#ok<*UNRCH>
        [medianBurstLengthMP(i),seBurstLengthMP(i)] = getMedianAndSE(burstLengthMP);
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XingDataX = 2*[0.025 0.05 0.075 0.1 0.125];
XingDataY = [0.115 0.122 0.14 0.155 0.165];

plot(XingDataX,XingDataY,'ks','linewidth',2); hold on;

colorNames = jet(numCGTSDList+3);
for i=1:numCGTSDList
    errorbar(burstLenList,medianBurstLengthCGT(:,i),seBurstLengthCGT(:,i),'color',colorNames(i,:));
end

errorbar(burstLenList,medianBurstLengthFeingold,seBurstLengthFeingold,'color',colorNames(numCGTSDList+1,:));

if showMPResults
    errorbar(burstLenList,medianBurstLengthMP,seBurstLengthMP,'color',colorNames(numCGTSDList+3,:));
end

% 
xs=0:0.1:1;
plot(xs,xs,'k');
axis([0 1 0 1]);
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