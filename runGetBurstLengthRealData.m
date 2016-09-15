electrodeNum = 86; %86-low, 83-medium, 29-high; 
subjectName = 'alpa'; expDate = '120316'; protocolName = 'GRF_001'; 
cVal=100; gammaFreqRangeHz=[40 60];
thresholdFraction=0.5;

[burstLengthCGT,burstLengthFeingold,burstLengthMP]=getBurstLengthRealData(subjectName,expDate,protocolName,electrodeNum,thresholdFraction,cVal,gammaFreqRangeHz);

numTrials=length(burstLengthCGT);

allBurstsCGT=[];
allBurstsFeingold=[];
allBurstsMP=[];

for i=1:numTrials
    allBurstsCGT=cat(2,allBurstsCGT,burstLengthCGT{i});
    allBurstsFeingold=cat(2,allBurstsFeingold,burstLengthFeingold{i});
    allBurstsMP=cat(1,allBurstsMP,burstLengthMP{i});
end