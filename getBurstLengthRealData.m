function [burstLengthCGT,burstLengthFeingold,burstLengthMP,diffPower]=getBurstLengthRealData(subjectName,expDate,protocolName,electrodeNum,thresholdFraction,cVal,gammaFreqRangeHz)

gridType = 'Microelectrode'; 
folderSourceString = '';
stimulusPeriodS=[0.5 2];
baselinePeriodS=[-1.5 0];

% CGT and Feingold
cgtGaborSD = 25/1000;
filterOrderFeingold=4;

% MP
maxIteration=50;
adaptiveDictionaryParam=0.9;
dictionarySize=2500000;

% Get Data
folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'segmentedData','LFP');
st = load(fullfile(folderName,['elec' num2str(electrodeNum) 'c' num2str(cVal) '.mat']));
analogData = st.analogData;
t  = load(fullfile(folderName,'lfpInfo.mat'));
timeVals = t.timeVals;

% Compute Threshold
diffPower = getChangeInPower(analogData,timeVals,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz);
thresholdFactor = diffPower*thresholdFraction;

% Get Burst lengths CGT and Feingold
burstLengthCGT = getBurstLengthCGT(analogData,timeVals,thresholdFactor,0,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,cgtGaborSD);
burstLengthFeingold = getBurstLengthFeingold(analogData,timeVals,thresholdFactor,0,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,filterOrderFeingold);  

% Save MP results 
folderNameMain = fullfile('data',subjectName,gridType,expDate,protocolName);
tag = ['elec' num2str(electrodeNum) 'c' num2str(cVal)];

% Output folder
outputFolder = fullfile(folderNameMain,'mpAnalysis',tag);
makeDirectory(outputFolder);

% Save gaborInfo as a mat file
gaborInfoFile = fullfile(outputFolder,['gaborInfo_G' num2str(gammaFreqRangeHz(1)) '-' num2str(gammaFreqRangeHz(2)) '_M' num2str(maxIteration) ...
    '_D' num2str(100*adaptiveDictionaryParam) '_R' num2str(dictionarySize) '_S' num2str(1000*stimulusPeriodS(1)) '-' num2str(1000*stimulusPeriodS(2)) '.mat']);
if exist(gaborInfoFile,'file')
    disp(['Opening saved file ' gaborInfoFile]);
    load(gaborInfoFile); 
    burstLengthMP = getBurstLengthMP(analogData,timeVals,sqrt(thresholdFactor),0,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,maxIteration,adaptiveDictionaryParam,dictionarySize,gaborInfo,header);
else
    [burstLengthMP,~,~,gaborInfo,header] = getBurstLengthMP(analogData,timeVals,sqrt(thresholdFactor),0,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,maxIteration,adaptiveDictionaryParam,dictionarySize);
    save(gaborInfoFile,'gaborInfo','header');
end
end