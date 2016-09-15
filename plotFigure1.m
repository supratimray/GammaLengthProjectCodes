function plotFigure1

subjectName = 'alpa'; expDate = '120316'; protocolName = 'GRF_001';
gridType = 'Microelectrode'; folderSourceString = '';
electrodeNum=83;cVal=100;trialNum=2;

folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'segmentedData','LFP');
st = load(fullfile(folderName,['elec' num2str(electrodeNum) 'c' num2str(cVal) '.mat']));
t  = load(fullfile(folderName,'lfpInfo.mat'));
timeVals = t.timeVals;

fRange = [0 100];
stimulusPeriod = [0.5 2];

% CGT
fRes=2; sd=0.025; 

% Feingold
gammaFreqRangeHz=[40 60];
filtOrder=4;

numTrials=size(st.analogData,1);
for i=1:numTrials
    [cgt(i,:,:),freqVals]=getCGT(st.analogData(i,:),timeVals,fRange,fRes,sd); %#ok<*AGROW>
    bpfPower(i,:)=getBPFPowerFeingold(st.analogData(i,:),timeVals,gammaFreqRangeHz,filtOrder);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fCenter=mean(gammaFreqRangeHz);
tPos = intersect(find(timeVals>=stimulusPeriod(1)),find(timeVals<stimulusPeriod(2)));

tRange = [-0.5 2];
cRange = [5 9];
colormap jet;

% CGT
subplot(321);
powerCGT = abs(cgt).^2;
pcolor(timeVals,freqVals,log10(squeeze(mean(powerCGT,1))));
shading interp;
axis([tRange fRange]);
caxis(cRange);
colorbar;

subplot(323);
pcolor(timeVals,freqVals,log10(squeeze(powerCGT(trialNum,:,:)))); 
shading interp;
axis([tRange fRange]);
caxis(cRange);
colorbar;

p=squeeze(powerCGT(:,freqVals==fCenter,tPos));
plotDataStats(subplot(3,4,9),p,trialNum);
plotDataStats(subplot(3,4,10),sqrt(p),trialNum);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Feingold
subplot(322);
plot(timeVals,mean(bpfPower,1));
xlim(tRange);

subplot(324);
plot(timeVals,bpfPower(trialNum,:)); 
xlim(tRange);

p=squeeze(bpfPower(:,tPos));
plotDataStats(subplot(3,4,11),p,trialNum);
plotDataStats(subplot(3,4,12),sqrt(p),trialNum);

end
function plotDataStats(plotHandle,powerVals,trialNum)

nBins = 100;
xLim = [0 5];
binStep = diff(xLim)/nBins;
histCenterVals = binStep:binStep:xLim(2);

powerVals = powerVals/mean(powerVals(:)); % Normalize to have a grand mean of 1

[h,c] = hist(powerVals(:),histCenterVals);
h=h/sum(h);
plot(plotHandle,c,h,'k');

hold(plotHandle,'on');
powerValsN=zeros(size(powerVals));
for i=1:size(powerVals,1)
    powerValsN(i,:) = powerVals(i,:)/mean(powerVals(i,:)); % Now mean in each trial is 1
end
[hN,cN] = hist(powerValsN(:),histCenterVals);
hN=hN/sum(hN);
plot(plotHandle,cN,hN,'r');

[h1,c1] = hist(powerVals(trialNum,:),histCenterVals);
h1=h1/sum(h1);
plot(plotHandle,c1,h1,'g');
xlim(plotHandle,xLim);
end