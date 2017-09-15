% CV is calculated for each trial separately, by calculating the std and
% mean of power or amplitude values across time.

function plotCVDifferentMethods

subjectName = 'alpa'; expDate = '120316'; protocolName = 'GRF_001';
gridType = 'Microelectrode'; folderSourceString = '';
electrodeNum=83;cVal=100;

folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'segmentedData','LFP');
st = load(fullfile(folderName,['elec' num2str(electrodeNum) 'c' num2str(cVal) '.mat']));
t  = load(fullfile(folderName,'lfpInfo.mat'));
timeVals = t.timeVals;
dt=timeVals(2)-timeVals(1);

fRange = [0 100];    
timePeriodS=[0.5 1.5];
timePeriodB=[-1 0];

colorNames{1} = 'b'; % CGT
colorNames{2} = 'g'; % BPF
colorNames{3} = 'm'; % Multi-taper
colorNames{4} = 'k'; % WT
colorNames{5} = 'r'; % HT

% CGT & WT
fRes=1; sd=0.025; 

% Feingold & Hilbert
bandHalfWidthHz = 10;
filtOrder=4;
bandCenterFreqHz=[11 20:10:100];

% Multi-taper
params.tapers=[1 1];
params.pad=-1;
params.Fs=round(1/dt);
params.fpass=fRange;
params.trialave=0;
movingwin=[0.1 dt];

numTrials=size(st.analogData,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute PSD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tPosS=intersect(find(timeVals>=timePeriodS(1)),find(timeVals<timePeriodS(2)));
tPosB=intersect(find(timeVals>=timePeriodB(1)),find(timeVals<timePeriodB(2)));
[psdS,freqValsS]=mtspectrumc(st.analogData(:,tPosS)',params);
[psdB,freqValsB]=mtspectrumc(st.analogData(:,tPosB)',params);

subplot(231);
plot(freqValsB,log10(mean(psdB,2)),'color','m'); hold on;
plot(freqValsS,log10(mean(psdS,2)),'color','r');

subplot(234);
plot(freqValsB,std(psdB,[],2)./mean(psdB,2),'color','m'); hold on;
plot(freqValsS,std(psdS,[],2)./mean(psdS,2),'color','r');

ampB=sqrt(psdB); ampS=sqrt(psdS);
plot(freqValsB,std(ampB,[],2)./mean(ampB,2),'color','m'); hold on;
plot(freqValsS,std(ampS,[],2)./mean(ampS,2),'color','r');

plot(freqValsB,1+zeros(1,length(freqValsB)),'k');
plot(freqValsB,0.53+zeros(1,length(freqValsB)),'k');
ylim([0 2]);

%%%%%%%%%%%%%%%%%%%%%% Compute Power Vs Time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:numTrials
    % CGT
    [cgt,freqValsCGT]=getCGT(st.analogData(i,:),timeVals,fRange,fRes,sd);
    powerCGT(i,:,:)=abs(cgt).^2; %#ok<*AGROW>
    
    %WT
    [cw1,freqValsWT]=getWavelet(st.analogData(i,:),timeVals,[fRange(1)+0.1 fRange(2)],fRes); % Wavelet function wont work with 0 frequency
    powerWT(i,:,:)=abs(cw1).^2;

    % BPF (Feingold)
    for j=1:length(bandCenterFreqHz)
        powerBPF(i,j,:)=getBPFPowerFeingold(st.analogData(i,:),timeVals,[bandCenterFreqHz(j)-bandHalfWidthHz bandCenterFreqHz(j)+bandHalfWidthHz],filtOrder);
    end
    
     % BPF (Hilbert)
    for j=1:length(bandCenterFreqHz)
        powerHT(i,j,:)=getHilbertPower(st.analogData(i,:),timeVals,[bandCenterFreqHz(j)-bandHalfWidthHz bandCenterFreqHz(j)+bandHalfWidthHz],filtOrder);
    end
end

% Multi-taper
[powerMT,timeValsMT,freqValsMT]=mtspecgramc(st.analogData',movingwin,params);

x{1}=powerCGT; tList{1}=timeVals; fList{1}=freqValsCGT;
x{2}=powerBPF; tList{2}=timeVals; fList{2}=bandCenterFreqHz;
x{3}=permute(powerMT,[3 2 1]); tList{3}=timeValsMT;fList{3}=freqValsMT;
x{4}=powerWT; tList{4}=timeVals; fList{4}=freqValsWT;
x{5}=powerHT; tList{5}=timeVals; fList{5}=bandCenterFreqHz;

for i=1:length(x)
    subplot(132)
    [mCV{i},semCV{i}]=getCVStats(x{i},tList{i},timePeriodS);
    plot(fList{i},mCV{i},'color',colorNames{i}); hold on;
    plot(fList{i},mCV{i}+semCV{i},'color',colorNames{i},'linestyle','--');
    plot(fList{i},mCV{i}-semCV{i},'color',colorNames{i},'linestyle','--');
    
    subplot(133)
    [mCV{i},semCV{i}]=getCVStats(sqrt(x{i}),tList{i},timePeriodS);
    plot(fList{i},mCV{i},'color',colorNames{i}); hold on;
    plot(fList{i},mCV{i}+semCV{i},'color',colorNames{i},'linestyle','--');
    plot(fList{i},mCV{i}-semCV{i},'color',colorNames{i},'linestyle','--');
end

subplot(132)
plot(fList{1},1+zeros(1,length(fList{1})),'k');
ylim([0 2]);
xlabel('Frequency (Hz)'); ylabel('CV of power');

subplot(133)
plot(fList{1},0.53+zeros(1,length(fList{1})),'k');
ylim([0 1]);
xlabel('Frequency (Hz)'); ylabel('CV of amplitude');
end

function [mCV,semCV]=getCVStats(x,timeVals,timePeriodS)

tPos=intersect(find(timeVals>=timePeriodS(1)),find(timeVals<timePeriodS(2)));
xShort=squeeze(x(:,:,tPos));
cv = std(xShort,[],3)./mean(xShort,3);
mCV = mean(cv);
semCV=std(cv)/sqrt(size(x,1));
end
