% To generate synthetic data, we use two conditions: real data when no
% stilumus is presented (contrast=0; called baseline data) and when a
% stimulus of contrast cVal is presented, which generates gamma
% oscillations (called stimulus data). The goal is to inject gamma bursts
% of a specified length (burstLen) in the baseline dataset such that the
% FFT of this synthetic data matches that of the stimulus data.

% The following strategy is used
% 1. All gamma bursts are injected such that they lie between the period
% denoted by the variable "stimulusPeriod (default: [0.5 2]). The centers
% of the gamma bursts lie between stimulusPeriod(1)+burstLen/2 and
% stimulusPeriod(2)-burstLen/2.

% 2. The program either injects a fixed number of bursts, specified by
% numBurstsPerTrial. If numBurstsPerTrial is empty, the average number of bursts per trial
% depends on the length of the gamma burst (burstLen) and is given by
% (stimulusPeriod/burstLen). For each trial, the number of bursts is
% obtained from a Poisson distribution with a mean of
% (stimulusPeriod/burstLen). However, we further have a condition that
% bursts do not overlap in time, so one of the overlapping bursts is
% deleted. Note that for large burstLengths, only a single burst can be
% accommodated.

% 3. For each gamma burst, the peak frequency is chosen between the
% variable "gammaRange" (default: [35 70]) following a uniform distribution

% 4. The amplitudes of the gamma bursts are obtained from a normal
% distribution with a mean given by the difference in the magnitude of the
% FFTs of the real data during stimulus and baseline periods at that
% frequency. The std is a free parameter (cvAmp*mean) because we don't
% know the actual physiological variability (as it is masked by the
% variability of the estimator). Increaaing cvAmp just makes the FFT of
% the synthetic data more noisy in the gamma range, but on average it
% should match the real data in the stimulus period.

% 5. After obtaining the gamma bursts and computing the FFTs, all the
% amplitudes are scaled by the same factor such that the FFTs of the
% synthetic and real data match in the gamma range

function [analogData,timeVals,analogData0] = generateBurstData(subjectName,expDate,protocolName,gridType,folderSourceString,electrodeNum,cVal,burstLen,cvAmp,displayFlag,synthColorName,stimulusPeriod,gammaRange,numBurstsPerTrial)

if ~exist('cvAmp','var');               cvAmp = 0.1;                    end
if ~exist('displayFlag','var');         displayFlag = 1;                end
if ~exist('synthColorName','var');      synthColorName = 'r';           end
if ~exist('stimulusPeriod','var');      stimulusPeriod = [0.5 1.5];     end
if ~exist('gammaRange','var');          gammaRange = [40 60];           end
if ~exist('numBurstsPerTrial','var');   numBurstsPerTrial=[];           end

% Get Real Data
folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'segmentedData','LFP');
st = load(fullfile(folderName,['elec' num2str(electrodeNum) 'c' num2str(cVal) '.mat']));
bl = load(fullfile(folderName,['elec' num2str(electrodeNum) 'c0.mat']));
t  = load(fullfile(folderName,'lfpInfo.mat'));
timeVals = t.timeVals;
numTrials = min(size(st.analogData,1),size(bl.analogData,1));

% Get time and frequency positions of interest
Fs = round(1/(timeVals(2)-timeVals(1)));
freqVals = 0:1/diff(stimulusPeriod):Fs-1/diff(stimulusPeriod);

tPos = intersect(find(timeVals>=stimulusPeriod(1)),find(timeVals<stimulusPeriod(2)));
fPos = intersect(find(freqVals>=gammaRange(1)),find(freqVals<gammaRange(2)));
numGammaPoints = length(fPos);

% Get Mean FFTs
mFFTst = mean(abs(fft(st.analogData(1:numTrials,tPos),[],2)));
mFFTbl = mean(abs(fft(bl.analogData(1:numTrials,tPos),[],2)));

% Get parameters for simulations
meanGammaAmps = (mFFTst(fPos)- mFFTbl(fPos));

burstSignal = zeros(numTrials,length(timeVals));
for i=1:numTrials
    burstCenterTimeList = getBurstTimes(numBurstsPerTrial,stimulusPeriod,burstLen);
    numBursts = length(burstCenterTimeList);
    
    if numBursts>0
        for j=1:numBursts
            burstCenterTime = burstCenterTimeList(j);
            fIndex = randi(numGammaPoints); % uniformly within gamma
            burstCenterFreq = freqVals(fPos(fIndex));
            burstAmplitude = meanGammaAmps(fIndex)*(1+cvAmp*randn);
            
            burst = burstAmplitude*cos(2*pi*burstCenterFreq*timeVals+2*pi*rand).*exp(-((timeVals-burstCenterTime)/(sqrt(2)*(burstLen/4))).^2);
            burstSignal(i,:) = burstSignal(i,:)+burst;
        end
    end
end

analogData0 = bl.analogData(1:numTrials,:);
analogData = analogData0 + burstSignal;
mFFTsynth = mean(abs(fft(analogData(:,tPos),[],2)));

% Rescale the burst sizes such that gamma band has the same energy
gammaPowerST = sum(mFFTst(fPos).^2);
gammaPowerSynth = sum(mFFTsynth(fPos).^2);
scalingFactor = sqrt(gammaPowerST/gammaPowerSynth);
analogData = analogData0 + scalingFactor*burstSignal;

if displayFlag
    plot(freqVals,log10(mFFTst),'k'); hold on;
    plot(freqVals,log10(mFFTbl),'g');
    plot(freqVals,log10(mean(abs(fft(analogData(:,tPos),[],2)))),'color',synthColorName);
    xlim([0 100]);
    legend('stimulus','baseline','synthetic');
end
end
function burstCenterTimeList = getBurstTimes(numBurstsPerTrial,stimulusPeriod,burstLen)

if isempty(numBurstsPerTrial)
    numBursts0 = poissrnd(diff(stimulusPeriod)/burstLen); % Intial number of bursts
else
    numBursts0 = numBurstsPerTrial;
end

burstCenterTimeRange=[(stimulusPeriod(1)+(burstLen/2)) (stimulusPeriod(2)-(burstLen/2))];    % Limits on the Mean of the Gabor signal to be injected

if numBursts0==0
    burstCenterTimeList = [];
elseif numBursts0==1
    burstCenterTimeList = burstCenterTimeRange(1) + diff(burstCenterTimeRange)*rand; %  take one time uniformly within allowed time range
else
    burstCenterTimeList0 = burstCenterTimeRange(1) + diff(burstCenterTimeRange)*rand(1,numBursts0); % create a list of numBurst0 burst times
    sortedBurstCenterTimeList0 = sort(burstCenterTimeList0); % Sort this list
    
    tPointer = sortedBurstCenterTimeList0(1);
    burstCenterTimeList(1) = tPointer; % first burst time is the earliest one
    
    while(tPointer<burstCenterTimeRange(2))
        nextGoodEntryPos = find(sortedBurstCenterTimeList0>tPointer+burstLen,1);
        if isempty(nextGoodEntryPos)
            tPointer=burstCenterTimeRange(2);
        else
            tPointer = sortedBurstCenterTimeList0(nextGoodEntryPos);
            burstCenterTimeList = cat(2,burstCenterTimeList,tPointer);
        end
    end
end
end