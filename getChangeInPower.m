% This program finds the increase in power in the stimulus period relative
% to baseline in the gamma band. This is used to find the appropriate
% threshold for each electrode.

function diffPower = getChangeInPower(analogData,timeVals,stimulusPeriod,baselinePeriod,gammaRange)

if diff(stimulusPeriod) ~= diff(baselinePeriod)
    error('Stimulus and baseline durations should be the same..');
end

tPosSt = intersect(find(timeVals>=stimulusPeriod(1)),find(timeVals<stimulusPeriod(2)));
tPosBl = intersect(find(timeVals>=baselinePeriod(1)),find(timeVals<baselinePeriod(2)));

% Get Mean FFTs
mFFTst = mean(abs(fft(analogData(:,tPosSt),[],2)));
mFFTbl = mean(abs(fft(analogData(:,tPosBl),[],2)));

Fs = round(1/(timeVals(2)-timeVals(1))); 
fAxSt=0:Fs/length(mFFTst):Fs-Fs/length(mFFTst);

% Gamma Increase
gammaPos = intersect(find(fAxSt>=gammaRange(1)),find(fAxSt<gammaRange(2)));
diffFFT = mFFTst(gammaPos)./mFFTbl(gammaPos);
diffPower = mean(diffFFT.^2);
end