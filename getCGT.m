% This function computes the CGT (Continuous Gabor Transform) of a signal
% Input - sig       : Signal in one row
%         timeVals  : Time Values in a row (seconds)
%         fRange    : Frequency range for which the CGT should be computed
%         fRes      : Frequency resolution
%         sd        : Standard deviation of the Gaussian window in seconds.

% Output - cgt      : CGT matrix (time in x direction and frequency in y)
%          freqVals : frequency values at which CGT is computed 

function [cgt,freqVals]=getCGT(sig,timeVals,fRange,fRes,sd)

if length(fRange)==1
    fRange = [0 fRange];
end

freqVals = fRange(1):fRes:fRange(2); % compute CGT at these frequencies
numFreqVals = length(freqVals);
numTimeVals = length(timeVals);
cgt = zeros(numFreqVals,numTimeVals);

for i=1:numFreqVals
    h=(1/(sd*sqrt(2*pi)))*((exp(-((timeVals-timeVals(numTimeVals/2))/(sqrt(2)*(sd))).^2 )).*(exp(1j*2*pi*freqVals(i)*timeVals)));
    cgt(i,:)=conv(sig,h,'same');
end
end