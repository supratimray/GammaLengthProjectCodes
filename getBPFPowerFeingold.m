function [smoothedPowerBpfSignal,bpfSignal] = getBPFPowerFeingold(signal,timeVals,gammaFreqRangeHz,filtOrder)
%BPF the signal
Fs=1/(timeVals(2)-timeVals(1));
normBand=gammaFreqRangeHz/(Fs/2);
[b,a]=butter(filtOrder,normBand,'bandpass');
bpfSignal=filtfilt(b,a,signal);
powerBpfSignal=bpfSignal.^2;
hannLen=round((2/gammaFreqRangeHz(1))*Fs);
h=hanning(hannLen);
h=h/sqrt(sum(h.^2));
smoothedPowerBpfSignalTemp=conv(powerBpfSignal,h,'same'); % Smoothed Power in the Gamma range
smoothedPowerBpfSignal=smoothedPowerBpfSignalTemp*(max(powerBpfSignal)/max(smoothedPowerBpfSignalTemp)); 
end