% This function computes the instantaneous power at gamma band using
% Hilbert transform. It uses a Butterworth filter to bandpass filter the
% signal in the gamma range

% Input - signal              : Signal in one row
%         timeVals            : Time Values in a row (seconds)
%         gammaFreqRangeHz    : Frequency range for which the Wavelet
%                               transform should be computed
%         filtOrder           : Order of the Butterworth filter

% Output - hilbertPower       : Instantaneous power calculated using
%                               Hilbert transform
%          bpfSignal          : Signal bandpass filtered in gamma band
%                               using Butterworth filter

function [hilbertPower,bpfSignal] = getHilbertPower(signal,timeVals,gammaFreqRangeHz,filtOrder)
%BPF the signal
Fs=1/(timeVals(2)-timeVals(1));
normBand=gammaFreqRangeHz/(Fs/2);
[b,a]=butter(filtOrder,normBand,'bandpass');
bpfSignal=filtfilt(b,a,signal);

%Find Inst. power using HT
hilbSignal=hilbert(bpfSignal);
hilbertPower=abs(hilbSignal).^2;
end