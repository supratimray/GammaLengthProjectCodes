% This function computes the scalogram (Wavelet Transform (WT)) of a signal
% Input - sig       : Signal in one row
%         timeVals  : Time Values in a row (seconds)
%         fRange    : Frequency range for which the Wavelet transform should be computed
%         fRes      : Frequency resolution

% Output - cw1      : Scalogram matrix (time in x direction and frequency in y)
%          freqVals : frequency values at which Scalogram is computed 

function [cw1,fListWavelet,tListWavelet]=getWavelet(data,timeVals,fRange,fRes)
    if length(fRange)==1
    fRange = [0 fRange];
    end
    
    Ts=timeVals(2)-timeVals(1);
    
    fListWaveletTemp = fRange(1):fRes:fRange(2); % compute WT at these frequencies
    fListWavelet = fliplr(fListWaveletTemp);
    
    % We choose the morlet wavelet
    wname = 'cmor4-1';
    cf = centfrq(wname);                         % centre frequency of morlet wavelet
    
    % We compute the scales for which the wavelet has the desired centre frequency.
    scaleList = cf./(Ts*fListWavelet);           % list of scales
    tListWavelet = timeVals;
    
    cw1 = cwt(data,scaleList,wname);             % Continuous Wavelet Transform
end
     
  
 
