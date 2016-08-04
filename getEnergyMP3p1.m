function [meanE,freqVals]=getEnergyMP3p1(gaborInfo,header,timeVals)
numTrials = size(gaborInfo,1);
Fs=round(1/(timeVals(2)-timeVals(1)));
N=length(timeVals);
freqVals = 0:Fs/N:Fs/2-Fs/N;
sumE = zeros(N/2,N);
for i=1:numTrials
    disp(['Computing Energy for trial ' num2str(i) ' of ' num2str(numTrials)]);
    E = mp2tf(squeeze(gaborInfo(i,:,:)),header(i,:));
    sumE = sumE+E;
end
meanE = sumE/numTrials;
end