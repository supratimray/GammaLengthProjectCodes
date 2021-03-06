function burstLengthS = getBurstLengthFeingold(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,filterOrder)

if ~exist('displayFlag','var');         displayFlag=1;                  end
if ~exist('stimulusPeriodS','var');     stimulusPeriodS=[0.5 1.5];      end
if ~exist('baselinePeriodS','var');     baselinePeriodS=[-1 0];         end
if ~exist('burstFreqRangeHz','var');    burstFreqRangeHz=[40 60];       end
if ~exist('filterOrder','var');         filterOrder=4;                  end

numTrials=size(analogData,1);

smoothedPowerBpfSignal=zeros(size(analogData));
bpfSignal=zeros(size(analogData));
for i=1:numTrials
    [smoothedPowerBpfSignal(i,:),bpfSignal(i,:)] = getBPFPowerFeingold(analogData(i,:),timeVals,burstFreqRangeHz,filterOrder);
end

mBL = getMedianBaseline(smoothedPowerBpfSignal,timeVals,baselinePeriodS);
if isempty(thresholdFactor)
    mST = getMedianBaseline(smoothedPowerBpfSignal,timeVals,stimulusPeriodS);
    thresholdFactor = mean(mST./mBL);
    disp(['Using threshold factor of: ' num2str(thresholdFactor)]);
end

threshold=thresholdFactor*mBL;

burstLengthS=cell(1,numTrials);
for i=1:numTrials
    if displayFlag
        clf;
        disp(['Trial : ' num2str(i) ' of ' num2str(numTrials)]);
        subplot(211);plot(timeVals,bpfSignal(i,:));xlim([stimulusPeriodS(1)-0.5 stimulusPeriodS(2)+0.5]);
    end
    burstLengthS{i} = getBurstLengthFeingoldSingleTrial(smoothedPowerBpfSignal(i,:),timeVals,threshold,stimulusPeriodS,displayFlag);
end
end

function burstLengthS = getBurstLengthFeingoldSingleTrial(bandPassPowerSingleTrial,timeVals,threshold,stimulusPeriodS,displayFlag)
stPos = intersect(find(timeVals>=stimulusPeriodS(1)),find(timeVals<stimulusPeriodS(2)));
st=timeVals(1,stPos);
stBandPassPowerSingleTrial=bandPassPowerSingleTrial(1,stPos);

seedPosList=find((stBandPassPowerSingleTrial/threshold)>1); % Seed point is any time point which exceeds the threshold
burstStartPosList = [];
burstEndPosList = [];
timePos=1;
timeLength=length(st);
while (timePos < timeLength) % As long as the time position is less than the search interval
    tPos = seedPosList(find(seedPosList>timePos,1)); % Find the position of the next seed
    minPos = find((stBandPassPowerSingleTrial(1:tPos)/(threshold/2))<1,1,'last');
    if ~isempty(tPos)
        maxPos = find((stBandPassPowerSingleTrial(tPos+1:timeLength)/(threshold/2))<1,1,'first');
        
        if isempty(minPos)
            tMinPos=1;
        else
            tMinPos=minPos;
        end
        if isempty(maxPos)
            tMaxPos=timeLength;
        else
            tMaxPos=tPos+maxPos;
        end
        timePos=tMaxPos;
        burstStartPosList=cat(1,burstStartPosList,stPos(tMinPos));
        burstEndPosList=cat(1,burstEndPosList,stPos(tMaxPos));
    else
        timePos=timeLength;
    end
end
burstLengthS = timeVals(burstEndPosList)-timeVals(burstStartPosList);

if displayFlag==1
    subplot(212);plot(timeVals,bandPassPowerSingleTrial);hold on;
    line([timeVals(1),timeVals(length(timeVals))],[threshold,threshold],'color','r');
    line([timeVals(1),timeVals(length(timeVals))],[threshold/2,threshold/2],'color','g');
    for j=1:length(burstStartPosList)
        line([timeVals(burstStartPosList(j)),timeVals(burstEndPosList(j))],[(threshold/2),(threshold/2)],'lineWidth',3,'color','k');
        plot(timeVals(burstStartPosList(j)),(threshold/2),'o','linewidth',3);
    end
    xlim([stimulusPeriodS(1)-0.5 stimulusPeriodS(2)+0.5]);
    pause;
    clf;
end
end
function mBL = getMedianBaseline(smoothedPowerBpfSignal,timeVals,baselinePeriodS)
blPos = intersect(find(timeVals>=baselinePeriodS(1)),find(timeVals<baselinePeriodS(2)));
blData = smoothedPowerBpfSignal(:,blPos);
mBL = median(blData(:));
end