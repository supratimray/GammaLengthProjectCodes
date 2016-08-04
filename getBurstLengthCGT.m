function [burstLengthS,burstFreqList] = getBurstLengthCGT(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,cgtGaborSDS,cgtFreqResolutionHz,searchRangeFreqHz,phaseThreshold,useMaxPowerSeedFlag)

if ~exist('thresholdFactor','var');     thresholdFactor=[];             end
if ~exist('displayFlag','var');         displayFlag=1;                  end
if ~exist('stimulusPeriodS','var');     stimulusPeriodS=[0.5 1.5];      end
if ~exist('baselinePeriodS','var');     baselinePeriodS=[-1 0];         end
if ~exist('burstFreqRangeHz','var');    burstFreqRangeHz=[40 60];       end
if ~exist('cgtGaborSDS','var');         cgtGaborSDS=0.025;              end
if ~exist('cgtFreqResolutionHz','var'); cgtFreqResolutionHz=2.5;        end
if ~exist('searchRangeFreqHz','var');   searchRangeFreqHz=2.5;          end
if ~exist('phaseThreshold','var');      phaseThreshold=pi/4;            end
if ~exist('useMaxPowerSeedFlag','var'); useMaxPowerSeedFlag=1;          end

mBL = getMeanBaseline(analogData,timeVals,baselinePeriodS,burstFreqRangeHz,cgtFreqResolutionHz,cgtGaborSDS);
if isempty(thresholdFactor)
    mST = getMeanBaseline(analogData,timeVals,stimulusPeriodS,burstFreqRangeHz,cgtFreqResolutionHz,cgtGaborSDS);
    thresholdFactor = mean(mST./mBL);
    disp(['Using threshold factor of: ' num2str(thresholdFactor)]);
end

threshold=thresholdFactor*mBL;

numTrials = size(analogData,1);
burstLengthS=cell(1,numTrials);
burstFreqList=cell(1,numTrials);
for i=1:numTrials
    if displayFlag
        disp(['Trial : ' num2str(i) ' of ' num2str(numTrials)]);
    end
    [burstLengthS{i},burstFreqList{i}] = getBurstLengthCGTSingleTrial(analogData(i,:),timeVals,threshold,displayFlag,cgtGaborSDS,stimulusPeriodS,burstFreqRangeHz,cgtFreqResolutionHz,searchRangeFreqHz,phaseThreshold,useMaxPowerSeedFlag,thresholdFactor);
end
end
function [burstLengthListS,burstFreqList] = getBurstLengthCGTSingleTrial(signal,timeVals,threshold,displayFlag,cgtGaborSDS,stimulusPeriodS,burstFreqRangeHz,cgtFreqResolutionHz,searchRangeFreqHz,phaseThreshold,useMaxPowerSeedFlag,thresholdFactor)

stPos = intersect(find(timeVals>=stimulusPeriodS(1)),find(timeVals<stimulusPeriodS(2)));

[cgt,freqVals]=getCGT(signal,timeVals,burstFreqRangeHz,cgtFreqResolutionHz,cgtGaborSDS);
stPower = abs(cgt(:,stPos)).^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Burst Detection %%%%%%%%%%%%%%%%%%%%%%%%%%
% Bursts are computed at "seed" values, at which power is more than the
% threshold. We assume that there is a single seed at a particular time.

% First, we find all time points at which the power at at least frequency
% bin exceeds the threshold
searchLengthFreq = size(stPower,1);
searchLengthTime = size(stPower,2);
seedFreqPos = zeros(1,searchLengthTime);

for i = 1:searchLengthTime % For each time position 
    bestFreqPos = find(stPower(:,i)==max(stPower(:,i))); % Find out the frequency at which power is highest
    if stPower(bestFreqPos,i)>threshold(bestFreqPos)
        seedFreqPos(i) = bestFreqPos;
    end
end
timePointsWithSeeds = find(seedFreqPos>0);

% Now we march through the time axis
timePos=1; burstCount=1;
searchRangeFreqPos = ceil(searchRangeFreqHz/cgtFreqResolutionHz);
seedPosList = [];
burstFreqList = [];
burstStartPosList = [];
burstEndPosList = [];

while (timePos < searchLengthTime) % As long as the time position is less than the search interval    
    tPos = timePointsWithSeeds(find(timePointsWithSeeds>timePos,1)); % Find the position of the next seed
    if ~isempty(tPos)
        fPos = seedFreqPos(tPos);
        if useMaxPowerSeedFlag % from this time point, keep going forward as long as power increases
            diffPower = diff(stPower(fPos,tPos:end));
            tPosAdvance = find(diffPower<0,1);
            if ~isempty(tPosAdvance) && ((tPos+tPosAdvance)<=searchLengthTime)
                tPos=tPos+tPosAdvance;
            end
        end
        
        seedPosList = cat(2,seedPosList,tPos);
        fPosList = max(1,fPos-searchRangeFreqPos):min(searchLengthFreq,fPos+searchRangeFreqPos);
        numFPosList = length(fPosList);
        
        tmpBurstLength = zeros(1,numFPosList);
        tmpBurstMinPos = zeros(1,numFPosList);
        tmpBurstMaxPos = zeros(1,numFPosList);
        for i=1:numFPosList
            [tmpBurstMinPos(i),tmpBurstMaxPos(i)] = getBurstLengthAtSeed(cgt(fPosList(i),stPos),tPos,timeVals(stPos),freqVals(fPosList(i)),phaseThreshold);
            tmpBurstLength(i)=tmpBurstMaxPos(i)-tmpBurstMinPos(i);
        end
        
        % Get Burst details
        maxPos = find(tmpBurstLength==max(tmpBurstLength),1);
        burstFreqList=cat(2,burstFreqList,freqVals(fPosList(maxPos))); %#ok<*NASGU>
        burstStartPosList=cat(2,burstStartPosList,tmpBurstMinPos(maxPos));
        burstEndPosList=cat(2,burstEndPosList,tmpBurstMaxPos(maxPos));
        timePos = tmpBurstMaxPos(maxPos);
    else
        timePos = searchLengthTime+1; % Terminate seacrh
    end
end
burstStartListS = timeVals(stPos(burstStartPosList));
burstEndListS   = timeVals(stPos(burstEndPosList));
burstLengthListS = burstEndListS-burstStartListS;
if displayFlag
    clf;
    subplot(211)
    plot(timeVals(stPos),signal(stPos)); 
    
    subplot(212)
    pcolor(timeVals(stPos),freqVals,stPower); 
    shading interp; colormap jet;
    hold on;
    
    for i=1:length(burstFreqList)
        tmpPos = burstStartPosList(i):burstEndPosList(i);
        plot(timeVals(stPos(tmpPos)),burstFreqList(i)+zeros(1,length(tmpPos)),'k','linewidth',3);
        plot(timeVals(stPos(seedPosList(i))),burstFreqList(i),'sk');
    end
    title(thresholdFactor);
    pause;
end
end
function [tMinPos,tMaxPos]=getBurstLengthAtSeed(x,tPos,t,f0,phaseThreshold)
tLen = size(x,2);
cgtPhase=unwrap(angle(x)); % unwraped time phase of cgt
linearPhase = 2*pi*f0*t; % Linear phase shift at the same frequency
cgtPhaseCorrected = (cgtPhase-linearPhase); % Fluctuations about the linear phase

phaseDiff = abs(cgtPhaseCorrected-cgtPhaseCorrected(tPos));
minPos = find(phaseDiff(1:tPos)>phaseThreshold,1,'last');
maxPos = find(phaseDiff(tPos+1:tLen)>phaseThreshold,1,'first');

if isempty(minPos)
    tMinPos=1;                       
else
    tMinPos=minPos;
end
if isempty(maxPos)
    tMaxPos=tLen;
else
    tMaxPos=tPos+maxPos;
end
end
function meanBL = getMeanBaseline(analogData,timeVals,baselinePeriodS,burstFreqRangeHz,cgtFreqResolutionHz,cgtGaborSDS)

numTrials = size(analogData,1);
blPos = intersect(find(timeVals>=baselinePeriodS(1)),find(timeVals<baselinePeriodS(2)));

for i=1:numTrials
    cgt(i,:,:)=getCGT(analogData(i,:),timeVals,burstFreqRangeHz,cgtFreqResolutionHz,cgtGaborSDS); %#ok<AGROW>
end
blPower = abs(cgt(:,:,blPos)).^2;

meanBL = squeeze(mean(mean(blPower,3),1));
end