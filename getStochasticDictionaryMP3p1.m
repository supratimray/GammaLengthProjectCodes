% This program uses MP version 3.1 by Piotr Durka's group. This is the
% original stochastic dictionary code used in their 2001 paper.

function [gaborInfo,header] = getStochasticDictionaryMP3p1(data,timeVals,maxIteration,adaptiveDictionaryParam,dictionarySize)

if ~exist('maxIteration','var');           maxIteration=50;             end
if ~exist('adaptiveDictionaryParam','var');adaptiveDictionaryParam=0.9; end
if ~exist('dictionarySize','var');      dictionarySize=[];              end

Fs=round(1/(timeVals(2)-timeVals(1)));          % Sampling Frequency
numTrials = size(data,1);
sigLen = size(data,2);

%%%%%%%%%%%%%%%%%%%%%_________MP Parameters _________%%%%%%%%%%%%%%%%%%%
recAc=100; % Pecentage of Maximum reconstructed Energy (Value should be less than 100)

%%%%%%%%%%%%%%%%%%_______ Run MP on single trials ________%%%%%%%%%%%%%%
gaborInfo = zeros(numTrials,maxIteration,7);
header = zeros(numTrials,8);

for i=1:numTrials
    disp(['Trial ' num2str(i) ' of ' num2str(numTrials)]);
    
    % Save data as ASCII files
    sig=data(i,:)'; %#ok<NASGU>
    save('sig.txt','sig','-ascii');
    clear sig;
    
    % Write Command File
    fp=fopen('commands.txt','w');
    if ~isempty(dictionarySize)
        fprintf(fp,['reinit -R ' num2str(dictionarySize) ' \n']);
    end
    fprintf(fp,['set -O ' num2str(sigLen) ' -M ' num2str(maxIteration) ' -E ' num2str(recAc) ' -F ' num2str(Fs) ' -D ' num2str(adaptiveDictionaryParam) ' \n']);
    fprintf(fp,'loadasc -O sig.txt\n');
    fprintf(fp,'mp\n'); % Run MP
    fprintf(fp,'save -S mpresults\n');
    fprintf(fp,'exit');        % Exit shell
    fclose(fp);
    
    % Run MP
    system([which('mp31.exe') ' < commands.txt'],'-echo');
    
    % Read Data
    [gaborInfo(i,:,:), header(i,:)]=readbook('mpresults',0);
    
    % Delete tmp files
    delete sig.txt;                         % Deletes the signal file
    delete commands.txt;                    % Deletes the commamnds file created
    delete mpresults;                       % Deletes mpresults
end
end