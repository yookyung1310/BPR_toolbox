function [pupil, gazeX, gazeY, bOnIdxAll, bOffIdxAll, bOffIdxInt, flagConta, isnBlinkOdd, fNum] = ...
    eyePrepro(pupil, conta_sec, sr, diffCrit, varargin)
% pupil, gazeX and gazeY should be vectors

% blink detection: use pa value + pa difference 
blinkPt = (pupil < 100) | vertcat((abs(diff(pupil)) > 50),0); 

isblink = xor(blinkPt,[blinkPt(2:end); 0]);
% store the point processed as blink from RAW
bIdx_noFixSep = find(isblink == 1);

% check the blink points are even number
if mod(length(bIdx_noFixSep),2)==1 % it means that the subject closed his eye at the start of recording
    isnBlinkOdd = true;
    bIdx_noFixSep = bIdx_noFixSep(2:end);
end

bOnIdxAll  = bIdx_noFixSep(1:2:end);
bOffIdxAll = bIdx_noFixSep(2:2:end);
[bOnIdxAll, bOffIdxAll] = FixSepBlinks(bOnIdxAll, bOffIdxAll, diffCrit*sr);

lenData = length(pupil);
flagConta = false(lenData,1);
for i = 1:length(bOffIdxAll)
    lIdx = bOnIdxAll(i)-conta_sec*sr; rIdx = bOffIdxAll(i)+conta_sec*sr;
    % mark +-conta_ms (ms) around the lost-pupil periods
    if lIdx < 1
        flagConta(1:rIdx) = true;
    elseif rIdx > lenData
        flagConta(lIdx:end) = true;
    else
        flagConta(lIdx:rIdx) = true;
    end
end
pupil(flagConta) = nan; gazeX(flagConta) = nan; gazeY(flagConta) = nan;

% linear interpolation
pupil = FixGaps(pupil);

% get intact blink 
if nargin > 4    
    wSize = varargin{1}; % window size to define intact blink
    blinkDiff = bOnIdxAll(2:end) - bOffIdxAll(1:end-1);
    isIntactB = [bOnIdxAll(1); blinkDiff] > wSize*sr  &  [blinkDiff; lenData-bOffIdxAll(end)] > wSize*sr;    
    bOffIdxInt = bOffIdxAll(isIntactB);
end


if nargin > 5
    if strcmp(varargin{2},'draw')    
    fhandle(1) = figure(); title('pupil area (au)'); hold on
    plot((0:(length(pupil)-1))/sr,pupil);
    xlabel('time (s)')
    fhandle(2) = figure(); title('gaze horizontal'); hold on
    plot((0:(length(gazeX)-1))/sr,gazeX);
    xlabel('time (s)')
    fhandle(3) = figure(); title('gaze vertical'); hold on
    plot((0:(length(gazeY)-1))/sr,gazeY);
    xlabel('time (s)')
    if PreRaw=='P'
        figure(fhandle(1)); ylabel('pupil area (au)')
        figure(fhandle(2)); ylabel('gaze x (pixel)')
        figure(fhandle(3)); ylabel('gaze y (pixel)')        
    else
        figure(fhandle(1)); ylabel('pupil area (mean %)')
    end
    fNum = get(fhandle,'Number');
    else
        fNum = [];
    end
end
 

end
% if count == Nframes
%     fprintf('run %d has equal number of flips\n',irun)
% end




