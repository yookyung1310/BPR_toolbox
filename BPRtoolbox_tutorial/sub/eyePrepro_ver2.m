function [pupil, bOnIdxAll, bOnIdxInt, bOffIdxAll, bOffIdxInt, bOffIdxArtAll, flagConta, isnBlinkOdd, fNum] = ...
    eyePrepro_ver2(pupil, sr, varargin)
% pupil, gazeX and gazeY should be vectors
%   pupil: pupil diameter (mm) or pupil size (au, from eyelink eyetracker) (vector)
%   gazeX: horizontal  (vector)
%   gazeY: indices of all blinks' onset in the time course (vector)

% It is highly recommend user to convert the pupil unit au to mm, because this code preprocess eyelid occlusion
% according to pupil diameter. It preprocess with rough value in au case but that might be inaccurate.
% conservative criteria for blink detection, and more liberal criteria for artifact defining.

%%% Kyoungwhan Choe's code modified by KY
%%% last update: 2019.05.30
%   Copyright 2019 Kyung Yoo
%   Contatct: yookyung1310@gmail.com

firstaxesinput = nargin-numel(varargin);
opts = parseinput(varargin, 'eyePrepro_ver2', firstaxesinput);


% blink detection: use pupil diameter value + high passed (for partial blinks)
fc = 10; % cutoff frequency (Hz)
order = 3;
[b,a] = butter(order, fc/(sr/2), 'high');
pd_hp = filtfilt(b, a, pupil-mean(pupil));

% NOTE: blinkPt = (pupil(:) < 1) | vertcat((abs(sm_pdDeriv) > 2),0); was changed to
%       blinkPt = (pupil(:) < 1) | (abs(pd_hp) > 0.25); because some BPR was larger than -2mm/s. 
%       
if opts.ismm
    blinkPt = (pupil(:) < 1) | (abs(pd_hp) > 0.25);
    derivCrit = 0.1; % criterion to detect artifact
else
%     if there is no configuration information 
    blinkPt = (pupil(:) < 1000) | (abs(pd_hp) > 250);
    derivCrit = 100; % criterion to detect artifact
end

isblink = xor(blinkPt,[blinkPt(2:end); 0]);
% store the point processed as blink from RAW
bIdx_noFixSep = find(isblink == 1);

% check the blink points are even number
if mod(length(bIdx_noFixSep),2)==1 % it means that the subject closed his eye at the start of recording
    isnBlinkOdd = true;
%     bIdx_noFixSep = bIdx_noFixSep(2:end); % ignore the first blink 
    bIdx_noFixSep = [1; bIdx_noFixSep]; % consider the start as blink onset
end

bOnIdxAll  = bIdx_noFixSep(1:2:end);
bOffIdxAll = bIdx_noFixSep(2:2:end);
[bOnIdxAll, bOffIdxAll] = FixSepBlinks(bOnIdxAll, bOffIdxAll, opts.sepCrit*sr);

bOffIdxArtAll = NaN(size(bOffIdxAll));
lenData   = length(pupil);
flagConta = false(lenData,1);

% to define artifact and eliminate it, add criteria for blink ARTIFACT offset
boxSize = 0.25; % (sec)
sm_pdDeriv = movmean(diff(pupil)/(1/sr), sr*boxSize);
for i = 1:length(bOffIdxAll)    
    lIdx = bOnIdxAll(i)-opts.conta(1)*sr; 
%     rIdx = bOffIdxAll(i)+opts.conta*sr;
    % detect artifact from offset with smoothed pupil derivative
    artDetectSpan = 0.5*sr;
    artIdx = artDetectSpan+bOffIdxAll(i);
    % if it is out of data length
    if artIdx > (lenData-1) % it was diff so the length is minus 1
       artIdx = (lenData-1); 
    end
    % Find the point where smoothed pupil derivative gets smaller than 0.1mm/s within the range_art, 
    % and regard it as an artifact offset.
    artDur = find( sm_pdDeriv(bOffIdxAll(i):(artIdx-1))>=derivCrit & sm_pdDeriv((bOffIdxAll(i)+1):artIdx)<derivCrit, 1 );    
    if isempty(artDur) % if it could not find artifact endpoint 
        artDur = opts.conta(3)*sr;
    end
    if artDur < opts.conta(2)*sr % if detected artifact endpoint is too near from blink offset, there was some artifact still. so endow at least 200ms
        artDur = opts.conta(2)*sr;
    end
    rIdx = bOffIdxAll(i)+artDur;
    % mark +-conta_ms (ms) around the lost-pupil periods
    if lIdx < 1
        flagConta(1:rIdx) = true;
    elseif rIdx > lenData
        flagConta(lIdx:end) = true;
    else
        flagConta(lIdx:rIdx) = true;
    end
    bOffIdxArtAll(i) = rIdx;
end
pupil(flagConta) = nan; 
% gazeX(flagConta) = nan; 
% gazeY(flagConta) = nan;

% linear interpolation
if ~all(isnan(pupil))
    pupil = FixGaps(pupil);
else
    flagConta(1:end) = true; % whole recording is invalid.
end

% get intact blink 
if ~isempty(opts.isolatedCrit) && ~isempty(bOffIdxAll)
    wSize =  opts.isolatedCrit; % window size to define intact blink
    blinkDiff  = bOnIdxAll(2:end) - bOffIdxAll(1:end-1);
    isIntactB  = [bOnIdxAll(1); blinkDiff] > wSize*sr  &  [blinkDiff; lenData-bOffIdxAll(end)] > wSize*sr;
    bOffIdxInt = bOffIdxAll(isIntactB);
    bOnIdxInt  = bOnIdxAll(isIntactB);
else
    bOffIdxInt = [];
    bOnIdxInt  = [];
end


if ~opts.suppressFig
    fhandle(1) = figure(); title('pupil area (au)'); hold on
    plot((0:(length(pupil)-1))/sr,pupil);
    xlabel('time (s)')
    fhandle(2) = figure(); title('gaze horizontal (pixel)'); hold on
    plot((0:(length(gazeX)-1))/sr,gazeX);
    xlabel('time (s)')
    fhandle(3) = figure(); title('gaze vertical (pixel)'); hold on
    plot((0:(length(gazeY)-1))/sr,gazeY);
    xlabel('time (s)')
    fNum = get(fhandle,'Number');
else
    fNum = [];
end
 

end
% if count == Nframes
%     fprintf('run %d has equal number of flips\n',irun)
% end


function [opts] = parseinput(input, funcname, inputoffset)
% default values
opts = struct('sepCrit',0.2, 'conta',[0.15 0.2 0.5], 'isolatedCrit',[], 'ismm',false, 'suppressFig',false);
% opts.Color = {'k','b','m','r'};

names = {'sepCrit','conta','ismm','suppressFig','isolatedCrit'};
inputlen = length(input);
for i = 1:2:inputlen
    
    name = validatestring(input{i},names);
    
    value = input{i+1};
    switch name  
        case 'sepCrit'
            validateattributes(value,{'numeric','logical'},{'real', ...
                'scalar','nonnegative','finite'}, funcname, 'sepCrit', i+1+inputoffset)
            opts.sepCrit = value;
        case 'conta'
            validateattributes(value,{'numeric','logical'},{'real', ...
                'vector','nonnegative','finite'}, funcname, 'conta', i+1+inputoffset)
            opts.conta = value;
        case 'isolatedCrit'
            validateattributes(value,{'numeric','logical'},{'real', ...
                'scalar','nonnegative','finite'}, funcname, 'isolatedCrit', i+1+inputoffset)
            opts.isolatedCrit = value;
        case 'ismm'
            validateattributes(value,{'numeric','logical'},{'scalar', ...
                'integer'}, funcname, 'ismm', i+1+inputoffset)
            opts.ismm = value;
        case 'suppressFig'
            validateattributes(value,{'numeric','logical'},{'scalar', ...
                'integer'}, funcname, 'suppressFig', i+1+inputoffset)            
            opts.suppressFig = value;
    end
end
end


