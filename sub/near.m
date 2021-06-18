function [raw_mod,idx] = near(raw,dest,varargin)
% find the nearest point.
% raw (single vector) to modify any nearest point of dest (single vector)
% NOTE1: Second inputs should be monotonically increasing.
% NOTE2: If first input is out of range of second input, it just ignore the out of
% range points.
% Made by KY
if nargin > 2
    wflag = varargin{1};
end
[~, tmpIdx] = histc(raw,dest);

if any(tmpIdx == 0)
    if ~exist('wflag','var')        
        warning('Be aware that the first input is out of range of the second input.')
    end
    raw     = raw(tmpIdx~=0);
    tmpIdx  = tmpIdx(tmpIdx~=0);    
end

raw_mod = zeros(size(raw));
idx     = zeros(size(raw));

for i = 1:numel(raw)    
    if tmpIdx(i) == numel(dest) % when the value is the exact last value
        raw_mod(i) = dest(tmpIdx(i));
        idx(i) = tmpIdx(i);
    elseif abs(raw(i) - dest(tmpIdx(i))) <= abs(raw(i) - dest(tmpIdx(i)+1))
        raw_mod(i) = dest(tmpIdx(i));
        idx(i) = tmpIdx(i);
    else
        raw_mod(i) = dest(tmpIdx(i)+1);
        idx(i) = tmpIdx(i)+1;
    end
end

