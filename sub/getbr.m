function [br, bPulse] = getbr(bOffIdxAllCat, pdCat, sr, varargin)
wSize_br = 0.2; % sliding window size (sec)
if nargin > 3
    wSize_br = varargin{1};
end

bPulse = cell(size(bOffIdxAllCat));
br     = cell(size(bOffIdxAllCat));
for i = 1:numel(bOffIdxAllCat)
    if ~isempty(pdCat{i})
        bPulse{i} = false(size(pdCat{i}));
        bPulse{i}(bOffIdxAllCat{i}) = true;
        br{i} = movmean(bPulse{i}, sr*wSize_br)*sr; % average window = 200ms
    end
end

