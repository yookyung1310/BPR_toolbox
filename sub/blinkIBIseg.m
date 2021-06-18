function [strictBidx, idx, segName, idx_dbl, dblName] = blinkIBIseg(bIdx, CritIdx, endIdx, varargin)

% get indices of grouped blinks on the pre,post IBI space.
% Optionally, it returns double blinks, the blinks that distance between them is close and the distance to other
% blinks are far (b---B-b---b)
segName = {'b-B---b','b-B-b','b---B---b','b---B-b'}; % the 4 patterns
dblName = {'b---B-b---b','b---b-B---b'}; % pre and post of the double blinks.

if ~isempty(bIdx)
    bOnIdx = bIdx(:,1); bOffIdx = bIdx(:,2);
    preCritIdx = CritIdx(1); postCritIdx = CritIdx(2);    
            
    blinkDiff = bOnIdx(2:end) - bOffIdx(1:end-1);
    preIBI  = [bOnIdx(1); blinkDiff];
    postIBI = [blinkDiff; endIdx-bOffIdx(end)];
    post2IBI = [postIBI(2:end); NaN]; % 2 blinks forward
    
    idx = false(numel(preIBI),4);
    idx(:,1) = preIBI<=preCritIdx & postIBI>postCritIdx;
    idx(:,2) = preIBI<=preCritIdx & postIBI<=postCritIdx;
    idx(:,3) = preIBI>preCritIdx  & postIBI>postCritIdx;
    idx(:,4) = preIBI>preCritIdx  & postIBI<=postCritIdx;
    
    strictBidx = idx(:,3); % the 'b---B---b' case is useful for intact blink selection so assign a variable.
    
    if nargin > 3
        % get indices of double blinks
        dbl.preCritIdx   = varargin{1}(1); % pre-IBIidx, post-IBIidx, 2nd post-IBIidx
        dbl.postCritIdx  = varargin{1}(2);
        dbl.post2CritIdx = varargin{1}(3);                
        
        idx_dbl = false(numel(preIBI),2);
        idx_dbl(:,1)  = preIBI>dbl.preCritIdx & postIBI<dbl.postCritIdx & post2IBI>dbl.post2CritIdx;
        idx_dbl(:,2) = [false; idx_dbl(1:end-1,1)];
    end
%     %check plot
%     figure(); hold on;
%     for i = 1:4
%        plot(preIBI(idx(:,i)), postIBI(idx(:,i)), '.');       
%        pause();
%     end
%     vline(preCritIdx); hline(postCritIdx);
%     xlabel('pre IBI'); ylabel('post IBI');
%     
%     figure(); hold on;
%     plot(1, preIBI(idx_dbl(:,1)), '.')
%     plot(2, postIBI(idx_dbl(:,1)), '.')
%     plot(3, post2IBI(idx_dbl(:,1)), '.')
else
    warning('blink index is empty')
    strictBidx = [];
    idx = [];
    idx_dbl = [];
end