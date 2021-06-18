function [BPRfree] = getBPRfree(ts, bOffIdxAll, bOnIdxAll, sr, critLen)
% This function preprocess the pupil time course data to get BPRfree periods

% It was improved form of getbf.m ('/Users/kyungyoo/Dropbox/KyungYoo/Experiments/PBPC/draftPlotCodes/sub/getbf.m')
% It does not resample anymore. Do resample outside this function.
% Last update: 2018.03.14 Kyung Yoo

% -- extract BPR-free periods who are longer than the window size (= 3 sec) --
% I will call the period BPRfree
if iscell(ts), nRun = numel(ts); else, nRun = size(ts,2); end

BPRfree.IdxAll = cell(nRun, 1); BPRfree.Idx = cell(nRun, 1);
tmpts     = cell(nRun, 1); % initialize

c = zeros(nRun,1);
for r = 1:nRun
    if iscell(ts), y = ts{r}; else, y = ts(:,r); end
    % indices in original sampling rate
    BPRfree.IdxAll{r}(:,1) = [1; bOffIdxAll{r}+critLen*sr]; % onset
    BPRfree.IdxAll{r}(:,2) = [bOnIdxAll{r}; length(y)]; % offset
    lengTmp = BPRfree.IdxAll{r}(:,2) - BPRfree.IdxAll{r}(:,1);
    idx = lengTmp > (critLen*sr+1); % The critetrion for SF period
    BPRfree.Idx{r}(:,1) = BPRfree.IdxAll{r}(idx,1); BPRfree.Idx{r}(:,2) = BPRfree.IdxAll{r}(idx,2);
    tmpts{r} = cell(size(BPRfree.Idx{r},1), 1); % initialize
    for b = 1:size(BPRfree.Idx{r},1)
        idxStart = BPRfree.Idx{r}(b,1)+1;
        idxEnd   = BPRfree.Idx{r}(b,2)-1;
        tmpts{r}{b} = y(idxStart:idxEnd); % get time course
    end
    c(r) = size(BPRfree.Idx{r},1);
end
% merge all the BPRfree spans across runs
BPRfree.Span = cell(sum(c),1); % nested cell because all BPRfree irf has different size...
for r = 1:nRun
    for b = 1:c(r)
        BPRfree.Span{b+sum(c(1:(r-1)))} = tmpts{r}{b};
    end
end

% -- segmentation: segment the BPRfree irf into critLen(= 3 sec)-long chunks --
% it is not IRF strictly speaking, but I will use the term.
BPRfree.tVec = 0:1/sr:(critLen-1/sr);
BPRfree.irf = []; BPRfree.slope = [];

for b = 1:size(BPRfree.Span,1)
    for seg = 1:floor((size(BPRfree.Span{b},1))/(critLen*sr+1))
        BPRfree.irf   = [BPRfree.irf BPRfree.Span{b}(((end-critLen*sr):end)-(seg-1)*critLen*sr)];
        BPRfree.slope = (BPRfree.irf(end,:)-BPRfree.irf(1,:))/(size(BPRfree.irf,1)/sr); % slope
    end
end


end