function [ed,P] = getblr(blrSpan, ed, bOffIdxAllCat, bOffIdxIntCat, sr, P)
% blink-locked time course of pre & post
% blrSpan = [-6 6];
% 2018.09.10 12secs from preparatory phase for the draft figure
% 2018.10.05 16secs. [-4 12] from preparatory phase to estimate BPR confound pattern at figure 3.
P.blrSpan   = blrSpan; % +- sec around an blink (sec) % 

[nRuns, nSbjs] = size(bOffIdxAllCat);

P.blrNsamples = diff(blrSpan)*sr+1;
P.bEvents = 2; % the kinds of blink (all/intact)
P.bEventNames = {'all','intact'};
ed.blr.bOff   = cell(P.bEvents,nSbjs);
ed.blr.tVec = blrSpan(1):1/sr:blrSpan(2);
BPRstartIdx = -blrSpan(1)*sr+1;
ed.tVecBPR = ed.blr.tVec(BPRstartIdx:BPRstartIdx+P.wSize1*sr); % time from blink offset to BPR range (P.wSize1)

for k = 1:P.bEvents
for sbj = 1:nSbjs
    % initialize
    if k == 1  % all blinks
        nBlink_tmp = numel(cat(1,bOffIdxAllCat{:,sbj}));
    elseif k == 2  % intact blinks
        nBlink_tmp = numel(cat(1,bOffIdxIntCat{:,sbj}));
    end
    ed.blr.bOff{k,sbj}     = NaN(P.blrNsamples,nBlink_tmp);
    count1 = 0; % counters for indexing    
    for r = 1:nRuns % session dimension was collapsed
        if isfield(ed,'pdCat')
            y = ed.pdCat{r,sbj};
        elseif isfield(ed,'pdcat')
            y = ed.pdcat{r,sbj};
        else
            error('no pupil time course field in the eye data structure')
        end            
        if ~isempty(y)                        
            lenData = length(y);
        
        if k == 1 % all blinks
            nBlinkRun_tmp = numel(bOffIdxAllCat{r,sbj});
        elseif k == 2 % intact blinks
            nBlinkRun_tmp = numel(bOffIdxIntCat{r,sbj});
        end
        for bl = 1:nBlinkRun_tmp                       
            count1 = count1+1;
            if k == 1 % all blinks
                eventIdx  = bOffIdxAllCat{r,sbj}(bl);
            elseif k == 2 % intact blinks
                eventIdx  = bOffIdxIntCat{r,sbj}(bl);
            end
            idxSpanL  = eventIdx+blrSpan(1)*sr;
            idxSpanR  = eventIdx+blrSpan(2)*sr;
            ed.blr.bOff{k,sbj}(:,count1)   = getelr(y,idxSpanL,idxSpanR);        
        end

        else
            try
                fprintf('%s run%d eyedata was lost\n',P.sbjNames{sbj},r)
            catch % if P.sbjNames does not exist
                fprintf('%dth run%d eyedata was lost\n',sbj,r)
            end
        end
    end
    [ed.blr.bOffAvg(:,k,sbj), ed.blr.bOffSem(:,k,sbj), ed.blr.bOffStd(:,k,sbj)]       = irfnanavgsem(ed.blr.bOff{k,sbj});    
    
    fprintf('- %dth subject processing done -\n',sbj)
end
end
% mean, std and sem across subject mean
[ed.blr.bOffAvg2, ed.blr.bOffSem2, ed.blr.bOffStd2]       = irfnanavgsem(ed.blr.bOffAvg);