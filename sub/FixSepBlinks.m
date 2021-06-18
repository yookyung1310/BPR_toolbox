function [blinkOn2, blinkOff2] = FixSepBlinks(blinkOn1, blinkOff1, diffCrit)
% Integrate separated blinks by eyetracker missing into a single blink

% 2016.11.02 
% blink offset is defined as the time when the eyetracker succeed to
% detect pupil RELIABLY. Visual insepction shows that it has some
% range that unreliably detect pupil, especially around the beggening
% and end of blinks. So we should integrate the unreliable
% detection preiod into a single segment
% NOTE: This cannot rule out short missing while eye open.

    blinkDiff = blinkOn1(2:end) - blinkOff1(1:end-1);
    if ~isempty(blinkDiff)
        % assume that blink cannot occur within specific span after blink offset
        bEndIdx   = find(blinkDiff > diffCrit);
        bStartIdx = bEndIdx+1;
        bEndIdx   = [bEndIdx; numel(blinkOff1)]; % index of the last blink
        bStartIdx = [1; bStartIdx]; % index of the first blink
        blinkOff2 = blinkOff1(bEndIdx);
        blinkOn2  = blinkOn1(bStartIdx);
    else % if the run has no blink at all...
        blinkOn2  = [];
        blinkOff2 = [];
    end