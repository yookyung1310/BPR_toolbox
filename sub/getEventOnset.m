function [sampleIdx, onsetTS, eventIdx] = getEventOnset(Eyedata, fieldString, tsword)
% This function get event onset from eyetracker data from eyelink
% fieldString: 'message' or 'codestring'. Case-insensitive.

% sampleIdx: index of the tsword event in time series sample
% onsetTS  : time of the tsword event 
% eventIdx : index of the tsword event in event vectrors. (nth event in events determined by eyelink)
% NOTE: fieldString is case-insensitive, but tsword is CASE-SENSITIVE
% FEVENT message
tp = Eyedata.FSAMPLE.time; % time points, or timestamps
onsetTS = [];
% field name
evalStr = sprintf('isequal(Eyedata.FEVENT(i).%s, ''%s'')', lower(fieldString), tsword);

count = 0;
eventIdx = [];
for i=1:length(Eyedata.FEVENT)
    if eval(evalStr) 
        count=count+1;
        onsetTS(count) = Eyedata.FEVENT(i).sttime;
        eventIdx = [eventIdx i];
    end
end
if count == 0
    warning(sprintf('The timestamp word is not Eyedata.FEVENT.%s at all, so all outputs are empty.', lower(fieldString)))
    sampleIdx = [];
    return
else
    % sample index when the message occured in the time course samples
    sampleIdx = zeros(1,count);
    for i=1:count
        [~, sampleIdx(i)] = min(abs(onsetTS(i)-double(tp)));
    end
end

end
