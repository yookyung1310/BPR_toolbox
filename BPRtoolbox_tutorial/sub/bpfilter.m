function [ts2] = bpfilter(fc, fs, order, ts, varargin)
% bandpass filtering with butterworth filter
ts2 = ts-mean(ts);
if nargin > 4
    [b,a] = butter(order, fc/(fs/2), varargin{1});
else    
    [b,a] = butter(order, fc/(fs/2));
end
ts2 = filtfilt(b,a,ts2);
ts2 = ts2 + mean(ts);
end

