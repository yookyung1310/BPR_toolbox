function [irf_rs] = resampleirf2(irf, sr, sr_rs, varargin)
%resamplets Resample time series properly by padding
%   irf: irf time series matrix. row: time, column: events
% only " irf_rs = cell(numel(irf),1);" has been changed to "irf_rs = cell(size(irf));"   2018.03.22 KY

padCoeff = 100; % default
if nargin > 3
    padCoeff = varargin{1};
end

% the beginning and end of the data is not zero, so padding
if iscell(irf)
    irf_rs = cell(size(irf));
    for r = 1:numel(irf)
        if ~isempty(irf{r})
            irfpad = [repmat(irf{r}(1,:), sr/sr_rs*padCoeff, 1); irf{r}; repmat(irf{r}(end,:), sr/sr_rs*padCoeff, 1)];
            irfpad_rs = resample(irfpad, sr_rs, sr);
            % now get rid of the pad
            irf_rs{r} = irfpad_rs((1+padCoeff):end-padCoeff,:);
        else % no irf case
            % search all runs to find out irf length
            for i = 1:numel(irf)
                if numel(irf{i})>0
                    irfLen = size(irf{i},1);
                    break;
                end
            end
            % empty vector
            irf_rs{r} = NaN(irfLen,0);            
        end
    end
else % matrix
    if ~isempty(irf)        
        irfpad = [repmat(irf(1,:), sr/sr_rs*padCoeff, 1); irf; repmat(irf(end,:), sr/sr_rs*padCoeff, 1)];
        irfpad_rs = resample(irfpad, sr_rs, sr);
        % now get rid of the pad
        irf_rs = irfpad_rs((1+padCoeff):end-padCoeff,:);
    else % no irf case
        disp('no blink in the condition...')
    end
end

end

