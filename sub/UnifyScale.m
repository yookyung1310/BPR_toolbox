function [xUniLim, yUniLim] = UnifyScale(h, varargin)
% UnifyScale(h) unifies Y scale of all axes of the figure handle h.
% UnifyScale(h, tgt) specifies target axes. tgt is a vector 

    % set target axes
    tgtAxes = findall(h,'type','axes');
    if nargin ~= 1
        tgtAxes = tgtAxes(varargin{1});
    end

    % Gather XY range of target axes
    for i = 1:numel(tgtAxes)
        xLims(i,:) = [min(tgtAxes(i).XLim) max(tgtAxes(i).XLim)];
        yLims(i,:) = [min(tgtAxes(i).YLim) max(tgtAxes(i).YLim)];
    end

    % Decide the unified limit
    % 1) just use min max of all axes
    xUniLim(1) = min(xLims(:,1)); xUniLim(2) = max(xLims(:,2));
    yUniLim(1) = min(yLims(:,1)); yUniLim(2) = max(yLims(:,2));

    % Apply it to the figure
    for i = 1:numel(tgtAxes)
        xlim(tgtAxes(i),xUniLim);
        ylim(tgtAxes(i),yUniLim);
    end
		
end