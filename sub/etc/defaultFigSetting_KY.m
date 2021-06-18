function defaultFigSetting_KY(varargin)
set(groot,'defaultAxesFontSize',20,'defaultAxesFontweight','demi')
set(groot,'defaultLineLineWidth',2)
set(groot, 'DefaultFigurePosition', [-1136 338 560 420]);
set(groot,'DefaultFigureWindowStyle','docked')

if nargin > 0
    set(groot, 'DefaultFigurePosition', varargin{1});
% set(groot, 'DefaultFigurePosition', [0 0 2560 1440]);

% % if you want to suppress figure pop up...
% set(groot, 'DefaultFigureVisible', 'off')

% 2016/11/18 added

end

end

