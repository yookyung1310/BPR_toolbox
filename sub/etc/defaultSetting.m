function [path, analysisDate, opts] = defaultSetting(varargin)
% path, default variable and figure setting
% figure setting is optional
firstaxesinput = nargin-numel(varargin);
opts = parseinput(varargin, 'defaultSetting', firstaxesinput);       

if ~isfield(opts,'mainpath')
    error('mainpath is necessary')
end

path.codepath  = opts.codepath;

addpath(genpath(path.codepath));

path.mainpath  = opts.mainpath;
path.datapath  = fullfile(path.mainpath,'data');
path.mergepath = fullfile(path.mainpath,'merged');
path.mergepath_post = fullfile(path.mergepath,'AnlsProcd');

analysisDate = datestr(now,'yymmddHHMM');
path.figpath   = fullfile(path.mainpath,'fig',analysisDate);

% figure setting
if opts.isFigSet
    global defaultFontSize
    defaultFontSize = opts.fontSize;
    set(groot, 'defaultAxesFontSize',defaultFontSize,'defaultAxesFontweight','normal')
%     set(groot, 'defaultAxesTextSize',defaultFontSize,'defaultAxesFontweight','normal')
    set(groot, 'defaultLinelineWidth', opts.lineLineWidth)
    set(groot, 'defaultAxeslineWidth', opts.axesLineWidth)
    set(groot, 'DefaultFigurePosition', opts.position);
    set(groot, 'DefaultFigurewindowStyle', opts.windowStyle)    
    set(groot, 'DefaultFigureVisible', 'on')
end
end

function [opts] = parseinput(input, funcname, inputoffset)
% default values
opts = struct('codepath',pwd, 'isFigSet',false, 'position',[-1136 338 560 420], 'fontSize',20, ...
    'lineLineWidth',1.5, 'axesLineWidth',1.5, 'windowStyle','docked');

names = {'mainpath', 'codepath', 'isFigSet', 'position', 'fontSize', 'axesLineWidth', 'lineLineWidth', 'windowStyle'};
inputlen = length(input);

for i = 1:2:inputlen
    name = validatestring(input{i},names);    
    value = input{i+1};
    switch name
        case 'mainpath'
            validateattributes(value,{'string','char'}, {}, funcname, 'mainpath', i+1+inputoffset);
            opts.mainpath = value;
        case 'codepath'
            validateattributes(value,{'string','char'}, {}, funcname, 'codepath', i+1+inputoffset);
            opts.codepath = value;
        case 'isFigSet'
            validateattributes(value,{'numeric','logical'},{'scalar','integer'}, funcname, 'isFigSet', i+1+inputoffset)
            opts.isFigSet = value;
        case 'position'
            validateattributes(value,{'double'}, {'vector'}, funcname, 'position', i+1+inputoffset);
            opts.position = value;
        case 'fontSize'
            validateattributes(value,{'double'}, {'scalar'}, funcname, 'fontSize', i+1+inputoffset);
            opts.fontSize = value;
        case 'axesLineWidth'
            validateattributes(value,{'double'}, {'scalar'}, funcname, 'axesLineWidth', i+1+inputoffset);
            opts.axesLineWidth = value;
        case 'lineLineWidth'
            validateattributes(value,{'double'}, {'scalar'}, funcname, 'lineLineWidth', i+1+inputoffset)            
            opts.lineLineWidth = value;
        case 'windowStyle'
            validateattributes(value,{'string','char'}, {}, funcname, 'windowStyle', i+1+inputoffset);
            opts.windowStyle = value;            
    end
end
end