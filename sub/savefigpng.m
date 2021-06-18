function [] = savefigpng(varargin)
figname = varargin{end};
if nargin == 2
    path1 = varargin{1};
elseif nargin == 3 % case including subfolder
    path1 = fullfile(varargin{1}, varargin{2});
else
    error('check the number of input arguments!')
end
% make fig and png folder that store .fig and .png file respectively
mkdir(fullfile(path1,'fig')); mkdir(fullfile(path1,'png'))
save_figureToPNG(gcf,fullfile(path1,'png',figname));
savefig(gcf,fullfile(path1,'fig',figname));
