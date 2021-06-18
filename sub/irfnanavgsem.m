function [avg_irf, sem_irf, std_irf] = irfnanavgsem(irf,varargin)
% Irf: a cell, or a matrix storing all Irfs of the condition
% varargin{1}: the dimension to mean, and std across
% NOTE: the first dimension should be time
% NOTE2: irfnanavgsem is excatly the same as irfavgsem except that it uses nanmean and nanstd
if iscell(irf)
    sizemat = size(irf{1});
else
    sizemat = size(irf);
end
dimmat = numel(sizemat);
dim = dimmat(end); % default: the last dimension
if nargin > 1
    dim = varargin{1};
end
    % initialize with IRF info of the first subject
    if iscell(irf)
        irfLen = size(irf{1},1);
        avg_irf = NaN([irfLen size(irf)]);
        sem_irf = NaN([irfLen size(irf)]);
        std_irf = NaN([irfLen size(irf)]);
        for i = 1:numel(irf)
            avg_irf(:,i) = squeeze(nanmean(irf{i},dim));
            std_irf(:,i) = squeeze(nanstd(irf{i},0,dim));
            sem_irf(:,i) = squeeze(std_irf(:,i)/sqrt(size(irf{i},dim)));
        end
    else
        avg_irf = squeeze(nanmean(irf,dim));
        std_irf = squeeze(nanstd(irf,0,dim));
        sem_irf = squeeze(std_irf/sqrt(size(irf,dim)));
    end
end