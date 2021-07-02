function [covmat, corrmat] = getAR1cov(params, nRow)
% params(1,cond) = sigma, params(2,cond) = rho 
    % make MLE cov mat with the fitted parameters
    
    % initialize
    sizeMat = size(params);
    nCond   = prod(sizeMat(2:end));
    covmat  = cell(sizeMat(2:end));
    corrmat = cell(sizeMat(2:end));

    for cond = 1:nCond
        corrmat{cond} = getCorrMat(params(2,cond),nRow);
        covmat{cond} = params(1,cond)^2*corrmat{cond};
    end
    % if the intput was not a cell, output is also not a cell
    if ~iscell(params)
        covmat  = covmat{cond};
        corrmat = corrmat{cond};
    end
end