function [nRow, nCol] = minsquare(nElement)
% find the minimum square (nRow x nMin) larger than nElemnt
nCol = ceil(sqrt(nElement));
nRow = ceil(nElement/ceil(sqrt(nElement)));



