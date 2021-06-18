function corrmat = getCorrMat(rho,Len)
corrmat = NaN(Len,Len);
idx = [Len-1:-1:1 0 1:Len-1];
for i = 1:Len
    corrmat(i,:) = rho.^idx(end+(-Len-i+2):end+(-i+1));
end
end