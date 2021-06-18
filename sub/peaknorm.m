function yNorm = peaknorm(y)
% peak normalize event locked response
% y = T x N
yNorm = y./abs(min(y)); 
end

