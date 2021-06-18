function negSLL = SLL_AR1(mu, sigma, rho, bfIrf_bc)
% This function is for MLE of parameters of spontaneous fluctuation (SF)
% It returns likelihood of mvnpdf(g(bfIrf),mu,Sigma) based on the two assumptions
% 1. SF_i ~ MVN(mu,Sigma), iid.
% 2. first-order auto regressive (AR(1)) (=Sigma can be expressed with the two terms, sigma and rho)
% NOTE: mu, sigma, rho of the function are parameters from TRANSFORMED data!

p = size(bfIrf_bc,1); % = bf.nSamples = sec*Hz
n = size(bfIrf_bc,2); % = # SF samples

covmat = NaN(p,p);
LL     = NaN(n,1);
idx = [p-1:-1:1 0 1:p-1];
for i = 1:p
    covmat(i,:) = sigma^2*rho.^idx(end+(-p-i+2):end+(-i+1));
end
for s = 1:n
    LL(s) = logmvnpdf(bfIrf_bc(:,s)', mu', covmat);
end

negSLL = -sum(LL);
end