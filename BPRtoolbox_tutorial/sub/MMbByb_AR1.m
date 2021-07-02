function negLL = MMbByb_AR1(params, bIrf, h, bc_avg, bc_Cov, lambda, g)
% it computes BPR scale for every single blink IRF 
Y = bIrf;
           
% data transformation with link function g and previously computed lambda
mvnP.mu    = bc_avg;
mvnP.cov   = bc_Cov;
               
x = Y-h*params;
% data transformation as well
x = g(x,lambda);
LL = logmvnpdf(x', mvnP.mu', mvnP.cov);

negLL = -LL;
end