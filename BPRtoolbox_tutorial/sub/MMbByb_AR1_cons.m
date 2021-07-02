function negLL = MMbByb_AR1_cons(params, bIrf, bOffIdx, dataLen, h, bc_avg, bc_Cov, lambda, g)

nConsB = numel(bOffIdx); % # consecutive blinks
% first, recover the full-length time course from the segments
bOffIdx = bOffIdx - bOffIdx(1) + 1; % set the first index as 1
BPRLen = length(bc_avg);

y              = NaN(dataLen,1);
SF_mu2_hats_bc = NaN(dataLen,nConsB); % matrix storing SF mu2 hat after each blink
% it computes BPR scale of the consecutive blinks at once
for i = 1:nConsB
    if i ~= nConsB
%         % Because index is approximated to the downsampled time course and h is 2.8 sec, Sometimes (e.g. the difference of indices = 2998)
%         if bOffIdx(i+1)-bOffIdx(i) <= size(bIrf,1)
            y(bOffIdx(i):bOffIdx(i+1)) = bIrf(1:bOffIdx(i+1)-bOffIdx(i)+1,i);
%         else
%             Y(bOffIdx(i):bOffIdx(i)+size(bIrf,1)-1) = bIrf;
%         end        
    else % if it is the last blink index
        y(bOffIdx(i):end) = bIrf(:,i);
    end
    
    % the length of the segment is different, so we need to expand the mu and sigma estimates.
    % NOTE1: I will expand the mu that assuming it is flat after first BPR lenghth, and covariance is just expanded.
    % NOTE2: the SF mu during consecutive blinks has been changed to mean of padded bc_avg (2019.06.10)
    SF_mu2_hats_bc(bOffIdx(i):(bOffIdx(i)+BPRLen-1),i) = bc_avg;
    % padding
    SF_mu2_hats_bc(1:bOffIdx(i),i) = bc_avg(1);
    SF_mu2_hats_bc((bOffIdx(i)+BPRLen-1):end,i) = bc_avg(end);        
end
% make BPRhat with convolution
BPRpulse = zeros(dataLen,1);
BPRpulse(bOffIdx) = params; % the amplitude of pulse is theta
BPRhat = conv(BPRpulse,h);
BPRhat = BPRhat(1:dataLen); % cut off the tail`


mvnP.mu    = mean(SF_mu2_hats_bc,2);
% just expand the covariance matrix with the sigma and rho already estimated
Sigma = bc_Cov(1); rho = bc_Cov(2)/Sigma;
mvnP.cov   = Sigma*getCorrMat(rho,dataLen);  

x = y-BPRhat;
% data transformation as well
x = g(x,lambda);
LL = logmvnpdf(x', mvnP.mu', mvnP.cov);

negLL = -LL;
end