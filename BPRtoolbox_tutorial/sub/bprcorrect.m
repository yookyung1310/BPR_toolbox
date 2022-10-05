function [ts_corrected, SF, BPR, BPRfree, BPRaff, BPRcor, summary] = ...
    bprcorrect(ts, bOffIdx, bOffIdxAll, bOnIdxAll, sr, varargin)
%BPRCORRECT corrects BPR from pupil diameter time series
%Please read '  Yoo et al. (PLOS one, in revision)'
%   [ts_corrected, SF, BPR, BPRfree, BPRaff, BPRcor] = bprcorrect(ts, bOffIdx, bOffIdxAll, bOnIdxAll, sr, varargin)
%   ts: pupil DIAMETER time series (mm, or au) (cell or matrix)
%   bOffIdx: indices of intact blinks' offset in the time course (vector)
%   bOffIdxAll: indices of all blinks' offset in the time course (vector)
%   bOnIdxAll: indices of all blinks' onset in the time course (vector)
%   sr: sampling rate of the time course (scalar)
%   varargin: flags (see the function, parseinput)
%   input should be a single subject, or single condition data
% NOTE: explanations on output should be amended later
%   ts_corrected: pupil DIAMETER time series after BPR correction (mm, or au) 
%   SF: a structure storing variables about Spontaneous Fluctuation 
%   BPR: a structure storing variables about BPR estimates 
%   BPRfree: a structure storing variables about BPR-free periods
%   BPRaff: a structure storing variables about BPR-affected periods
%   BPRCor: a structure storing variables about BPR-corrected periods

%  Last update: 2022.10.05 Kyung Yoo
 
%   Copyright 2021 Kyung Yoo
%   Contatct: yookyung1310@gmail.com

% We recommned you to use a subject's own BPR profile, because the shape of BPR is slightly different across subject and
% background corneal flux density. If it is not possible to get an average BPR profile from your data 
% (not enough intact blinks, or other event is highly cross-correlated to blink, etc...), you can use canonical BPR from
% fixation task experiment of ours. Set BPRUseCano 'Yes' to apply canonical BPR to all subjects, or 'Auto' to apply 
% specific subjects who does not have enough samples or has so different BPR profile estimated.
% Yoo et al. (PLOS one, in revision) Please read on the paper for details.    

assert(mod(sr, 1) == 0, 'Sampling rate should be an integer') 

firstaxesinput = nargin-numel(varargin);
opts = parseinput(varargin, 'bprcorrect', firstaxesinput);
opts.path1 = fullfile(opts.SavePath,opts.CondName);
% flag to save figures
if ~isempty(opts.SavePath)
    isFigSave = true; % save figures
else
    isFigSave = false;
end
global ylabel1
xlabel1 = 'time from blink offset (sec)';
    
% range to get BPR profile
BPR.ConfRange  = [0 3];   % time span of BPR confounded (sec)
% BPR.ConfRange  = [-1 2];   % time span of BPR confounded (sec)
BPR.LockRange  = [-5 5];   % time span to show blink locked response (sec)
BPR.baseRange = [-0.1 0]; % time span of baseline (sec)
sr_rs = 5; % resample rate = 5Hz
BPR.tVec_rs = 0:1/sr_rs:BPR.ConfRange(2);    
BPR.crit.nBlinks  = 50;  % the number of blinks (all, or intact) should be equal or higher than this to estimate 
% BPR.crit.pdist    = 0.5; % cosine distance should be equal or lower than this to estimate 
% NOTE: Is it so controversial...? so changed at 20190507 by KY
BPR.crit.pdist    = 1.5; % cosine distance should be equal or lower than this to estimate 
BPR.crit.nBPRfree = 50; % the number of BPR-free spans should be equal or higher than this to estimate, if not use Lpass 

if iscell(ts), nRun = numel(ts); else, nRun = size(ts,2); end

%% Step 1: Get intact blink-locked response
% collect and concatenate all BPRs
BPR.intLockCat = []; % concatenate intact blink locked pupil time series
BPR.nIntBRun   = NaN(nRun,1);
for r = 1:nRun
    if iscell(ts), y = ts{r}; else, y = ts(:,r); end
    assert(numel(bOffIdxAll{r}) == numel(bOnIdxAll{r}), 'The number of all blinks'' onset (bOnIdxAll{run}) should be equal to all blinks'' offset (bOffIdxAll{run})') 
%     elrbInt_run = getIRF(y, bOffIdx{r}, BPR.LockRange, sr, 'NaN');
    elrbInt_run = getIRF(y, bOffIdx{r}, BPR.LockRange, sr, 'Last'); % changed from NaN to last value padding because NaN padding made the blr mean like stepwise (2019.05.28)
    BPR.intLockCat = [BPR.intLockCat elrbInt_run];
    BPR.nIntBRun(r) = size(elrbInt_run,2);    
end
BPR.nIntB = sum(BPR.nIntBRun);
fprintf('the number of intact blinks = %d\n',BPR.nIntB)

% all blink as well (2018.05.28)
BPR.allLockCat = []; % concatenate intact blink locked pupil time series
BPR.nAllBRun   = NaN(nRun,1);
BPR.ibiRun     = cell(nRun,1);
BPR.ibi = []; % concatenate interblink-interval (IBI)
for r = 1:nRun
    if iscell(ts), y = ts{r}; else, y = ts(:,r); end
%     elrbAll_run = getIRF(y, bOffIdxAll{r}, BPR.LockRange, sr, 'NaN');
    elrbAll_run = getIRF(y, bOffIdxAll{r}, BPR.LockRange, sr, 'Last'); % changed from NaN to last value padding because NaN padding made the blr mean like stepwise (2019.05.28)
    BPR.allLockCat = [BPR.allLockCat elrbAll_run];
    BPR.nAllBRun(r) = size(elrbAll_run,2);
    % get interblink interval. It was defined as span between blink offsets
    BPR.ibiRun{r} = diff(bOffIdxAll{r})/sr;
    BPR.ibi = [BPR.ibiRun{r}; BPR.ibi];
end
BPR.nAllB = sum(BPR.nAllBRun);
fprintf('the number of all blinks = %d\n',BPR.nAllB)

baseIdx = 1-BPR.LockRange(1)*sr;
if BPR.nIntB~=0
    % get intact blink-locked response
    BPR.iblr = nanmean(BPR.intLockCat,2); 
    BPR.iblrBcorr = BPR.iblr - BPR.iblr(baseIdx); % subtract baseline
else
    BPR.iblr = NaN(diff(BPR.LockRange)*sr+1,1);
    BPR.iblrBcorr = NaN(diff(BPR.LockRange)*sr+1,1);
    warning('no intact blink!')
end
% all blink as well
BPR.allblr = nanmean(BPR.allLockCat,2);
BPR.allblrBcorr = BPR.allblr - BPR.allblr(baseIdx); % subtract baseline
%     BPR.iblr(1:(-BPR.LockRange*esr)) = 0; % make 0 before blink offset
BPR.sr = sr;
if min(BPR.iblrBcorr)>=0 % If BPR profile is positive in all points, it's definetely weired
    error('Positive minimum of the profile. Please check the profile and consider use the known profile.')
end
    
% selected blink-locked response used to estimate parameters
if any(strcmpi(opts.LockIsl,{'Yes','Y'})) % if use intact-lock in all subjects
    BPR.selblr = BPR.iblr;
    BPR.nSelB = BPR.nIntB;
    BPR.strWhichLock = 'intact';
elseif any(strcmpi(opts.LockIsl,{'No','N'})) % if use all-lock in all subjects
    BPR.selblr = BPR.allblr;
    BPR.nSelB = BPR.nAllB;
    BPR.strWhichLock = 'all';
elseif any(strcmpi(opts.LockIsl,{'Auto','A'})) % if select automatically according to subject
    if BPR.nIntB >= BPR.crit.nBlinks  % if enough intact blinks
        BPR.selblr = BPR.iblr;
        BPR.nSelB = BPR.nIntB;
        BPR.strWhichLock = 'intact';
    else % if not enough
        BPR.selblr = BPR.allblr;
        BPR.nSelB = BPR.nAllB;
        BPR.strWhichLock = 'all';
    end
end
BPR.selblrBcorr = BPR.selblr - BPR.selblr(baseIdx); % subtract baseline

BPR.tVec_selblr = BPR.LockRange(1):1/BPR.sr:BPR.LockRange(2);
BPR.tVec    = 0:1/BPR.sr:BPR.ConfRange(2);
    
% show the number of intact/all blink, blink rate, and interblink-interval (IBI)
if iscell(ts)
    nSamplesRun = cellfun(@numel,ts); % # samples each run
    nSamples = sum(nSamplesRun(:)); % # samples total
else
    nSamples = numel(ts); % # samples total
end
BPR.brIntB = BPR.nIntB/nSamples*sr; % # blink/sec
BPR.brAllB = BPR.nAllB/nSamples*sr;

if ~opts.SuppressFig
    figure(99); clf;
    subplot(2,1,1); hold on;
    bar([0 1], [BPR.brAllB BPR.brIntB])
    xlim([-0.5 1.5]); xticks([0 1]); 
    xticklabels({'all','intact'})
    ShowText(sprintf('# all blinks: %d', BPR.nAllB), 0);
    ShowText(sprintf('# intact blinks: %d', BPR.nIntB), 1);
    title('blink rate')
    
    subplot(2,1,2); hold on;
    histogram(BPR.ibi);    
    
    
if isFigSave, savefigpng(opts.path1, 'profile', 'Blink properties'); end
end

% check the BPR profile and h
if ~opts.SuppressFig
    if BPR.nIntB~=0
        figure(100); clf;
        plot(BPR.tVec_selblr, BPR.selblr);
        ylabel(ylabel1)
    end
% show how many blinks are used for profile
if any(strcmpi(opts.UseCano,{'No','N'}))
    ShowText(sprintf('# %s blinks: %d', BPR.strWhichLock, BPR.nSelB));
end
title('Blink-locked mean')
if isFigSave, savefigpng(opts.path1, 'profile', 'Blink-locked mean'); end
end


%% Step 2: estimate parameters for SF (SF.lambda, mu_BPRfree, Sigma)
%%% Step 2-1: chunk and collect BPR-free data
critLen      = 3; % BPR is considered as remaining for 3 secs
nSamples_SF  = critLen*sr_rs+1; % # samples (length) in a single BPR or SF. If critLen=3 and sr_rs=5, # sample = 16 (0:0.2:3)
nSamples_org = critLen*sr+1;

BPRfree = getBPRfree(ts, bOffIdxAll, bOnIdxAll, sr, critLen);

% if a seesion could not record pupil at all, the pupil time course is full of zeros, and it was not processed at the
% preprocessing step. Thus, no blink, full of zeros pupil time course. It was fixed at 20190403, such that it has on
% e blink (the onset is start time and the offset is end time of the recording). 
% For the past data, it just throw out the zeros cases in BPR free irf time courses.
BPRfree.irf   = BPRfree.irf(:,~all(BPRfree.irf == 0));
BPRfree.slope = BPRfree.slope(:,~all(BPRfree.irf == 0));

% If the number of BPRfree span samples are too small, use low passed whole time course as a second-best solution.
% (2019.03.23 added, see Milestone 2019.01.16 slide for detail)
if size(BPRfree.irf,2) < BPR.crit.nBPRfree 
    % Store the orginial BPRfree, and make new BPRfree assuming that it can be approximated to low passed time course to
    % build the SF distribution later.
    BPRfree_org = BPRfree;
    % low-pass
    passCF    = 0.1; % (Hz)
    passOrder = 3;
    if iscell(ts), ts_lp = cell(size(ts)); else, ts_lp = NaN(size(ts)); end
    for r = 1:nRun
        if iscell(ts)
            ts_lp{r}   = bpfilter(passCF, sr, passOrder, ts{r}, 'low'); 
        else
            ts_lp(:,r) = bpfilter(passCF, sr, passOrder, ts(:,r), 'low');                 
        end
    end
%     figure(1); hold on; % plot to check
%     plot(ts); plot(ts_lp)
    
    % get BPRfree as blink did not occur at all
    BPRfree = getBPRfree(ts_lp, cell(size(ts_lp)), cell(size(ts_lp)), sr, critLen); 
    
    % if a seesion could not record pupil at all, the pupil time course is full of zeros, and it was not processed at the
    % preprocessing step. Thus, no blink, full of zeros pupil time course. It was fixed at 20190403, such that it has on
    % e blink (the onset is start time and the offset is end time of the recording).
    % For the past data, it just throw out the zeros cases in BPR free irf time courses.
    BPRfree.irf   = BPRfree.irf(:,~all(BPRfree.irf == 0));
    BPRfree.slope = BPRfree.slope(:,~all(BPRfree.irf == 0));
    
    % flag that it is low-passed and neglect blinks
    BPRfree.isLpassed = true;        
    
    fprintf('the number of BPRfree = %d\n ', size(BPRfree.irf,2))
    fprintf('the number criterion  = %d\n ', BPR.crit.nBPRfree)
    fprintf('It is not enough to estimate SF parameters, so use low passed time course as BPRfree\n')
end

% downsampling
% the beginning and end of the data is not zero, so padding and downsample
BPRfree.irf_rs = resampleirf(BPRfree.irf, sr, sr_rs); 


if ~opts.SuppressFig
figure(200); clf; 
subplot(2,1,1); hold on; 
plot(BPR.tVec, BPRfree.irf)
plot(BPR.tVec_rs, BPRfree.irf_rs)
subplot(2,1,2); hold on; 
plot(BPR.tVec, mean(BPRfree.irf,2))
plot(BPR.tVec_rs, mean(BPRfree.irf_rs,2))
if isFigSave, savefigpng(opts.path1, 'resampling', 'BPR-free irf resample'); end
end

%%% Step 2-2: box-cox transformation (estimate SF.lambda)
g = @(x,lambda) (x.^lambda-1)/lambda; % link function g
g_inv = @(x,lambda) exp(log(lambda*x+1)/lambda);
pupils = BPRfree.irf_rs(:); % all points
[~, SF.lambda] = boxcox(pupils);
BPRfree.irfbc_rs = g(BPRfree.irf_rs,SF.lambda); % box-cox transformed
pupilsbc = BPRfree.irfbc_rs(:); % all points

% BPR-free baseline histogram 
if ~opts.SuppressFig
figure(210); clf; 
histogram(pupils)
xlabel(ylabel1); ylabel('# samples');
titleStr = sprintf('pupil histogram'); title(titleStr);
if isFigSave, savefigpng(opts.path1, 'link function', titleStr); end
end

% the link function and box-cox data
xVec_bctmp = min(pupils):0.001:max(pupils);
yVec_bctmp = g(xVec_bctmp,SF.lambda);
if ~opts.SuppressFig
figure(211); clf; hold on;
plot(xVec_bctmp,yVec_bctmp)
xlabel(ylabel1); ylabel('transformed');
titleStr = sprintf('link function transformation'); title(titleStr);
if isFigSave, savefigpng(opts.path1, 'link function', titleStr); end
end

% normal distribution approximation to AFTER data transformation
% I will approximate normpdf to box-cox transformed data, and return to original space with the inverse link function
% see how different the approximation and data acquired are
if ~opts.SuppressFig
figure(212); clf; hold on;
% data
histogram(BPRfree.irfbc_rs,'Normalization','pdf')
% approximation (marginal normal distributions)
bVec_tmp = linspace(min(pupilsbc),max(pupilsbc),100);
plot(bVec_tmp,normpdf(bVec_tmp,mean(pupilsbc),std(pupilsbc)))
axis tight;
xlabel('transformed PD'); ylabel('pdf');
titleStr = sprintf('transformed data and normal fit'); title(titleStr);
if isFigSave, savefigpng(opts.path1, 'link function', titleStr); end
end

%%% Step 2-3: estimate mu_BPRfree and Sigma (sigma & rho)
% Two core assupmtions: 1. g(SF_i) ~ MVN(mu_BPRfree, Sigma)  
%                       2. g(SF_i(j+1)) = sqrt(rho)*g(SF_i(j)) + e + mu_BPRfree(j+1)-mu_BPRfree(j), e~N() (AR(1))
% other assumptio.ns: SF_i is iid., SF_i and e are indepedent
% NOTE: mu_BPRfree and Sigma is parameters for G-TRANSOFORMED BPR-free data
SF.mu_hat = mean(BPRfree.irfbc_rs,2); % sample mean is just equal to mean estimate

% - cov mat estimation -
% NOTE: cov matrix is acquired after making the distribution normal by applying one's own g, 
% so each condition has different scale of course
BPRfree.Scov_bc  = cov(BPRfree.irfbc_rs'); % sample covariance matrix in the transfromed space
BPRfree.Scorr_bc = corrcoef(BPRfree.irfbc_rs'); % sample correlation matrix

             % sigma                     rho
initParams = [sqrt(BPRfree.Scov_bc(1,1)) 0.95]';
minParams  = [0                          0   ]';
maxParams  = [inf                        1   ]';

% estimate covariance matrix from gen BPR-free data with the 2 assumptions
SLL_AR1_f1 = @(p) SLL_AR1(SF.mu_hat, p(1), p(2), BPRfree.irfbc_rs);
%             options = optimset('MaxFunEvals',100000, 'MaxIter',3000);
fprintf('-- Covariance matrix (Sigma) estimation start --\n')
[SF.AR1covP.params, SF.AR1covP.fitRes.fval, SF.AR1covP.fitRes.exitflag, SF.AR1covP.fitRes.output] = ...
    fminsearchbnd(SLL_AR1_f1, initParams, minParams, maxParams);
fprintf('-- Covariance matrix (Sigma) estimation end --\n')
[SF.cov_hat, SF.corr_hat] = getAR1cov(SF.AR1covP.params, nSamples_SF);

% check the result
% compare sample cov mat and AR(1) constrained MLE of cov mat
clim1 = [min(BPRfree.Scov_bc(:)) max(BPRfree.Scov_bc(:))]; % set color limit
if ~opts.SuppressFig
figure(220); clf; 
subplot(1,2,1); % sample cov mat
imagesc(BPRfree.Scov_bc, clim1);
colormap('Gray'); colorbar;
subplot(1,2,2); % AR(1) MLE of cov mat
imagesc(SF.cov_hat, clim1);
colormap('Gray'); colorbar;
titleStr = sprintf('sample cov VS AR1 MLE cov'); mtit(titleStr);
if isFigSave, savefigpng(opts.path1, 'covmat', titleStr); end
end

% corr mat as well
% NOTE: it was computed from TRANSOFRMED data
if ~opts.SuppressFig
figure(221); clf; 
clim1 = [min(BPRfree.Scorr_bc(:)) max(BPRfree.Scorr_bc(:))]; % set color limit
subplot(1,2,1); % sample cov mat
imagesc(BPRfree.Scorr_bc, clim1);
colormap('Gray'); colorbar;
subplot(1,2,2); % AR(1) MLE of cov mat
imagesc(SF.corr_hat, clim1);
colormap('Gray'); colorbar;
titleStr = sprintf('sample corr VS AR1 MLE corr'); title(titleStr);
if isFigSave, savefigpng(opts.path1, 'covmat', titleStr); end
end

%%% Step 2-4: estimate mu_BPRaff and BPR kernel (h)
% Assume that covariance of SF is the same whether blink occur or not, but only mean time course is different.
% 
% fit the mu_BPRaff and h to intact, or all-blink locekd mean (profile)

x = BPR.LockRange(1):1/BPR.sr:BPR.LockRange(end);
y = BPR.selblr;
startIdx = (-BPR.LockRange(1)+BPR.ConfRange(1))*sr+1; 
endIdx   = (-BPR.LockRange(1)+BPR.ConfRange(2))*sr+1; % point at 0 and P.wSize1 (3) sec
% modified at 20190507 by KY
% contaSec = 0.15; % consider 0~0.15 is contaminated by artifact.
% startTime = near(contaSec,x); % the neareast point from P.conta
startTime = near(0,x); % the neareast point from P.conta
% startTime = near(-0.15,x); % get offset before positive peak
isoffSet = x==startTime; 
yOffset = y(isoffSet);

% get slopes of two points
tmp = diff(x); xtic = tmp(1);
% Let's call the slopes at P.conta and P.wSize1 b1 and b2
% method 1-2: use the slope at reliable point (e.g. -1 sec) as f'(0)
yprev1 = BPR.selblr(startIdx-1*sr+1); % just to get the derivative at -1 sec
yprev2 = BPR.selblr(startIdx-1*sr);
yend1  = BPR.selblr(endIdx);
yend2  = BPR.selblr(endIdx-1);
b1 = (yprev1-yprev2)/xtic;
b2 = (yend1-yend2)/xtic;

% initial parameter from Exp1. result
Exp1Res = load('cano/cano_p','cano_p'); 
% sbj 1-14 + mean across the subjects (=canonical h)
initParams = [Exp1Res.cano_p mean(Exp1Res.cano_p,2)];
initParams(4,1:15)  = 0;
initParams(4,16:30) = 0.1;
initParams(4,31:45) = 0.2;
initParams(4,46:60) = 0.3;
if any(strcmpi(opts.Ismm, {'Yes','Y'}))
    minParams  = [0;      0;      0.01;    0];
    maxParams  = [inf;    0.4;    3;       0.3];
elseif any(strcmpi(opts.Ismm, {'No','N'}))
    minParams  = [0;      0;      1;       0];
    maxParams  = [inf;    0.4;    30000;   0.3];    
% if not defined, regard it as mm unit if pupil value mean < 10
elseif any(strcmpi(opts.Ismm, {'Auto','A'})) 
    if nanmean(BPR.selblr) < 10
        BPR.regard_as_mm = 'Yes';
        minParams  = [0;      0;      0.01;    0];
        maxParams  = [inf;    0.4;    3;       0.3];
    else
        BPR.regard_as_mm = 'No';
        minParams  = [0;      0;      1;       0];
        maxParams  = [inf;    0.4;    30000;   0.3];
    end
end
nipm = size(initParams,2);
npm  = size(initParams,1);



options = optimset('MaxFunEvals',numel(initParams)*200, 'MaxIter',numel(initParams)*200); % roughly estimate
% options = optimset('MaxFunEvals',numel(initParams)*2000, 'MaxIter',numel(initParams)*2000); % x10 estimate

% options = optimset('MaxFunEvals',numel(initParams)*2000, 'MaxIter',numel(initParams)*2000, 'Display','iter'); 
% options = optimset('MaxFunEvals',numel(initParams)*2000, 'MaxIter',numel(initParams)*2000, 'Display','iter', ...
%     'TolFun',1.e-13,'TolX',1.e-13); 

% % old canonical h
% load canonical BPR data from the fixation experiment
% load(fullfile(getNASdir,'/people/KY/PBPC_2017/merged/profileTS6.mat')) % change the pathway later
% prepare canonical BPR, variables, and functions for fitting
% dblgamma = @(p) p(5)*(p(6)*gampdf(x2,p(1),p(2)) - gampdf(x2,p(3),p(4)));
% BPR.h_cano = dblgamma(pf.cano.params);
% new canonical h (2019.06.05)
sglgamma = @(p) p(3)* -gampdf(BPR.tVec_rs-p(4),p(1),p(2));
BPR.h_cano = sglgamma([mean(Exp1Res.cano_p,2); 0]);

% % use downsampled data and estiamte to decrease burden of computation 
% % cut the span to begin at offset timing for h
% y = BPR.selblr(startIdx:endIdx);
% y1 = y(end);    
% y2 = resampleirf2(y, sr, sr_rs);
% cov_hat2 = SF.cov_hat/BPR.nIntB; % it is mean time course so variance is decreased 
% dblgamma = @(p) p(5)*(p(6)*gampdf(x2,p(1),p(2)) - gampdf(x2,p(3),p(4)));
% getSFmeth1_f1 = @(p) getSFmeth1(x2, y1, p(1), p(2), yOffset);
% % minus log-likelihood of kernel and mu_BPRaff
% negLLhSF2 = @(p) -logmvnpdf(g(y2'-dblgamma(p(1:6)), SF.lambda), g(getSFmeth1_f1(p(7:8)), SF.lambda), cov_hat2); 

% use larger sampling rate (than resample rate, sr_rs) to avoid weird sharp h estimation
sr_hSF2 = 20; % (Hz)
BPR.sr_hSF2 = sr_hSF2; % save it
nSamples_hSF2 = critLen*sr_hSF2+1;
BPR.tVec_hSF2 = 0:1/BPR.sr_hSF2:BPR.ConfRange(2);
x3 = BPR.tVec_hSF2;
y = BPR.selblr(startIdx:endIdx);
y1 = y(end);    
y3 = resampleirf2(y, sr, sr_hSF2);

% cubic spline to substitue SF mu2
x = BPR.LockRange(1):1/BPR.sr:BPR.LockRange(end);
y = BPR.selblr;
y2 = y;
bprSpan = ((-BPR.LockRange(1)+BPR.ConfRange(1))*BPR.sr):((-BPR.LockRange(1)+BPR.ConfRange(2))*BPR.sr); % BPR lasts 3 sec
y2(bprSpan) = NaN; % make NaN the BPR affected period
yy = spline(x,y2,x);
SF.mu2_hat_org      = yy(bprSpan);
SF.mu2_hat_hEst     = resampleirf2(SF.mu2_hat_org', sr, sr_hSF2);
SF.mu2_hat          = resampleirf2(SF.mu2_hat_org', sr, sr_rs);
% save it as box-cox transformed as convention
SF.mu2_hat_org      = g(SF.mu2_hat_org, SF.lambda);
SF.mu2_hat_hEst     = g(SF.mu2_hat_hEst, SF.lambda);
SF.mu2_hat          = g(SF.mu2_hat, SF.lambda);

% Sigma and rho, rho should be converted according to the the new smapling rate (hEst)
covParams = NaN(2,1);
covParams(1) = SF.AR1covP.params(1); % sigma
covParams(2) = SF.AR1covP.params(2)^(1/(sr/sr_hSF2)); % rho converted
cov_hat3 = getAR1cov(covParams, nSamples_hSF2)/BPR.nSelB;
sglgamma = @(p) p(3)* -gampdf(x3-p(4),p(1),p(2));

% minus log-likelihood of kernel and mu_BPRaff
negLLhSF2 = @(p) -logmvnpdf(g(y3'-sglgamma(p(1:4)), SF.lambda), SF.mu2_hat_hEst', cov_hat3); 

% modified at 20190507 by KY
% negLLhSF2_f1 = @(p) negLLhSF2(p, y3, SF.lambda, cov_hat3, g, getSFmeth1_f1, dblgamma, pf.cano.params, BPR.tVec, BPR.crit.pdist);
% negLLhSF2_f1 = @(p) negLLhSF2(p, y3, SF.lambda, cov_hat3, g, getSFmeth1_f1, dblgamma, pf.cano.params, BPR.tVec, 0.3);



% fminCount = 0; 
% shape1nStep = 5; shape2nStep = 5;
sglgamma = @(p) p(3)* -gampdf(BPR.tVec-p(4),p(1),p(2));
% dblgamma = @(p) p(5)*(p(6)*gampdf(BPR.tVec,p(1),p(2)) - gampdf(BPR.tVec,p(3),p(4)));

clearvars tmp
if any(strcmpi(opts.UseCano,{'No','N'}))
    % SF + single gamma fitting  
    for ipm = 1:nipm
%         [BPR.hSF2.params(:,ipm), BPR.hSF2.fitRes.fval(:,ipm),  BPR.hSF2.fitRes.exitflag(:,ipm), BPR.hSF2.fitRes.output(:,ipm)] = ...
%             fminsearchbnd(negLLhSF2, initParams(:,ipm), minParams, maxParams, options);
               
        [BPR.hSF2.params(:,ipm), BPR.hSF2.fitRes.fval(ipm),  BPR.hSF2.fitRes.exitflag(ipm), ...
            BPR.hSF2.fitRes.output(ipm), BPR.hSF2.fitRes.history{ipm}] = ...
            myproblem(negLLhSF2, initParams(:,ipm), minParams, maxParams, options);
        
%         h_hat = sglgamma(BPR.hSF2.params(1:npm,ipm))';
    end
    [~,sIdx] = sort(BPR.hSF2.fitRes.fval);
    h_hat = sglgamma(BPR.hSF2.params(1:npm,sIdx(1)))';
    
    BPR.isCanoUsed = false;
elseif any(strcmpi(opts.UseCano,{'Yes','Y'}))    
    BPR.isCanoUsed = true;
elseif any(strcmpi(opts.UseCano,{'Auto','A'}))
    if BPR.nSelB >= BPR.crit.nBlinks % if enough selected (intact or all) blinks for estimation
        % SF + single gamma fitting
        for ipm = 1:nipm
            [BPR.hSF2.params(:,ipm), BPR.hSF2.fitRes.fval(ipm),  BPR.hSF2.fitRes.exitflag(ipm), ...
                BPR.hSF2.fitRes.output(ipm), BPR.hSF2.fitRes.history{ipm}] = ...
                myproblem(negLLhSF2, initParams(:,ipm), minParams, maxParams, options);
            
%             h_hat = sglgamma(BPR.hSF2.params(1:npm,ipm))';
        end
        [~,sIdx] = sort(BPR.hSF2.fitRes.fval);
        h_hat = sglgamma(BPR.hSF2.params(1:npm,sIdx(1)))';
    
%         % compute pairwise cosine distance, and throw away h hat if too different. 
%         BPR.h_cano = dblgamma(pf.cano.params)';
%         cosdist = pdist2(BPR.h_cano', h_hat','cosine');
%         BPR.cosdist_beforeLine = cosdist;
%         if cosdist < BPR.crit.pdist
%             BPR.isCanoUsed = false;
%             BPR.isSFline = false;
%         else                        
% %             BPR.isCanoUsed = true;
% %             disp('h is too different from canonical BPR') % old
%             disp('h is too different from canonical BPR, re-estimate h assuming that SFmu_2 is linear.')
%             BPR.isSFline = true;
%             
%             [BPR.hSF2.params, BPR.hSF2.fitRes.fval, BPR.hSF2.fitRes.exitflag, BPR.hSF2.fitRes.output] = ...
%                 fminsearchbnd(negLLhSF2_line, initParams(1:npm), minParams(1:npm), maxParams(1:npm), options);
% %             [BPR.hSF2.params, BPR.hSF2.fitRes.fval, BPR.hSF2.fitRes.exitflag, BPR.hSF2.fitRes.output] = ...
% %                 fminunc(negLLhSF2_line, initParams(1:npm), options);
%             
%             SFline_up = SFline_coeff.slope*BPR.tVec + SFline_coeff.intercept;
%             SF.mu2_hat_org = g(SFline_up, SF.lambda);
%             h_hat = sglgamma(BPR.hSF2.params(1:npm))';
%             
%             cosdist = pdist2(BPR.h_cano', h_hat','cosine');
%             if cosdist < BPR.crit.pdist
%                 BPR.isCanoUsed = false;
%             else
%                 BPR.isCanoUsed = true;
%             end
%             
%         end
    else        
        BPR.isCanoUsed = true;
        fprintf('# %s blinks %d is less than criterion %d, so replace h with canonical BPR, and SF_mu2 hat is flat', ...
            BPR.strWhichLock, BPR.nSelB, BPR.crit.nBlinks)
    end
end
% if canonical BPR is used, give up estimating SF so just consider it flat.
if BPR.isCanoUsed    
    SF.mu2_hat_org = mean(SF.mu_hat)*ones(length(BPR.tVec),1);
%     BPR.h_org = sglgamma(pf.cano.params)'; % the kernel is now canonical BPR         
    BPR.h_org = BPR.h_cano;
else    
    BPR.h_org = h_hat;
end
fprintf('SF in BPR and h estimation finished\n')
% scaling the h such that minimum of h to be -0.1mm
BPR.h = peaknorm(BPR.h_org)*0.1;
normalizer = min(BPR.h)/min(BPR.h_org);
% downsampled h for computation
sglgamma = @(p) p(3)* -gampdf(BPR.tVec_rs-p(4),p(1),p(2));
% dblgamma = @(p) p(5)*(p(6)*gampdf(BPR.tVec_rs-p(4),p(1),p(2)) - gampdf(BPR.tVec_rs,p(3),p(4)));

if ~BPR.isCanoUsed
    BPR.h_rs    = sglgamma(BPR.hSF2.params(1:npm))'*normalizer; % set the amplitude same as h
    
else
%     BPR.h_rs    = dblgamma(pf.cano.params)'*normalizer; % set the amplitude same as h
end

if ~opts.SuppressFig
figure(101); clf; hold on;
plot(BPR.tVec, BPR.h);
plot(BPR.tVec_rs, BPR.h_rs);
titleStr = sprintf('kernel resample'); title(titleStr);
legend('kernel (h)','kernel (h) resampled')
if isFigSave, savefigpng(opts.path1, 'resampling', 'kernel resample'); end
end

if ~opts.SuppressFig
    figure(102); clf; hold on;
    obj1 = plot(BPR.tVec_selblr, BPR.selblr, 'k');    
    if ~BPR.isCanoUsed
        x2 = BPR.tVec_rs;
%         getSFmeth1_f1 = @(p) getSFmeth1(BPR.tVec, y1, p(1), p(2), yOffset);
%         SF_mu_ginved = getSFmeth1_f1(BPR.hSF2.params(4:5))';
%         obj2 = plot(BPR.tVec,SF_mu_ginved);
        
        h_hat_rs = g_inv(SF.mu2_hat_hEst, SF.lambda);
        SF_mu_ginved = g_inv(SF.mu2_hat_hEst, SF.lambda);
        obj3 = plot(BPR.tVec_hSF2, SF_mu_ginved, '--b'); % downsampled SF_mu which was used to estimate SF_mu2 and h
                
        offset1 = SF_mu_ginved(1); % offset just to overlap
        obj4 = plot(BPR.tVec, BPR.h_org+offset1, 'r');
        
        startIdx = -BPR.LockRange(1)*sr+1; endIdx = (-BPR.LockRange(1)+BPR.ConfRange(2))*sr+1; % point at 0 and P.wSize1 (3) sec
        y = BPR.selblr(startIdx:endIdx);
        
        obj5 = plot(BPR.tVec, y-BPR.h_org);
        %     xlim([-6 6]);
        xlim([-2 4]);
        %     xlim([0 3]);
        
        hline(offset1)
        xlabel(xlabel1); ylabel(ylabel1);
        legend([obj1 obj3 obj4 obj5], ...
            'data','SF mu2 hat (g inv)','h hat (pre norm)','data-h hat')
        title('h and SF_mu2(g inversed) fit to blink-locked time course')
        
%         % data - SF_mu vs BPR.h_org
%         figure(1020); clf; hold on;
%         plot(BPR.tVec, BPR.selblr(startIdx:endIdx)-BPR.selblr(startIdx), 'k');
%         plot(BPR.tVec, BPR.selblr(startIdx:endIdx)-SF_mu_ginved);
%         plot(BPR.tVec, BPR.h_org);
%         legend('data (baselined)','data-SF mu hat','h hat (pre norm)')
    else
        legend(obj1, 'data')
        title('blink-locked time course')
    end
    if isFigSave, savefigpng(opts.path1, 'profile', 'h and SF mean estimates'); end
end

%% Step 3: estimate BPR amplitudes (theta) blink-by-blink
%%% Step 3-1: get BPR-affected periods and downsample as well
BPRaff.irf = cell(nRun,1);  BPRaff.irf_rs = cell(nRun,1);
for r = 1:nRun
    if iscell(ts), y = ts{r}; else, y = ts(:,r); end
    BPRaff.irf{r} = getIRF(y, bOffIdxAll{r}, [0 critLen], sr, 'Last');
%     nSamples_run = length(y); % # samples in the run        
%     for b = 1:numel(bOffIdxAll{r})
%         % Out of boundary case is very rare, but if it is, just pad the last value
%         if bOffIdxAll{r}(b)+sr*critLen <= nSamples_run 
%             BPRaff.irf{r}(:,b) = y(bOffIdxAll{r}(b):bOffIdxAll{r}(b)+sr*critLen);
%         else % Out of boundary case
%             padLen = sr*critLen - (nSamples_run - bOffIdxAll{r}(b));
%             BPRaff.irf{r}(:,b) = [y(bOffIdxAll{r}(b):nSamples_run); y(nSamples_run)*ones(padLen,1)];
%         end 
%     end
    if ~isempty(bOffIdxAll{r})
        BPRaff.irf_rs{r} = resampleirf(BPRaff.irf{r}, sr, sr_rs);
    else % no blink in the run
        BPRaff.irf_rs{r} = NaN(nSamples_SF,0);
    end
end

%%% Step 3-2: estimate BPR amplitudes (theta)
initParams = abs(min(BPR.selblrBcorr))*10; % the adjusted scale from profile (= mean BPR-affected)
minParams  = 0;
maxParams  = inf;
    
BPRaff.theta_hat = cell(nRun,1);
nBlinks_run      = NaN(nRun,1);
for r = 1:nRun
    BPRaff.theta_hat{r} = zeros(1,size(BPRaff.irf{r},2));
    nBlinks_run(r) = size(BPRaff.irf_rs{r},2);
    isProcessed = false(nBlinks_run(r),1); % flag whether the blink is processed
    isAffecting = [(diff(bOffIdxAll{r})/sr) < critLen; 0]; % Does this blink's BPR affect its subsequent blink(s)?
    consStartorEnd = find(xor(isAffecting, [0; isAffecting(1:end-1)]));
    BPR.consStartB{r} = consStartorEnd(1:2:end);
    BPR.consEndB{r}   = consStartorEnd(2:2:end);
    BPR.isSep{r}      = false(nBlinks_run(r),1);
    BPR.isCons{r}     = false(nBlinks_run(r),1); % separated & consecutive blinks are exactly opposite, but record them at the moment
    BPR.computeTime{r}  = NaN(nBlinks_run(r),1);
    
    consCount = 0;
    for b = 1:nBlinks_run(r)        
        % 2018.04.03 KY added. Correct the code to deal with consecutive blinks cases
        if ~isProcessed(b)
            % if the blinks are separated enough
            if ~ismember(b,BPR.consStartB{r})
                tic; % compute time (start)
                
                MMbByb_AR1_f1 = @(p) MMbByb_AR1(p, BPRaff.irf_rs{r}(:,b), BPR.h_rs, ...
                    SF.mu2_hat, SF.cov_hat, SF.lambda, g);
                [BPRaff.theta_hat{r}(b), BPRaff.fitRes{r}.fval(b), BPRaff.fitRes{r}.exitflag(b), BPRaff.fitRes{r}.output(b)] = ...
                    fminsearchbnd(MMbByb_AR1_f1, initParams, minParams, maxParams);
                isProcessed(b) = true;      
                
                BPR.isSep{r}(b) = true;                                
                BPR.computeTime{r}(b) = toc; % compute time (end)
            else                                
                % treat consecutive blinks
                % LL of the consecutive blinks should be computed at once.
                % b is the first blink of the conescuetive blinks, and b2 is all of the consecutive blinks                
                tic; % compute time (start)
                
                consCount = consCount+1;                
                b2 = BPR.consStartB{r}(consCount):BPR.consEndB{r}(consCount); % the consecutive blinks
                BPR.nCons{r}(b2) = numel(b2);
                % (critLen-1/sr_rs) = the lentgh of h (2.8sec)
                bOffIdx_rs = round((bOffIdxAll{r}(b2)-1)/sr*sr_rs); % approximation of the blink offset timings because of the downsampling
                dataLen    = bOffIdx_rs(end)-bOffIdx_rs(1)+1+critLen*sr_rs;
                
                MMbByb_AR1_f2 = @(p) MMbByb_AR1_cons(p, BPRaff.irf_rs{r}(:,b2), bOffIdx_rs, dataLen, BPR.h_rs, ...
                    SF.mu2_hat, SF.cov_hat, SF.lambda, g);
                [BPRaff.theta_hat{r}(b2), fval_tmp, exitflag_tmp, output_tmp] = ...
                    fminsearchbnd(MMbByb_AR1_f2, initParams*ones(1,numel(b2)), minParams*ones(1,numel(b2)), maxParams*ones(1,numel(b2)));
                
                BPRaff.fitRes{r}.fval(b2)     = fval_tmp;
                BPRaff.fitRes{r}.exitflag(b2) = exitflag_tmp;
                BPRaff.fitRes{r}.output(b2)   = output_tmp;
                
                isProcessed(b2) = true;
                
                BPR.isCons{r}(b2) = true;
                
                BPR.computeTime{r}(b2(end)) = toc; % compute time (end)
            end
        end
    end
    BPR.consCount(r) = consCount;
    
    fprintf('run%d theta estimation finished\n',r)
end


%% Step 4: Correction with the BPR amplitude estimates
% replace the BPR-affected periods with the residuals
if iscell(ts), ts_corrected = cell(size(ts)); else, ts_corrected = NaN(size(ts)); end    
for r = 1:nRun
    if iscell(ts), y = ts{r}; else, y = ts(:,r); end    
    % make BPRhat with convolution
    dataLen = size(y,1);
    BPRpulse = zeros(dataLen,1);
    BPRpulse(bOffIdxAll{r}) = BPRaff.theta_hat{r}; % the amplitude of pulse is theta
    BPRhat = conv(BPRpulse,BPR.h);
    BPRhat = BPRhat(1:dataLen); % cut off the tail
    if iscell(ts), ts_corrected{r} = y-BPRhat; else, ts_corrected(:,r) = y-BPRhat; end
end



% compare the BPR correction result!
for r = 1:10:nRun % every 10 runs    
    if ~opts.SuppressFig
    figure(400+r); clf; hold on;
    if iscell(ts), y = ts{r}; y2 = ts_corrected{r}; ...
    else, y = ts(:,r); y2 = ts_corrected(:,r); end
    if ~isempty(y)
        tVecRun = 0:1/sr:length(y)/sr-1/sr;
        obj1 = plot(tVecRun,y);
        obj2 = plot(tVecRun,y2);
        if ~isempty(bOffIdxAll{r})
            obj3 = vline(bOffIdxAll{r}/sr);
            legend([obj1, obj2, obj3(1)], 'pre-correction','post-correction','blink (offset)');
        else % no blink case
            legend([obj1, obj2], 'pre-correction','post-correction');
        end
    end
    xlabel('time (sec)'); ylabel(ylabel1);
    title(sprintf('run %d',r))
    if isFigSave, savefigpng(opts.path1, fullfile('BPR correction (time series)',sprintf('run%d',r)), titleStr); end
    end
end

% get corrected irf blink-by-blink
BPRcor.irf = cell(nRun,1);  BPRcor.irf_rs = cell(nRun,1);
for r = 1:nRun
    if iscell(ts_corrected), y = ts_corrected{r}; else, y = ts_corrected(:,r); end
    BPRcor.irf{r} = getIRF(y, bOffIdxAll{r}, [0 critLen], sr, 'Last');
    if ~isempty(bOffIdxAll{r})
        BPRcor.irf_rs{r} = resampleirf(BPRcor.irf{r}, sr, sr_rs);
    else % no blink in the run
        BPRcor.irf_rs{r} = NaN(nSamples_SF,0);
    end
end

% treat empty matrices in BPRcor
for r = 1:nRun
    if isempty(BPRcor.irf_rs{r})
        BPRcor.irf_rs{r} = []; BPRcor.irf{r} = [];
    end
end
r = 1;
if ~opts.SuppressFig
figure(300); clf; hold on;
plot(BPR.tVec, BPRcor.irf{r})
plot(BPR.tVec_rs, BPRcor.irf_rs{r})
titleStr = sprintf('residual resample (run %d)',r); title(titleStr);
xlabel(xlabel1); ylabel(ylabel1);
if isFigSave, savefigpng(opts.path1, 'resampling', 'residual resample'); end
end

% get avg, sem, and std of them
BPRaff_irfCat = cat(2,BPRaff.irf{:}); % the concatenated variables are redundant so it would be only used temporally.
BPRaff_irf_rsCat = cat(2,BPRaff.irf_rs{:});
BPRcor_irfCat = cat(2,BPRcor.irf{:});
BPRcor_irf_rsCat = cat(2,BPRcor.irf_rs{:});

[BPRaff.irf_Avg,  BPRaff.irf_Sem,  BPRaff.irf_Std]  = irfavgsem(BPRaff_irfCat);
[BPRfree.irf_Avg, BPRfree.irf_Sem, BPRfree.irf_Std] = irfavgsem(BPRfree.irf); % it was already concatenated
[BPRcor.irf_Avg,  BPRcor.irf_Sem,  BPRcor.irf_Std]  = irfavgsem(BPRcor_irfCat);
[BPRaff.irf_rsAvg,  BPRaff.irf_rsSem,  BPRaff.irf_rsStd]  = irfavgsem(BPRaff_irf_rsCat);
[BPRfree.irf_rsAvg, BPRfree.irf_rsSem, BPRfree.irf_rsStd] = irfavgsem(BPRfree.irf_rs); % it was already concatenated
[BPRcor.irf_rsAvg,  BPRcor.irf_rsSem,  BPRcor.irf_rsStd]  = irfavgsem(BPRcor_irf_rsCat);
% get avg, sem and std with baselined time courses
[BPRaff.bl_irf_Avg,  BPRaff.bl_irf_Sem,  BPRaff.bl_irf_Std]  = irfavgsem(BPRaff_irfCat - BPRaff_irfCat(1,:));
[BPRfree.bl_irf_Avg, BPRfree.bl_irf_Sem, BPRfree.bl_irf_Std] = irfavgsem(BPRfree.irf - BPRfree.irf(1,:)); % it was already concatenated
[BPRcor.bl_irf_Avg,  BPRcor.bl_irf_Sem,  BPRcor.bl_irf_Std]  = irfavgsem(BPRcor_irfCat - BPRcor_irfCat(1,:));
[BPRaff.bl_irf_rsAvg,  BPRaff.bl_irf_rsSem,  BPRaff.bl_irf_rsStd]  = irfavgsem(BPRaff_irf_rsCat - BPRaff_irf_rsCat(1,:));
[BPRfree.bl_irf_rsAvg, BPRfree.bl_irf_rsSem, BPRfree.bl_irf_rsStd] = irfavgsem(BPRfree.irf_rs - BPRfree.irf_rs(1,:)); % it was already concatenated
[BPRcor.bl_irf_rsAvg,  BPRcor.bl_irf_rsSem,  BPRcor.bl_irf_rsStd]  = irfavgsem(BPRcor_irf_rsCat - BPRcor_irf_rsCat(1,:));

% - check the result -
% blink-by-blink correction result
for r = [1 nRun] % first and last run only
    if nBlinks_run(r) > 0
        blinkSel = unique(floor(linspace(1,nBlinks_run(r),10))); % plot 10 blinks within the run
        for b = blinkSel
            y_data  = BPRaff.irf_rs{r}(:,b);
            y_BPRhat = BPR.h_rs*BPRaff.theta_hat{r}(b);
            
            if ~opts.SuppressFig
                figure(30000+b); clf; hold on;
                plot(BPR.tVec_rs, y_data, 'k');
                plot(BPR.tVec_rs, y_data-y_BPRhat, 'r');
                plot(BPR.tVec_rs, y_data(1)+y_BPRhat, 'r--', 'linewidth', 1);
                xlabel(xlabel1); ylabel(ylabel1);
                legend('data','data-BPRhat');
                titleStr = sprintf('data and BPR-corrected (%dth blink)',b); title(titleStr);
                if isFigSave, savefigpng(opts.path1, fullfile('BPR correction (blinks)',sprintf('run%d',r)), titleStr); end
            end
        end        
    end
    %         close all
end

% put correction result summary into summary structure (2019.07.03)
% fields added (2019.08.28)
summary.tVec_selblr     = BPR.tVec_selblr;
summary.selblr          = BPR.selblr;
summary.tVec            = BPR.tVec;
summary.mu2_hat         = SF.mu2_hat;
summary.mu2_hat_org     = SF.mu2_hat_org;
summary.lambda          = SF.lambda;
summary.h               = BPR.h;
summary.h_org           = BPR.h_org;
summary.h_rs            = BPR.h_rs;
summary.BPRaff.irfAvg   = BPRaff.irf_Avg;   summary.BPRaff.bl_irfAvg   = BPRaff.bl_irf_Avg; 
summary.BPRaff.irfSem   = BPRaff.irf_Sem;   summary.BPRaff.bl_irfSem   = BPRaff.bl_irf_Sem; 
summary.BPRaff.irfStd   = BPRaff.irf_Std;   summary.BPRaff.bl_irfStd   = BPRaff.bl_irf_Std;
summary.BPRcor.irfAvg   = BPRcor.irf_Avg;   summary.BPRcor.bl_irfAvg   = BPRcor.bl_irf_Avg; 
summary.BPRcor.irfSem   = BPRcor.irf_Sem;   summary.BPRcor.bl_irfSem   = BPRcor.bl_irf_Sem; 
summary.BPRcor.irfStd   = BPRcor.irf_Std;   summary.BPRcor.bl_irfStd   = BPRcor.bl_irf_Std;
summary.nBlinks         = BPR.nAllB; % all blink only
summary.recTime         = nSamples/sr; % whole recording time
summary.blinkRate_const = BPR.brAllB; % = summary.nBlinks/summary.recTime;

if ~opts.SuppressFig
figure(301); clf; hold on;
h1 = shadedErrorBar(BPR.tVec, BPRaff.irf_Avg,  BPRaff.irf_Sem,  'k',1);
h2 = shadedErrorBar(BPR.tVec, BPRcor.irf_Avg,  BPRcor.irf_Sem,  'r',1);
% BPR-free period time course DATA
h3 = shadedErrorBar(BPR.tVec, BPRfree.irf_Avg, BPRfree.irf_Sem, 'b',1);
h4 = plot(BPR.tVec, BPR.h+BPRaff.irf_Avg(1), '--');
xlabel(xlabel1); ylabel(ylabel1);
% legend([h1.mainLine h2.mainLine, h3.mainLine], 'data','BPR corrected','BPR-free', 'Location', 'southeast');
legend([h1.mainLine h2.mainLine, h3.mainLine, h4], 'data','BPR corrected','BPR-free','h', 'Location', 'southeast');
if isFigSave, savefigpng(opts.path1, 'BPR correction (summary)', 'pre-post-free'); end
end

if ~opts.SuppressFig
figure(302); clf; hold on;
h1 = shadedErrorBar(BPR.tVec_rs, BPRaff.irf_rsAvg,  BPRaff.irf_rsSem,  'k',1);
h2 = shadedErrorBar(BPR.tVec_rs, BPRcor.irf_rsAvg,  BPRcor.irf_rsSem,  'r',1);
% BPR-free period time course DATA
h3 = shadedErrorBar(BPR.tVec_rs, BPRfree.irf_rsAvg, BPRfree.irf_rsSem, 'b',1);
h4 = plot(BPR.tVec_rs, BPR.h_rs+BPRaff.irf_rsAvg(1), '--');
xlabel(xlabel1); ylabel(ylabel1);
% legend([h1.mainLine h2.mainLine, h3.mainLine], 'data','BPR corrected','BPR-free', 'Location', 'southeast');
legend([h1.mainLine h2.mainLine, h3.mainLine, h4], 'data','BPR corrected','BPR-free','h', 'Location', 'southeast');
if isFigSave, savefigpng(opts.path1, 'BPR correction (summary) (resampled)', 'pre-post-free'); end
end


end

function irf = getIRF(ts, bOffIdx, range, esr, padWhat)
if isempty(padWhat) || strcmpi(padWhat,'NaN')
    padWhat = 'NaN'; % default
elseif strcmpi(padWhat,'Last')
    padWhat = 'Last';
else
    error('padding should be ''NaN'', or ''Last''')
end
    rangeVec = (range(1)*esr:range(2)*esr);
    irf = NaN(numel(rangeVec), numel(bOffIdx));
    for b = 1:numel(bOffIdx)
        idxVec = bOffIdx(b) + rangeVec;
        nIdx = diff(range)*esr+1;
        if idxVec(1)<1 % first trials, left padding
            idxVec2 = 1:idxVec(end);
            if strcmpi(padWhat,'NaN')
                irf_tmp = [NaN(nIdx-numel(idxVec2),1); ts(idxVec2)];
            elseif strcmpi(padWhat,'Last')
                irf_tmp = [ts(1)*ones(nIdx-numel(idxVec2),1); ts(idxVec2)];
            end
        elseif idxVec(end)>numel(ts) % last trals, right padding
            idxVec2 = idxVec(1):numel(ts);
            if strcmpi(padWhat,'NaN')
                irf_tmp = [ts(idxVec2); NaN(nIdx-numel(idxVec2),1)];
            elseif strcmpi(padWhat,'Last')
                irf_tmp = [ts(idxVec2); ts(end)*ones(nIdx-numel(idxVec2),1)];
            end
        else
            irf_tmp = ts(idxVec);
        end
        irf(:,b) = irf_tmp;
    end
end

function [opts] = parseinput(input, funcname, inputoffset)
% default values
opts = struct('Ismm','Auto','UseCano','Auto','LockIsl','Auto','CondName','','SavePath',[],'SuppressFig',false);
% opts.Color = {'k','b','m','r'};

names = {'Ismm','UseCano','LockIsl','FigSave','CondName','SavePath','SuppressFig'};
inputlen = length(input);
for i = 1:2:inputlen
    
    name = validatestring(input{i},names);
    
    value = input{i+1};
    switch name
        case 'Ismm'
            opts.Ismm = validatestring(value, {'Auto','A' 'No','N' 'Yes','Y'}, funcname, 'Ismm', i+1+inputoffset);
        case 'UseCano'
            opts.UseCano = validatestring(value, {'Auto','A' 'No','N' 'Yes','Y'}, funcname, 'UseCano', i+1+inputoffset);
        case 'LockIsl'
            opts.LockIsl = validatestring(value, {'Auto','A' 'No','N' 'Yes','Y'}, funcname, 'LockIsl', i+1+inputoffset);    
        case 'CondName'
            validateattributes(value,{'string','char'}, {}, funcname, 'CondName', i+1+inputoffset);
            opts.CondName = value;
        case 'SavePath'
            validateattributes(value,{'string','char'}, {}, funcname, 'SavePath', i+1+inputoffset);
            opts.SavePath = value;
        case 'SuppressFig'
            validateattributes(value,{'numeric','logical'},{'scalar', ...
                'integer'}, funcname, 'SuppressFig', i+1+inputoffset)            
            opts.SuppressFig = value;
    end
end
end

function [params, fval,  exitflag, output, history] = myproblem(func, initParams, minParams, maxParams, options)
% save output of every iteration
history = struct('p',[], 'funccount',[], 'fval',[], 'iteration',[], 'procedure',[]);
options = optimset('OutputFcn', @myoutput);

[params, fval,  exitflag, output] = fminsearchbnd(func, initParams, minParams, maxParams, options);
 
    function stop = myoutput(x,optimvalues,state)
        stop = false;
        if isequal(state,'iter')            
            history.p = [history.p x];
            history.funccount = [history.funccount optimvalues.funccount];
            history.fval = [history.fval optimvalues.fval];
            history.iteration = [history.iteration optimvalues.iteration];
            history.procedure{numel(history.funccount)} = optimvalues.procedure;
        end
    end
 
end