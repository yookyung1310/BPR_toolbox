% Last update: 2021.06.27 by Kyung Yoo
% This code shows how to apply the BPR correction toolbox with example raw data from eyelink eyetracker.
% You may modify the part reading eyetracking files if you use other maker's eyetracker.
% It will show the preprocessing steps one-by-one, with two exmaplary subjects.
% Data used here is from Fixation task of Yoo, Ahn & Lee (2021) (IN REVISION)
% Please see Yoo, Ahn & Lee (2021) (IN REVISION) for detailed explanation, algorithm, and validation.

% This code requires 3 following Toolboxes
% 1) Signal Processing Toolbox
% 2) Statistics and Machine Learning Toolbox
% 3) Financial Toolbox

clear all
close all

% figure setting
defaultAxesFontSize  = 14;
defaultAxeslineWidth = 1;
defaultlineLineWidth = 2;
set(groot, 'defaultAxesFontSize', defaultAxesFontSize, 'defaultAxesFontweight', 'normal')
set(groot, 'defaultAxeslineWidth', defaultAxeslineWidth)
set(groot, 'defaultLinelineWidth', defaultlineLineWidth)

fileDir = mfilename('fullpath');
folderDir = fileparts(fileDir);

funcDir  = fullfile(folderDir, 'sub'); % functions to be used
addpath(genpath(funcDir));
dataDir  = fullfile(folderDir, 'data'); % input eyedata files ('.mat')
procdDir = fullfile(folderDir, 'processed'); % output eyedata file ('.mat')
if ~exist(procdDir,'dir'), mkdir(procdDir); end
figDir   = fullfile(folderDir, 'fig'); % output figure files ('fig','png', ...)
if ~exist(figDir,'dir'), mkdir(figDir); end
cd(folderDir)

sbjnames = dir(fullfile(dataDir, 'sbj*'));

%% Parameters & initialization
% experiment parameters
P.Nsbj     = length(sbjnames);
P.Nbg      = 4; % # backgrounds
P.Nrun     = 12; % # runs for each session (minimum/dark/medium/light BG condition has 3 runs per session, respectively)
P.Nsess    = 3; % # sessions
P.NbgRun   = P.Nrun / P.Nbg; % # runs who have the same background in each session

P.sr  = 500; % sampling rate (Hz)
P.dur = 140; % duration of single run (sec)

% analysis parameters
P.ismm      = 1; % convert pupil area (au) to pupil diameter (mm)

% Pupil has artifact around blink. It is conventionally processed as throwing away 150-200ms around blink.
% We use little bit more sophisticated method here. Blink offset artifact tend to be longer than onset artifact
% and somewhat vary across blinks  It will detect the blink offset artifact with pupil derivative blink-by-blink
% the artifact is assumed to be in 200 - 500ms range
% Note that this artifact preprocessing is not relavant to BPR correction. 
% P.conta(1) = artifact period before blink onset
% P.conta(2) = minimal artifact period after blink offset
% P.conta(3) = maximal artifact period after blink offset
P.conta     = [0.15 0.2 0.5]; % (sec)  

% NOTE: (2021.04.02) check overshooting by take long isolated period
P.wSize    = 3; % (sec) time window size to decide isolated blink. you don't have to consider about it if you will use 

P.sepCrit  = 0.2; % (sec) integrate separated blinks if it is 
P.durCrit  = [0.03 1]; % Assume that if it is real blink, duration can't be shorter than 30ms or longer than 1s
P.pickCrit = 25; % the number of necessary blinks to be picked (cirterion for pSbjs). the number of blinks should be larger than this across all BG

% Here, we will apply only low-pass filter (<4Hz) to remove electric noise
P.isLpass   = 1;
P.isHpass   = 0; 

P.eye = 1; % use left eye here. You may choose other eye.

%% preprocessing
% This code will load all the data one by one, sort eye data by
% background luminance, preprocess, and save it as merged data named 'ed'.
% The preprocessing steps are as follows.
% 1. unit conversion: area (AU) to diameter (mm) (optinoal)
% 2. detect and treat artifact by blink
% 3. bandpass filterting
%
% BPR correction will be applied after these preprocessing steps.

% check if the raw data is already processed. Load processed data if you have
loadDataFlag = input('Load a processed data (Y/N)? Type N if you run this code first time: ','s');
if upper(loadDataFlag) == 'N' % process from raw data
    for sbj = 1:P.Nsbj
        disp(sbjnames(sbj).name)
        sbjDir = [dataDir '/' sbjnames(sbj).name];
        % file that has experiment information, such as background level used at each run
        ParaFile = dir(fullfile(sbjDir, 'Para*')); 
        
        for iSess = 1:P.Nsess
            load(fullfile(sbjDir, ParaFile(iSess).name)); % load parameter info
            % Note that this .mat data was acquired by converting .edf file
            % with edfmex.m (https://github.com/iandol/opticka/blob/master/communication/edfmex.m). 
            % Proceed the conversion first if you have .edf file.
            % Here, we will read .mat file already converted from .edf file
            EyeFiles = dir(fullfile(sbjDir, ['EyeData_' Para.recTime{1}(1:4) '*']));            
            % get experiment start time by using the first sample of the first run
            load(fullfile(sbjDir, EyeFiles(1).name))
            fprintf('session%d has been loaded\n', iSess)
            for i = 1:P.Nbg
                % bring the runs from the same background level
                sameBgRuns = find(Para.runInfo == i);
                for r = (iSess-1)*P.NbgRun+1:iSess*P.NbgRun
                    
                    EyeFile_thisRun = EyeFiles(sameBgRuns(r-(iSess-1)*P.NbgRun)).name;
                    load(fullfile(sbjDir, EyeFile_thisRun)); % load P.eye data
                    
                    switch P.eye
                        case {1,2} % 1:left / 2:right
                            pupilTS = double(Eyedata.FSAMPLE.pa(P.eye,:));
                        case 3 % 3: mean of both eyes
                            pupilTS = mean(double(Eyedata.FSAMPLE.pa),1);
                    end
                    
                    startword = 'Run_onset';
                    [sampleIdx, onsetTS, eventIdx] = getEventOnset(Eyedata, 'message', startword);
                    pupilTS = pupilTS(sampleIdx:sampleIdx+P.sr*P.dur-1);
                    % 1. unit conversion: area (AU) to diameter (mm) (optinoal)
                    % sqrt(au) to mm coefficent measured in room649 setting (See mmconversion.m for detail)
                    if P.ismm
                        eye1000pluscoeff = 0.082734375;
                        pupilTS2 = unitconversion(pupilTS, eye1000pluscoeff, 47, 47, 44.5421);
                    end
                    ed.pdRawCat(:,i,r,sbj) = pupilTS2;
                    
                    % 2. detect and treat artifact by blink
                    [pupilTS3, bOnIdxAll, bOnIdxIsl, bOffIdxAll, bOffIdxIsl, bOffIdxArtAll, flagConta] = ...
                        eyePrepro(pupilTS2', P.sr, ...
                        'sepCrit',P.sepCrit, 'conta',P.conta, 'isolatedCrit',P.wSize, 'ismm',true, 'suppressFig',true);
                    
                    % 3. bandpass filterting
                    if P.isLpass && ~P.isHpass
                        fc = 4; order = 3;
                        bpStr = 'low';
                    elseif P.isLpass && P.isHpass
                        fc = [0.02 4]; order = 3;
                        bpStr = 'bandpass';
                    end
                    if P.isLpass || P.isHpass
                        pupilTS4 = bpfilter(fc, P.sr, 3, pupilTS3, bpStr);
                    else % if no bandpasss filtering at all.
                        pupilTS4 = pupilTS3;
                    end
                    
                    % concatenate runs across session
                    ed.pdCat(:,i,r,sbj) = pupilTS4;
                    ed.bOnIdxAllCat{i,r,sbj}     = bOnIdxAll;       % onset index (all blinks)
                    ed.bOnIdxIslCat{i,r,sbj}     = bOnIdxIsl;       % onset index (isolated blinks)
                    ed.bOffIdxAllCat{i,r,sbj}    = bOffIdxAll;      % offset index (all blinks)
                    ed.bOffIdxIslCat{i,r,sbj}    = bOffIdxIsl;      % offset index (isolated blinks)
                    ed.bOffIdxArtAllCat{i,r,sbj} = bOffIdxArtAll;   % blink offset artifact index
                    ed.flagContaCat{i,r,sbj}     = flagConta;       % mark true if the time point is confounded by blink or blink artifact
                    
                    % visualize the interim and final results.
                    % too many plots, so just show first run of 3rd BG.
                    if i == 3 && r == 1
                        fh = figure(100+sbj); hold on
                        fh.Position(3) = fh.Position(3)*2;
                        title(sprintf('pupil time course (%s, bg%d, run%d)', sbjnames(sbj).name, i, r))
                        x = (0:(numel(pupilTS)-1))/P.sr;
                        if P.ismm
                            plot(x, pupilTS2, 'LineWidth', 2)
                            str1 = 'mm converted (step1, optional)';
                            unitStr1 = 'diameter (mm)';
                        else
                            plot(x, pupilTS, 'LineWidth', 2)
                            str1 = 'raw';
                            unitStr1 = 'size (AU)';
                        end
                        plot(x, pupilTS3, 'LineWidth', 2)
                        plot(x, pupilTS4, 'LineWidth', 2)
                        xlabel('time (sec)')
                        ylabel(['pupil ' unitStr1])
                        legend({str1, 'treat arfitacts by blinks (step2)', 'low-pass filtered (step3, final)'}, ...
                            'location', 'southeast');
                        ylim([min(pupilTS4)-0.2 max(pupilTS4)+0.2])
                        
                        
                        
                        % show zoomed output
                        fh = figure(200+sbj); hold on
                        title(sprintf('zoomed pupil time course (%s, bg%d, run%d)', sbjnames(sbj).name, i, r))
                        x = (0:(numel(pupilTS)-1))/P.sr;
                        if P.ismm
                            plot(x, pupilTS2, 'LineWidth', 2)
                        else
                            plot(x, pupilTS, 'LineWidth', 2)
                        end
                        plot(x, pupilTS3, 'LineWidth', 2)
                        plot(x, pupilTS4, 'LineWidth', 2)
                        xlabel('time (sec)')
                        ylabel(['pupil ' unitStr1])
                        legend({str1, 'treat arfitacts by blinks (step2)', 'low-pass filtered (step3, final)'}, ...
                            'location', 'southeast');
                        ylim([min(pupilTS4)-0.2 max(pupilTS4)+0.2])
                        if sbj == 1
                            xlim([90 110]); ylim([5 5.9]);                                                        
                        elseif sbj == 2
                            xlim([68 74]); ylim([4.5 5.6]);       
                        end
                    end
                end
            end
        end
    end

%% Save processed data Exp1
% save except for few large matrix
clearvars Eyedata Stim ans
P.processedDate = datestr(now,'yyyymmdd_HHMM');
save(fullfile(procdDir,P.processedDate));

elseif upper(loadDataFlag) == 'Y' % load processed data
    procdData = input('Type the name of processed (merged) data to be loaded: ','s');
    fprintf('Load data. Directory: %s\n',fullfile(procdDir,procdData)); 
    try
    load(fullfile(procdDir,procdData));  
    catch
        error('file %s does not exist', procdData)
    end    
else
    error('type Y or N')
end

%% BPR correction
% Algorithm of BPR correction toolbox probabilistically infer, so need much time. 
% empirically 10-20 hours to correct 25-75 hours-long eye data if we used single CPU (3.4 GHz Intel Core i5)
disp('----------------- BPR correction starts -----------------')
disp('It will take appproximately 10-30 minutes...')
for sbj = 1:P.Nsbj
    for i = 1:P.Nbg        
        % collapse run dimension, but not background (= corneal flux density) dimension because
        % BPR shape is somewhat dependent with not only subject, but also background.                
        pre_pupil  = squeeze(ed.pdCat(:,i,:,sbj));
        bOffIdxIsl = squeeze(ed.bOffIdxIslCat(i,:,sbj));
        bOffIdxAll = squeeze(ed.bOffIdxAllCat(i,:,sbj));
        bOnIdxAll  = squeeze(ed.bOnIdxAllCat(i,:,sbj));
        
        % We recommned you to use a subject's own BPR profile, because the shape of BPR is slightly different across subject and
        % background corneal flux density. If it is not possible to get an average BPR profile from your data
        % (not enough isolated blinks, or other event is highly cross-correlated to blink, etc...), you can use canonical BPR from
        % fixation task experiment of ours. Set BPRUseCano 'Yes' to apply canonical BPR to all subjects, or 'Auto' to apply
        % specific subjects who does not have enough samples or has so different BPR profile estimated.
        % Yoo et al. (PLOS one, in revision) Please read on the paper for details.        
        
        % option of using canonical BPR shape. BPR shape is variable across
        % subjects, so we recommend not to use canonical BPR, but use it
        % only if only if estimating BPR is hard (e.g. too small number of blinks)
        edBC.UseCano = 'No';   % Use subject's own BPR shape (recommended)
%         edBC.UseCano = 'Yes';  % Use canonical BPR shape acquired from our fixation experiment. choose this option only if estimating BPR is hard.
        
        % option of using isolated blinks or all blinks to estimate BPR profile, h. 
        % isolated blinks are purer, but use all blinks if the number of isolated blink
        % is not enough. 'Auto' option will do the decision automatically
%         edBC.LockIsl = 'No'; % Use all blinks
%         edBC.LockIsl = 'Yes'; % Use isolated blinks only
        % Use isolated blinks only to estimate BPR profile if possible, but
        % use all blinks if there is not enough isolated blinks.(recommended) 
        edBC.LockIsl = 'Auto'; 
        
        % Lastly, it corrected BPR.
        [edBC.pdCat{i,sbj}, edBC.SF{i,sbj}, edBC.BPR{i,sbj}, edBC.BPRfree{i,sbj}, ...
            edBC.BPRaff{i,sbj}, edBC.BPRcor{i,sbj}, edBC.summary{i,sbj}] = ...
            bprcorrect(pre_pupil, bOffIdxIsl, bOffIdxAll, bOnIdxAll, P.sr, ...
            'UseCano',edBC.UseCano, 'LockIsl',edBC.LockIsl, 'SuppressFig',true, ...
            'SavePath',fullfile(figDir,'BPRresult'), 'CondName','BPRinfo');

        fprintf('sbj%d, bg%d BPR correction finished\n', sbj,i)
    end
end
AnalysisDate = datestr(now,'yymmddHHMM');
save(fullfile(procdDir,sprintf('edBC%s_%s',AnalysisDate,[edBC.UseCano edBC.LockIsl])),'-struct','edBC','-v7.3');

disp('----------------- BPR correction done -----------------')

% save the BPR-corrected data as previous 
save(fullfile(procdDir,P.processedDate),'-regexp','^(?!edBC$)\w');

%% plot the BPR correction result
% show the blink-locked time course to visualize if the shape of BPR is
% removed successfully.

%
str_sub = cell(P.Nsbj*P.Nbg);
for sbj = 1:P.Nsbj
    for i = 1:P.Nbg
        str_sub{4*(sbj-1)+i} = sprintf('%s bg%d', sbjnames(sbj).name, i);
    end
end

[fh1, axg1] = bprcorrSummary(edBC, str_sub);

% add title
for i = 1:numel(fh1)
    fh1(i).Position(3:4) = [1020 420];
end

% compare pre and post BPR correction result
% just show first run of 3rd BG as well
color_pre = 'k';
color_post = [0 158 115]/255; % bluish green

i = 3; r = 1;
for sbj = 1:P.Nsbj
    fh2 = figure(300+sbj); clf; hold on
    fh2.Position(3) = fh2.Position(3)*2;
    
    y_pre  = ed.pdCat(:,i,r,sbj);
    y_post = edBC.pdCat{i,sbj}(:,r);
    
    title(sprintf('pupil time course (%s, bg%d, run%d)', sbjnames(sbj).name, i, r))
    x = (0:(numel(pupilTS)-1))/P.sr;
    if P.ismm        
        unitStr1 = 'diameter (mm)';
    else
        unitStr1 = 'size (AU)';
    end
    l1 = plot(x, y_pre, 'LineWidth', 2, 'color', color_pre);
    l2 = plot(x, y_post, 'LineWidth', 2, 'color', color_post);
    
    xlabel('time (sec)')
    ylabel(['pupil ' unitStr1])    
    ylim([min(y_pre)-0.2 max(y_pre)+0.2])
    
    v1 = vline(ed.bOffIdxAllCat{i,r,sbj}/P.sr, 'k:');
    for v = 1:numel(v1)
        v1(v).LineWidth = 1;
        v1(v).LineWidth = 1;
    end
    
    legend([l1, l2, v1(1)], {'Uncorrected', 'Corrected (BPR toolbox)', 'blink offset'}, ...
        'location', 'southeast');
end


