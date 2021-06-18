function ed = process_preBPR(pa, startIdx, P)
% It preprocess the merged pupil area, pa before BPR correction. The preprocessing is as follow.
% (1) detect blinks, get rid of artifacts and interpolate the missings by blinks 
% (2) unit conversion: area (AU) to diameter (mm) [optional]
% (3) bandpass filterting                        

ed.startIdx   = startIdx; % just store the start index to the structure ed

ed.pd         = cell(P.nRuns,P.nSbjs);
ed.pd_tonic   = cell(P.nRuns,P.nSbjs);
ed.bOnIdxAll  = cell(P.nRuns,P.nSbjs);
ed.bOffIdxAll = cell(P.nRuns,P.nSbjs);
ed.bOffIdxInt = cell(P.nRuns,P.nSbjs);

nValidRunVec = sum(~cellfun(@isempty,pa));

for sbj = 1:numel(P.sbjNames)
    for r = 1:nValidRunVec(sbj) % each run
        if iscell(pa), y = pa{r,sbj}; else, y = pa(:,r,sbj); end
        
        % 2. treat missing by blink
        % [pupil, gazeX, gazeY, tp, onsetTS, varargout] = EyeFirstProcess_forPCaB(Eyedata,dp,conta_ms,varargin);
        [pupilTS2, ~, ~, bOnIdxAll, bOffIdxAll, bOffIdxInt, flagConta] = eyePrepro(y, P.conta, P.sr, P.diffCrit, P.wSize1);
        
        % 3. unit conversion: area (AU) to diameter (mm) [optional]
        if P.ismm            
            pupilTS3 = unitconversion(pupilTS2, P.convertCoeff, P.distance_prev, P.distance, P.base);
        else % just convert size (au) to diameter (au) assuming that it is a circle
            pupilTS3 = 2*sqrt(pupilTS2/pi);
        end
        
        % 4. bandpass filterting        
        if P.isLpass
            pupilTS4 = bpfilter(P.passCF(end), P.sr, P.passOrder, pupilTS3, 'low');
        else
            pupilTS4 = pupilTS3;
        end
        if P.isHpass
            pupilTS5 = bpfilter(P.passCF(1), P.sr, P.passOrder, pupilTS4, 'high');
            % remember the low frequency (tonic) pupil time course
            pupilTS_tonic = pupilTS4-pupilTS5;             
        end        
                        
        ed.pd{r,sbj}         = pupilTS5;
        ed.pd_tonic{r,sbj}   = pupilTS_tonic;
        ed.bOnIdxAll{r,sbj}  = bOnIdxAll;
        ed.bOffIdxAll{r,sbj} = bOffIdxAll;
        ed.bOffIdxInt{r,sbj} = bOffIdxInt;
        ed.flagConta{r,sbj}  = flagConta;    
        
        ed.nBlink(r,sbj,1) = numel(ed.bOffIdxAll{r,sbj});
        ed.nBlink(r,sbj,2) = numel(ed.bOffIdxInt{r,sbj});
    end
    fprintf('Subject %s (%dth) preBPR processing done\n', P.sbjNames{sbj}, sbj)
end
