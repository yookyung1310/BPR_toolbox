function [pa, startIdx] = mergePupil_eyelink(eyeFilepath, P, varargin)
% NOTE: This function is for eyelink formate .edf (or converted .mat) data. Customize it if you use eyetracker from
% other suppliers.
% It converts .edf in eyeFilepath to.mat, loads the eyetracker data , and merge them for the next preprocessings.
% eyeFilepath: MxN cell. Each element has eyetracker file path string. If it ends with .mat or without extension, it
% just load and merge the data. If the string ends with .edf, it converts the file first, and then load and merge it.
% M is the number of runs and N is the number of subjects
% P: parameter variable. P.startword is the stamp word indicating the start of the run. 
% It should be predefined (e.g. P.startword = 'Recording Start';). Otherwise startIdx is 1, meaning it considers the run
% is started at the first time of the recording.


% decide how to pick an eye to be used, and get fieldname if the eyetracker info is in a structure
firstaxesinput = nargin-numel(varargin);
opts = parseinput(varargin, 'mergePupil', firstaxesinput);      

pa       = cell(P.nRuns,numel(P.sbjNames));
startIdx = zeros(P.nRuns,numel(P.sbjNames));

for sbj = 1:numel(P.sbjNames)
    for r = 1:P.nRuns % each run
        if ~isempty(eyeFilepath{r,sbj}) % get the inofrmation if only the path is not empty
            if strcmpi(eyeFilepath{r,sbj}(end-3:end),'.edf')
                edfdata = edfmex(eyeFilepath{r,sbj}); % ed = Eye Data
                eyeFilepath{r,sbj}(end-3:end) = '.mat'; % change the extension to load the converted
                save(eyeFilepath{r,sbj},'-struct','edfdata');
            end
            raw_ed = load(eyeFilepath{r,sbj}); % load data
            
            if isfield(opts,'fieldName') % if the info is inside a field.
                raw_ed = raw_ed.(opts.fieldName);
            end
            
            % 1. pick the eye
            % NOTE: We can pick an eye manually or by criteria, but it will use mean of both eyes at the moment.
            if opts.isPickAuto
                eyeIdx = raw_ed.RECORDINGS(1).eye;
            else
                eyeIdx = eyeMat(r,sbj);
            end
            switch eyeIdx
                case {1,2} % pick an eye
                    tmp = double(raw_ed.FSAMPLE.pa(eyeIdx,:));
                    pa{r,sbj} = tmp(:);% vectorize
                case 3 % mean of both eyes
                    tmp = mean(double(raw_ed.FSAMPLE.pa),1);
                    pa{r,sbj} = tmp(:);
            end
            % get the index when the run exactly started
            if isfield(P,'startWord')
                [sampleIdx, onsetTS, eventIdx] = getEventOnset(raw_ed, 'message', P.startword);
                startIdx(r,sbj) = sampleIdx;
            else
                startIdx(r,sbj) = 1;
            end
        end
    end
    fprintf('Subject %s (%dth) merge done\n', P.sbjNames{sbj}, sbj)
end

end

function [opts] = parseinput(input, funcname, inputoffset)
% default values
opts = struct('isPickAuto',true);
% opts.Color = {'k','b','m','r'};

names = {'isPickAuto','fieldName'};
inputlen = length(input);
for i = 1:2:inputlen    
    name = validatestring(input{i},names);
    
    value = input{i+1};
    switch name
        case 'isPickAuto'
            validateattributes(value,{'numeric','logical'},{'scalar', ...
                'integer'}, funcname, 'isPickAuto', i+1+inputoffset)
            opts.isPickAuto = value;
        case 'fieldName'
            validateattributes(value,{'string','char'}, {}, funcname, 'fieldName', i+1+inputoffset);
            opts.fieldName = value;
    end
end
end