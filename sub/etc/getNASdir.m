function [NASdir] = getNASdir
currentDir = pwd;
tmpDir = dir('/Volumes/Data_CSNL*');
for i=1:numel(tmpDir)
    try
        dirStr = fullfile('/Volumes',tmpDir(i).name);
        % use cd for permission check 
        cd(dirStr)
        % if cd works, it means you are using the directory
        NASdir = dirStr;
        break;
    catch        
    end
end
cd(currentDir); % go back
end