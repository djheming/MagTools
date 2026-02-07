function setup()

% Set workspace to include any necessary dependencies.

% Get the root of this MagTools installation.
magToolsPath = fileparts(mfilename('fullpath'));

% Add MagTools (and any internal /libs) to the path.
add_path_no_git(magToolsPath);

% Check if the dependency (BaseTools) is satisfied.
if exist('BaseTools', 'class')
    fprintf('MagTools and BaseTools are ready.\n');
else
    % Here, we haven't found BaseTools yet, but we can go look for it in a
    % sibling folder.
    parentDir = fileparts(magToolsPath);
    siblingBaseTools = fullfile(parentDir, 'BaseTools');
    if exist(siblingBaseTools, 'dir')
        add_path_no_git(siblingBaseTools);
    end
    if exist('BaseTools', 'class')
        fprintf('MagTools ready. BaseTools found in sibling folder.\n');
    else
        warning('BaseTools NOT FOUND.');
        fprintf('MagTools requires BaseTools to function.\n');
    end
end

end



function add_path_no_git(targetPath)

% Get the raw path string
rawPath = genpath(targetPath);

% Split into individual directory strings
% MATLAB uses ';' on Windows and ':' on Mac/Linux
pathList = strsplit(rawPath, pathsep);

% Filter out directories containing '.git' or 'resources' (common in Projects)
% We use 'filesep' to ensure we only match actual folder names
isBad = contains(pathList, [filesep '.git']) | ...
    contains(pathList, [filesep 'resources']) | ...
    cellfun(@isempty, pathList);
cleanPathList = pathList(~isBad);

% Join them back and add to path
if ~isempty(cleanPathList)
    addpath(strjoin(cleanPathList, pathsep));
end

end

