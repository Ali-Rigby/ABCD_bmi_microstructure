% Add unified BMI imaging pipeline paths.

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo_root, 'scripts'));
addpath(genpath(fullfile(repo_root, 'utilities', 'matlab_utils')));
fprintf('MATLAB paths: %s/scripts + utilities/matlab_utils\n', repo_root);
