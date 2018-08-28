restoredefaultpath
filePath = mfilename('fullpath');
[sourceDir,~,~] = fileparts(filePath);

if (~exist(sourceDir,'dir'))
	error('  Error determining the paths. Exiting...');
end

warning off
rmpath(genpath(sourceDir));
warning on

addpath(sourceDir);
addpath(genpath(fullfile(sourceDir,'common')));
addpath(genpath(fullfile(sourceDir,'CWL_Segmentation_3D')));
addpath(fullfile(sourceDir,'Evaluation'));
addpath(fullfile(sourceDir,'FGT_Segmentation_3D'));
addpath(fullfile(sourceDir,'imageIO'));

matlabVersion=version('-release');