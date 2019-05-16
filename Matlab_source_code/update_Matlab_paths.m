function update_Matlab_paths
% get current directory
sTmp=what;
addpath(genpath(sTmp.path));
sMsg=sprintf('Directory : %s and all subdirectories added to Matlab path',sTmp.path);
disp(sMsg);
