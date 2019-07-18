function update_Matlab_paths
% get current directory
sTmp=what;
addpath(genpath(sTmp.path));
sMsg=sprintf('Directory : %s and all subdirectories added to Matlab path',sTmp.path);
disp(sMsg);

% check for development folder in directory above
% add all folders that contain "development_folder" in the name
cd ../
folder_list= string(ls);

dev_folder=(strfind(folder_list, 'development_folder'));
for i = 1:length(dev_folder)
    
    if cell2mat(dev_folder(i))==1
        full_name=fullfile(pwd,folder_list(i));
        
        addpath(genpath(full_name));
        sMsg=sprintf('Directory : %s and all subdirectories added to Matlab path',full_name);
        disp(sMsg);
    end
end

cd(sTmp.path)



