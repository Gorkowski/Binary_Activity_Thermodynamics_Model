function BAT_logo=load_bat_logo
%%
% Created by Kyle Gorkowski [LAPTOP-A4QKFAC8] on 2019-May-16  1:35 PM
% Copyright 2019 Kyle Gorkowski 
% 
% Matlab function to load BAT logo image.
%   the BAT_logo.mat data file must be in the same directory as the
%   funciton load_bat_logo
%%

name='load_bat_logo';


Lname=length(name);

mpath=mfilename('fullpath'); % finds funciton path
Lpath=length(mpath);
folderpath=mpath(1:Lpath-Lname); % finds file path
database_path=[folderpath, 'BAT_logo.mat']; % makes full path to data file
load(database_path); % loads data file, containing BAT_logo image data

end