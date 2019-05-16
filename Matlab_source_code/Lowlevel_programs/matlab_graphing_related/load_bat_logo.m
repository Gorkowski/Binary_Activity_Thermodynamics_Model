function BAT_logo=load_bat_logo

name='load_bat_logo';



Lname=length(name);

mpath=mfilename('fullpath');
Lpath=length(mpath);
folderpath=mpath(1:Lpath-Lname);
database_path=[folderpath, 'BAT_logo.mat'];
load(database_path);

end