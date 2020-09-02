function save_library(script,new_dir,lib_name);

%------------------------------------------------------------------------
% Save all functions in a script to a new toolbox location.
% 
% INPUT:
%   script = full directory link to script of interest
%       e.g. /link/to/script/script.m
%   new_dir = directory to which the toolbox / library will be saved
%       e.g. /new/directory/location
%   lib_name = name of new toolbox / library 
%       e.g. my_toolbox
% 
% OUTPUT:
%   check the new_dir!
% 
% Last updated: June 2020
% R. Izett, rizett@eoas.ubc.ca
% UBC Oceanography
%------------------------------------------------------------------------

%--- obtain list of all functions in script of interest
    [fList] = matlab.codetools.requiredFilesAndProducts(script);
    clc; fList'
        
%--- Create a new folder for the toolbox / library
    cdir = cd;
    cd(new_dir);
    mkdir(lib_name)
    cd(cdir);
    
    if ispc
        new_dir = [new_dir,'\',lib_name];
    else
        new_dir = [new_dir,'/',lib_name];
    end
    
%--- Copy all scripts to new location    
    for kk = 1:numel(fList);
        copyfile(fList{kk},new_dir)
    end    
    
return