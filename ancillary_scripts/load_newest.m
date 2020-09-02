function [data,fname] = load_newest(fabbrev);

%--------------------------------------------------------------------------% 
% Function to load youngest file, in current directory, with matching file
% abbreviation (fabrev), if version control / file date stamping used.
%
% USAGE: fname = load_newest(fabrev);
%
% INPUT:
%   fabbrev = file abbreviation to search within current directory
% 
% OUTPUT:
%   data = loaded data
%   fname = name of file that is loaded
% 
% R. Izett (rizett@eoas.ubc.ca)
% UBC Oceanography
% Last modified: July 2019
%--------------------------------------------------------------------------

%--- Search current directory for file with specified abbreviation
    f = dir(['*',fabbrev,'*']);
    
%--- get rid of any ".", zip or image files
    rm = [];
    for kk = 1:numel(f)
        if strcmp(f(kk).name(1),'.');
            rm = [rm;kk];
        elseif strcmp(f(kk).name(end-2:end),'zip');
            rm = [rm;kk]; 
        elseif ~isempty(strfind(f(kk).name,'tif')) | ~isempty(strfind(f(kk).name,'jpg')) | ~isempty(strfind(f(kk).name,'jpeg')) | ~isempty(strfind(f(kk).name,'png'))
            rm = [rm;kk]; 
        end        
    end
    f(rm) = [];
    clear rm
    
%--- Get "youngest" file (most recent date)
    if numel(f) > 1
        for kk = 1:numel(f)
            fd(kk) = f(kk).datenum;
        end
        fi = find(fd == max(fd));
        if numel(fi)>1; fi = fi(1); end
        data = load(f(fi).name);
        fname = f(fi).name;
    else
        data = load(f.name);
        fname = f.name;
    end

display([fname,' loaded as:'])
data

return
