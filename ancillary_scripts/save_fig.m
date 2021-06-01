function fname = save_fig(dir,name,type,res)

%---------------------------------------------------------------------------------
% Function to save a figure to a specified directory (dir) in .tif, .jpg, ,.png, 
% .fig and high-resolution .jpg formats.
%
% USAGE:
%   save_fig(dir,name)
% 
% INPUTS:
%   dir = directory where figure will be saved
%   name = figure name (without extension)
%   type = {'t','j','p','f',hr'} for tiff, jpeg, png, .fig and high-res.
%   res = desired resolution (commonly, 300, 500, 600, 800 ... dpi);
%   default = 500;
% 
% OUTPUTS:
%   fname = figure name and directory
%   figure saved to specified directory
% 
% Last updated: Aug. 2019
% R. Izett, rizett@eoas.ubc.ca
% UBC Oceanography
%---------------------------------------------------------------------------------

%--- set save type if unspecified
    if nargin < 3
        type = {'t'};
    end
    
%--- Create figure name:
    if ispc %check if PC or mac
        fname = [dir,'\',name]; %PC
    else
        fname = [dir,'/',name]; %mac
    end

%--- save figure
    for kk = 1:numel(type)
        if strcmp(type{kk},'t');      
            saveas(gcf,[fname,'.tif']);
        elseif strcmp(type{kk},'j')
            saveas(gcf,[fname,'.jpg']);
        elseif strcmp(type{kk},'p')
            saveas(gcf,[fname,'.png']);
        elseif  strcmp(type{kk},'f')            
            saveas(gcf,[fname,'.fig']);
        elseif  strcmp(type{kk},'hr')
            % set resolution
            if nargin < 4
                res = '-r500';
            else
                res = ['-r',num2str(res)];
            end
            print(gcf,[fname,'_high-res.jpg'],'-djpeg',res);
        elseif strcmp(type{kk},'hr_tif')
            if ispc
                % set resolution
                if nargin < 4
                    res = 300;
                end
                exportgraphics(gcf,[fname,'_high-res.tif'],'Resolution',res)
            else
                % set resolution
                if nargin < 4
                    res = '-r500';
                else
                    res = ['-r',num2str(res)];
                end
                print(gcf,[fname,'_high-res.tif'],'-dtiff',res);
            end
        end
            
    end
    
    