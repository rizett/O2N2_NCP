function zML = depML(rho,pres,ud,plt);

%-------------------------------------------------------------------------
% Function to calculate mixed layer depth.
% 
% usage: [z] = depML(rho,pres,ud,plt);
% 
% INPUT:
% rho, pres = profiles of density [kg/m3] and pressure [db]
% ud = string defining upward or downwards profile {'up','down'}; default is 'up'
% plt = string defining whether the user wants to plot density profile {'Y' or 'N'; default is N
% 
% OUTPUT:
% zML = structure containing mixed layer depths with anomalies of 0.125
% kg/m3 and 0.05 kg/m3
% 
% R. Izett, 2015
% UBC Oceanography
%-------------------------------------------------------------------------

%--- extract only profile in only one direction
    if ~exist('ud','var')
        ud = 'nan';
    end
    
    %remove nans
    no = find(isnan(rho)); 
    rho(no) = [];
    pres(no) = [];
    
    %find bottom of profile
    maxi = find(pres == max(pres)); %index of bottom depth
        
    if strcmp(ud,'up')==1
        up = length(pres):-1:maxi; %take index values for upward cast only
        
        pres = pres(up);
        rho = rho(up);
    elseif strcmp(ud,'down') == 1
        dwn = 1:maxi; %take index values for downward cast only
        
        pres = pres(dwn);
        rho = rho(dwn);
    end
    
rho(rho<1016)=nan; %remove any data with density less than 1016 kg/m3

%--- duplicate and smooth
    rho2 = rho; %rho = nanmoving_average(rho,2);    

%--- define surface density
    surf =find(pres<=10); %avg over top 10 metres 
%     surf =find(pres>10 & pres<=15); %avg over top 10-15 metres 
    if isempty(surf); surf =find(pres<=15); end %avg over top `5 metres 
    surf_rho = nanmean(rho2(surf));
    if isnan(surf_rho); surf = find(pres<=20); surf_rho = nanmean(rho2(surf)); end %avg over top 5 metres 
    
%--- find where density becomes greater than surface density plus the anomaly    
    if ~isempty(surf)
        a125 = find((rho > (surf_rho + 0.125) & pres>5),1,'first');
%         a125 = find(rho > (surf_rho + 0.125),1,'first'); %use 'first' to get only the first index value
        a05 = find((rho > (surf_rho + 0.05) & pres>5),1,'first');
%         a05 = find(rho > (surf_rho + 0.05),1,'first');
        a03 = find((rho > (surf_rho + 0.03) & pres>5),1,'first');
%         a05 = find(rho > (surf_rho + 0.05),1,'first');
        
        zML.a125 = pres(a125); 
        zML.a05 = pres(a05);
        zML.a03 = pres(a03);
        
    else
        a125 = 1:length(rho);
        a05 = 1:length(rho);
        a03 = 1:length(rho);
        zML.a125 = nan;
        zML.a05 = nan;
        zML.a03 = nan;
    end
    
%plot density profile with zML
if exist('plt','var') %if the plotting feature is specified
    if strcmp(plt,'Y')==1
        figure; hold on
        plot(rho2,pres,'k-'); 
        plot(rho,pres,'b-');
        p1=plot(rho2(a125),zML.a125,'or','markerfacecol','r');
        p2=plot(rho2(a05),zML.a05,'og','markerfacecol','g');
        p3=plot(rho2(a03),zML.a03,'om','markerfacecol','m');
        
        axis ij
        
        xlabel('\rho [kg m^-^3]','fontweight','bold')
        ylabel('P [dbar]','fontweight','bold')
        
        legend([p1 p2 p3],{'0.125 kg m^-^3 anomaly'; '0.05 kg m^-^3 anomaly'; '0.03 kg m^-^3 anomaly'},'location','southwest')
    end
end


end