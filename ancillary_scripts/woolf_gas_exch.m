function [f_tot, k0, kb, dele] = woolf_gas_exch(S,T,u10,slp,C,gas)

%-------------------------------------------------------------------------
% Calculate bubble-mediated air-sea gas exchange based on Woolf 1997
%
% USAGE:
%   [f_tot, k0, kb, dele] = woolf_gas_exch(S,T,u10,C,gas,u_factor);
%
% INPUT:
%   S = salinity (PSU)
%   T = temperature (C)
%   u10 = wind speed (m/s)
%   slp = sea level pressure (mbar)
%   C = gas concentration (mol/m3)
%   gas = one of 'o2', 'n2' or 'ar'
%   u_factor = multiplication factor to change Woolf/Thorpe 1991 ux values
%   (default = 2, as proposed by Vagle et al., 2020)
%
% OUTPUT:
%   f_tot = Total air-sea flux; positive out of ocean
%   k0 = non-bubble mediated piston velocity; m/d
%   kb = bubble mediated piston velocity; m/d
%   dele = bubble-induced supersaturation; %/100   
%
% R. Izett; rizett@eoas.ubc.ca
% UBC Oceanography
% Updated: March 2020
%-------------------------------------------------------------------------
    
% Calculate Ceq (units X)
    Ceq = gasmoleq(S,T,gas) .* slp./1013.25; %mol/m3
    if numel(Ceq) == 1
        Ceq = repmat(Ceq, size(u10));
    end
    if numel(C) == 1
        C = repmat(C, size(u10));
    end
    
% Calculate bubble-induced supersaturation, dele
    if strcmp(gas,'o2')
        ux = 9.0;
    elseif strcmp(gas,'n2')
        ux = 7.2;
    elseif strcmp(gas,'ar')
        ux = 9.6;
    else
        error('This parameterization is only for O2, N2 or Ar')
    end

    dele = 0.01 .* (u10 ./ ux) .^2; %Modified based on Woolf & Thorpe, 1991
    
% Calculate schmidt
%     sc = calc_schmidt(S,T,gas);
    [o2sc,n2sc,~,arsc] = sw_Schmidt_v3(S,T);
    if strcmp(gas,'o2')
        sc = o2sc;
    elseif strcmp(gas,'n2')
        sc = n2sc;
    elseif strcmp(gas,'ar')
        sc = arsc;
    end
    
% Calcualte air-side friction velocity (same units as u10)
    %Drag coefficient
    cd = nan(size(u10));
        cd(u10<11) = 0.0012;
        cd(u10>=11 & u10<=20) = (0.49 + 0.065 * u10(u10>=11 & u10<=20)) * 1e-3;
        cd(u10>20) = 0.0018;

    ua = sqrt(cd) .* u10; 
    
%Diffusive piston velocity [m/s]
    k0 = 1.57e-4 .* ua .* (sc./600).^(-0.5);
    k0 = k0 .* 3600 .* 24; %m/s --> m/d
    
%Bubble-mediated piston velocity [m/s]  
    betax = gasBunsen(S,T,gas) .* (T+273.15)./(273.15) .* 0.987; %Otswald solubility
    brakt = 1 + (14.*betax .* sc.^(-0.5)).^(-1/1.2);
    kb = 9.4e-3 .* u10.^3.41 .* betax .^ (-1) .* brakt .^(-1.2); %cm/hr
    kb = kb .* 24 ./ 100; %cm/hr --> m/d
    
% Total flux
    %mol/m2/d
    f_tot = -(k0 + kb) .* ((1 + dele) .* Ceq - C); %positive out of ocean
    
return   

