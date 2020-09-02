function [f_tot, k0, kb, dele] = vagle_gas_exch(S,T,u10,slp,C,gas,u_factor)

%-------------------------------------------------------------------------
% Calculate bubble-mediated air-sea gas exchange based on Vagle et al.,
% 2010 and Steiner et al., 2007
%
% USAGE:
%   [f_tot, k0, kb, dele] = vagle_gas_exch(S,T,u10,C,gas,u_factor);
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
%
% R. Izett; rizett@eoas.ubc.ca
% UBC Oceanography
% Updated: Jan. 2020
%-------------------------------------------------------------------------

% Errors
    if ~strcmp(gas,'o2') & ~strcmp(gas,'n2') & ~strcmp(gas,'ar')
        error('Vagle et al., 2010 applies best for O2, N2 and Ar')
    end
%     if strcmp(gas,'ar')
%         warning('Vagle et al., 2010 applies best for O2 and N2')
%     end
    
% Calculate Ceq (units X)
    Ceq = gasmoleq(S,T,gas) .* slp./1013.25; %mol/m3
    if numel(Ceq) == 1
        Ceq = repmat(Ceq, size(u10));
    end
    if numel(C) == 1
        C = repmat(C, size(u10));
    end
    
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

% Calculate non-bubble mediated piston velocity based on Liss & Merlivat 86
    a = [0; 0.17; 0; 0]; %u10 < 3.6
    b = [-9.65; 2.85; 0; 0]; %u10 > 3.6 & u10 < 13
    
    if numel(sc) == 1
        sc = repmat(sc, size(u10));
    end
    
    k0 = (a(1) + a(2) .* u10 + a(3) .* u10.^2 + a(4) .* u10.^3) .* ( ( sc ./ 600 ) .^ -2/3 ); %units cm hr-1

    ui = find(u10 >= 3.6);
        k0(ui) = (b(1) + b(2) .* u10(ui) + b(3) .* u10(ui).^2 + b(4) .* u10(ui).^3) .* ( ( sc(ui) ./ 600 ) .^ -1/2 ); %units cm hr-1
    
    clear a b ui
    
    k0 = k0 .* 24 ./ 100; %units m/d
    
% Calculate bubble exchange transfer rate (Eq. 16 in Vagle et al., 2010)
    v0 = 2450; %cm/hr
    e = 14;
    n = 1.2;
    alph = gasBunsen(S,T,gas) .* (T+273.15)./(273.15) .* 0.987; %Otswald solubility
    %%%%WOOLF 1993    
    %convert L gas / (L soln * atm) to L / (L * bar)
    
    w = 3.84 .* 10^-6 .* (u10 .^ 3.41);
    
    kb = v0 .* w ./ alph .* (1 + (e .* alph .* (sc).^-0.5) .^ (-1./n)) .^ -n;
    
    kb = kb .* 24 ./ 100; %cm/hr --> m/d
    
% Calculate delta-e
    if strcmp(gas,'o2')
        ux = 9.0;
    elseif strcmp(gas,'n2')
        ux = 7.2;
    elseif strcmp(gas,'ar')
        ux = 9.6;
    else
        error('This parameterization is only for O2, N2 or Ar')
    end

    %adjust ux based on Steiner et al., 2007 and Vagle et al., 2010
        %Vagle: factor of 2
        %Steiner: factor of 1.7-2
        if nargin < 7
            u_factor = 2;
        elseif isempty(u_factor)
            u_factor = 2;
        end
        ux = ux .* u_factor;
    
    dele = 0.01 .* (u10 ./ ux) .^2; %Modified based on Woolf & Thorpe, 1991
       
% Calculate total flux (Eq. 17 in Vagle et al., 2010)    
    %mol/m2/d
    f_tot = -(k0 + kb) .* ((1 + dele) .* Ceq - C); %positive out of ocean
%     f_tot(u10<13) = -k0(u10<13) .* (Ceq(u10<13) - C(u10<13));
    
return   

