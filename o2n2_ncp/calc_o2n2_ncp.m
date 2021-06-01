function [ncp,kw_o2,tro2,wt_t] = calc_o2n2_ncp(do2n2,S,T,mld,u10mat,Ndays,dt,mod);

%-------------------------------------------------------------------------
% Calculate uncorrected NCP from underway O2/N2 (or O2/N2', O2/Ar) 
% measurements, following:
%       NCP = (dO2/N2') * [O2]eq * k_O2;
% where [O2]eq is the equilibrium saturation concentration of O2, and k_O2
% is the weighted piston velocity of O2.
% 
% USAGE: [ncp,kw_o2,tro2,wt_t] = calc_o2n2_ncp(do2n2,S,T,mld,u10mat,Ndays,dt,mod);
% 
% INPUT:
%     do2n2,S,T,mld
%               = arrays for delO2/N2 [decimal %], sea surface salinity
%                   [PSU], sea surface temperature [C] and mixed layer
%                   depth [m]. Note that all arrays should be the same 
%                   size [1xN].
%     u10mat    = wind speed (at 10m) matrix [m/s] of size [MxN] where M
%                   corresponds with the number of observations backwards 
%                   in time. M = Ndays / dt + 1;
%                   This matrix is used to calculate a weighted O2 piston velocity.
%           NOTE: u10mat(:,M) = wind speed corresponding with time of
%           underway observations. u10mat(:,1) = most "historic" wind
%           speed.
%     Ndays     = number of days over which to perform piston velocity
%                   weighting (recommend 30 or 60)
%     dt        = time increment, in days, of u10mat and mldmat
%     mod       = gas exchange model; one of:
%               'n00' = Nightingale et al. 2005;
%               'w92' = Wanninkhof 1992.
%               'h06' = Ho et al., 2006
%               'w14' = Wanninkhof 2014 (updated / revisited) (DEFAULT)
%               'bm16' = Butterworth & Miller, 2016 (eddy covariance / Southern Ocn.)
%               'bm16b' = Butterworth & Miller, 2016 forced through origin
%               'lm86' = Liss & Merlivat 1986
% 
% OUTPUT:
%     ncp    = mmol O2 / m2 / day
%     kw_o2  = weighted piston velocity of O2 [m/day]
%     tro2   = O2 residence time (mld/ko2) [days]
%     wt_t   = weighting time period [days]
% 
% MODEL / SCRIPT AUTHOR
%   R. Izett
%   rizett@eoas.ubc.ca
%   UBC Oceanography
%   Last modified: Aug. 2020
%-------------------------------------------------------------------------

%--- Default gas exchange model
    if nargin<8
        mod = 'w14';
    end

%--- Reshape all variables
    do2n2   = reshape(do2n2,length(do2n2),1);
    S       = reshape(S,length(do2n2),1);
    T       = reshape(T,length(do2n2),1);
    mld     = reshape(mld,length(do2n2),1);
    u10mat  = reshape(u10mat,length(do2n2),numel(u10mat)/length(do2n2));   
    
%--- Create backwards time variable
    back_t  = -Ndays:dt:0;
    lag     = 1:numel(back_t);

%--- Calculate O2 saturation concentration (mmol O2/m3)
    o2eq    = O2sol(S,T); %umol/kg
    o2eq    = o2eq .* sw_dens(S,T,5) ./ 1000; %mmol/m3 = umol/kg * kg/m3 * 1 mmol/1000 umol
    
%--- Calculate Schmidt number of O2
    %Create T, S and MLD matrices - assume all are constant back in time
    T       = repmat(T,1,length(lag)); 
    S       = repmat(S,1,length(lag)); 
    mldmat  = repmat(mld,1,length(lag));
    
    sc_o2   = nan(size(T));
    for kk = 1:length(do2n2) %for all time (i.e. underway) observations
        for ii = 1:length(lag) %for all lag (i.e. backwards) observations
            [sc_o2(kk,ii)] = sw_Schmidt_v3(S(kk,ii),T(kk,ii));
    end,end; clear kk ii

%--- Calculate piston velocity [m/d] matrix
    k_o2 = piston_velocity(u10mat, sc_o2, mod);

%--- Calculate Weighted k
    wt_t    = repmat(Ndays,size(do2n2));
    wt_o2   = nan(size(k_o2));
    fr      = nan(size(k_o2));
    kw_o2   = nan(size(do2n2));
    tro2    = nan(size(do2n2));
        
    for kk = 1:length(do2n2)
        [kw_o2(kk),wt_o2(kk,end-wt_t(kk)*(1/dt):end),fr(kk,end-wt_t(kk)*(1/dt):end)] = ...
            kw_weighting(k_o2(kk,end-wt_t(kk)*(1/dt):end), ...
                dt, wt_t(kk), mldmat(kk,end-wt_t(kk)*(1/dt):end));
           
        %Estimate O2 residence time
            wtmld    = nansum(wt_o2(kk,end-wt_t(kk)*(1/dt):end) .* mldmat(kk,end-wt_t(kk)*(1/dt):end)) ./ nansum(wt_o2(kk,end-wt_t(kk)*(1/dt):end));
            
            tro2(kk) = wtmld ./ kw_o2(kk);
        
        clear ti wtmld
    end
    
%--- Calculate NCP!
    ncp = do2n2 .* o2eq .* kw_o2;
    
return
