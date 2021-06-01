function [kw_wt,wt,fr] = kw_weighting(kw,dt,ndays,zml);

%-------------------------------------------------------------------------
% ABOUT:
% Calculate the weighted piston velocity for historic wind data.
% 
% INPUT:
% kw = vector of piston velocities [m/day], from a single location. NOTE:
% time is forward -- kw(1) is the oldest piston velocity; kw(end) is "today"
% dt = time interval [day]
% ndays = total number of days over which to weight the piston velocity
% zml = mixed layer depth (vector, same size as kw), or scalar
% 
% OUTPUT:
% kw_wt = weighted kw [m/day]
% wt = weighting 
% fr = fraction of mld ventilated in each step
% 
% NOTE: Kw weighting based on Ruer et al., 2007, with corrections from
% Teeter et al. (2014).
% 
% Script created by:
% R. Izett 
% Last updated: Oct. 2018
% Adapted from time_weigted_kw_kaiser by M. Long (2007).
%-------------------------------------------------------------------------

%--- adjust length of zml
    if numel(zml) == 1
        zml = repmat(zml,size(kw));
    end

%--- Trim input
    tot = ndays/dt;
    kw = kw(end-tot:end);
    zml = zml(end-tot:end);

%--- Calculate fraction of mld ventilated in each step
    fr = kw .* dt ./ zml; %piston velocity [m/d] * t step [d] / mld [m]
    fr(fr >= 1) = 0.9999;
    
%--- set weighting for most recent point to 1
    wt = nan(size(kw));
    wt(end) = 1;
    
%--- Calculate weights, going backwards in time, descending as a function of 
%--- the weight for the day following a particular day
    for ff = (ndays/dt): -1: 1
        wt(ff) = wt(ff+1) .* (1-fr(ff+1));
    end
        
%--- Calculate weighted kw
    prod = kw.*wt; %the product of kw * weight for each observation
    
    kw_wt = (sum(prod)) / (sum(wt)); 
        
end
