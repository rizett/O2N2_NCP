function [diff_CFC11,diff_CFC12,diff_SF6] = gas_diffusion_CFC(S,T)

% gas_diffusion_CFC   Diffusion coefficients for various CFC type gases in fresh/sea water
%=========================================================================
% gas_diffusion_CFC Version 1.0 16 July 2016
%          Author: Roberta C. Hamme (University of Victoria)
%
% USAGE:  [diff_CFC11,diff_CFC12,diff_SF6] = gas_diffusion_CFC(S,T)
%
% DESCRIPTION:
%    Diffusion coefficients of CFC-type gases in fresh/sea water
%
% INPUT:  (if S and T are not singular they must have same dimensions)
%   S = salinity    [PSS-78]
%   T = temperature [degree C]
%
% OUTPUT:
%   diff_CFC11 = diffusion coefficient of CFC11  [cm^2/s or 10^-4 m^2/s] 
%   diff_CFC12 = diffusion coefficient of CFC12  [cm^2/s or 10^-4 m^2/s] 
%   diff_SF6   = diffusion coefficient of SF6    [cm^2/s or 10^-4 m^2/s] 
% 
% AUTHOR:  Roberta Hamme (rhamme@uvic.ca)
%
% REFERENCE:
%    CFC-11 and CFC12 values from M. Zheng, W. J. De Bruyn, and E. S. Saltzman.
%       "Measurement of the dffusion coeffients of CFC-11 and CFC-12 in pure water and seawater"
%       J. Geophys. Res., 103(C1), 1375-1379.
%    SF6 values from D. B. King and E. S. Saltzman
%       "Measurement of the diffusion coefficient of sulfur hexafluoride in water"
%       J. Geophys. Res., 100(C4) 7083-7088.
%
% DISCLAIMER:
%    This software is provided "as is" without warranty of any kind.  
%=========================================================================

% CALLER: general purpose
% CALLEE: none


%----------------------
% Check input parameters
%----------------------
if nargin ~=2
   error('gas_diffusion_CFC.m: Must pass 2 parameters')
end %if

% CHECK S,T dimensions and verify consistent
[ms,ns] = size(S);
[mt,nt] = size(T);

  
% Check that T&S have the same shape or are singular
if ((ms~=mt) || (ns~=nt)) && (ms+ns>2) && (mt+nt>2)
   error('gas_diffusion_CFC.m: S & T must have same dimensions or be singular')
end %if

%------
% BEGIN
%------

% Calculate temperature in Kelvin
TK = T + 273.15;

% Apply Eyring equation and appropriate coefficients for each gas in freshwater
diff_CFC11 = 0.015 * exp(-18.1./(8.3144621e-3 * TK));         
diff_CFC12_0sal = 0.036 * exp(-20.1./(8.3144621e-3 * TK));         
diff_SF6 = 0.029 * exp(-19.3./(8.3144621e-3 * TK));         

% Apply salinity correction of 7.2% per 35 ppt salinity for CFC12 only
diff_CFC12 = diff_CFC12_0sal .* (1 - 0.072 * S / 35);


return
