function [conc_O2] = O2sol(S,T)

% O2sol   Solubility of O2 in sea water
%=========================================================================
% O2sol Version 1.1 4/4/2005
%          Author: Roberta C. Hamme (Scripps Inst of Oceanography)
%
% USAGE:  concO2 = O2sol(S,T)
%
% DESCRIPTION:
%    Solubility (saturation) of oxygen (O2) in sea water
%    at 1-atm pressure of air including saturated water vapor
%
% INPUT:  (if S and T are not singular they must have same dimensions)
%   S = salinity    [PSS]
%   T = temperature [degree C]
%
% OUTPUT:
%   concO2 = solubility of O2  [umol/kg] 
% 
% AUTHOR:  Roberta Hamme (rhamme@ucsd.edu)
%
% REFERENCE:
%    Hernan E. Garcia and Louis I. Gordon, 1992.
%    "Oxygen solubility in seawater: Better fitting equations"
%    Limnology and Oceanography, 37, pp. 1307-1312.
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
   error('O2sol.m: Must pass 2 parameters')
end %if

% Check S,T dimensions and verify consistent
[ms,ns] = size(S);
[mt,nt] = size(T);

  
% Check that T&S have the same shape or are singular
if ((ms~=mt) | (ns~=nt)) & (ms+ns>2) & (mt+nt>2)
   error('O2sol: S & T must have same dimensions or be singular')
end %if

%------
% BEGIN
%------

% convert T to scaled temperature
Ts = log((298.15 - T)./(273.15 + T));

% constants from Table 1 of Garcia & Gordon for the fit to Benson and Krause (1984)
A0 = 5.80871; 
A1 = 3.20291;
A2 = 4.17887;
A3 = 5.10006;
A4 = -9.86643e-2;
A5 = 3.80369;
B0 = -7.01577e-3;
B1 = -7.70028e-3;
B2 = -1.13864e-2;
B3 = -9.51519e-3;
C0 = -2.75915e-7;

% Corrected Eqn (8) of Garcia and Gordon 1992
conc_O2 = exp(A0 + A1*Ts + A2*Ts.^2 + A3*Ts.^3 + A4*Ts.^4 + A5*Ts.^5 + S.*(B0 + B1*Ts + B2*Ts.^2 + B3*Ts.^3) + C0*S.^2);
