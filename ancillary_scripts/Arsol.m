function [conc_ar, Kh_ar,deriv,derivS] = Arsol(S,T)

% Arsol   Solubility of Ar in sea water
%=========================================================================
% Arsol Version 1.2 4/4/2005
%          Author: Roberta C. Hamme (Scripps Inst of Oceanography)
%
% USAGE:  concAr = Arsol(S,T)
%
% DESCRIPTION:
%    Solubility (saturation) of argon (Ar) in sea water
%    at 1-atm pressure of air including saturated water vapor
%
% INPUT:  (if S and T are not singular they must have same dimensions)
%   S = salinity    [PSS]
%   T = temperature [degree C]
%
% OUTPUT:
%   concAr = solubility of Ar  [umol/kg] 
% 
% AUTHOR:  Roberta Hamme (rhamme@ucsd.edu)
%
% REFERENCE:
%    Roberta Hamme and Steve Emerson, 2004.
%    "The solubility of neon, nitrogen and argon in distilled water and seawater."
%    Deep-Sea Research I, 51(11), p. 1517-1528.
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
   error('Arsol.m: Must pass 2 parameters')
end %if

% CHECK S,T dimensions and verify consistent
[ms,ns] = size(S);
[mt,nt] = size(T);

  
% Check that T&S have the same shape or are singular
if ((ms~=mt) | (ns~=nt)) & (ms+ns>2) & (mt+nt>2)
   error('Arsol: S & T must have same dimensions or be singular')
end %if

%------
% BEGIN
%------

% convert T to scaled temperature
temp_S = log((298.15 - T)./(273.15 + T));

% constants from Table 4 of Hamme and Emerson 2004
A0_ar = 2.79150;
A1_ar = 3.17609;
A2_ar = 4.13116;
A3_ar = 4.90379;
B0_ar = -6.96233e-3;
B1_ar = -7.66670e-3;
B2_ar = -1.16888e-2;

% Eqn (1) of Hamme and Emerson 2004
conc_ar = exp(A0_ar + A1_ar*temp_S + A2_ar*temp_S.^2 + A3_ar*temp_S.^3 + S.*(B0_ar + B1_ar*temp_S + B2_ar*temp_S.^2));
Kh_ar=conc_ar/1e6/(0.93/100);

% Calculate derivative
    deriv = conc_ar .* (A1_ar.*dT(T,1) + A2_ar.*dT(T,2) + A3_ar.*dT(T,3) + S.*(B1_ar.*dT(T,1) + B2_ar.*dT(T,2)));
    
    derivS = conc_ar .* (B0_ar + B1_ar*temp_S + B2_ar*temp_S.^2);
    
end

function der=dT(T,e)
    der = (e.*(log(298.15-T) - log(273.15+T)).^(e-1)).*(-1./(298.15-T)-1./(273.15+T));
end