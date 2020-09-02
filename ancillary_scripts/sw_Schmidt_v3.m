function [Schmidt_O2,Schmidt_N2,Schmidt_Ne,Schmidt_Ar,Schmidt_Kr,Schmidt_Xe,Schmidt_He,Schmidt_CFC11,Schmidt_CFC12,Schmidt_SF6] = sw_Schmidt_v3(S,T)

% sw_Schmidt_v3   Solubility of various gases in seawater
%=========================================================================
% sw_Schmidt_v3 8 Mar 2017
%          Author: Roberta C. Hamme (University of Victoria)
%
% USAGE:  [Schmidt_O2,Schmidt_N2,Schmidt_Ne,Schmidt_Ar,Schmidt_Kr,Schmidt_Xe,Schmidt_He,Schmidt_CFC11,Schmidt_CFC12,Schmidt_SF6] = sw_Schmidt_v3(S,T)
%
% DESCRIPTION:
%    This function calculates the Schmidt numbers in seawater of 
%    O2, N2, Ne, Ar, Kr, Xe, He, CFC-11, CFC-12, and SF6
%
% INPUT:  (if S and T are not singular they must have same dimensions)
%   S = salinity    [PSS-78]
%   T = temperature [degree C]
%
% OUTPUT:
%   Schmidt_O2 = Schmidt number of O2 [dimensionless]
%   Schmidt_N2 = Schmidt number of N2 [dimensionless]
%   Schmidt_Ne = Schmidt number of Ne [dimensionless]
%   Schmidt_Ar = Schmidt number of Ar [dimensionless]
%   Schmidt_Kr = Schmidt number of Kr [dimensionless]
%   Schmidt_Xe = Schmidt number of Xe [dimensionless]
%   Schmidt_He = Schmidt number of He [dimensionless]
%   Schmidt_CFC11 = Schmidt number of CFC11 [dimensionless]
%   Schmidt_CFC12 = Schmidt number of CFC12 [dimensionless]
%   Schmidt_SF6 = Schmidt number of SF6 [dimensionless]
% 
% AUTHOR:  Roberta Hamme (rhamme@uvic.ca)
%
% REFERENCE:
%    see help file for gas_diffusion.m and gas_diffusion_CFC.m for references
%       for the diffusion coefficients
%    Viscosity comes from M. H. Sharqawy, J. H. Lienhard V, and S. M. Zubair,
%       Desalination and Water Treatment, 16, 354-380, 2010.
%
% DISCLAIMER:
%    This software is provided "as is" without warranty of any kind.  
%=========================================================================

% CALLER: general purpose
% CALLEE: gas_diffusion, gas_diffusion_CFC, SW_Viscosity, sw_dens0

%----------------------
% Check input parameters
%----------------------
if nargin ~=2
   error('sw_Schmidt_v3.m: Must pass 2 parameters')
end %if

% CHECK S,T dimensions and verify consistent
[ms,ns] = size(S);
[mt,nt] = size(T);

  
% Check that T&S have the same shape or are singular
if ((ms~=mt) || (ns~=nt))
   error('sw_Schmidt_v3: S & T must have same dimensions')
end %if


% calculate viscosity
visc = SW_Viscosity(T,'C',S,'ppt') ./ sw_dens0(S,T) * 1e4;                        % kinematic viscosity at inputed temperature, salinity and 0 dbar in cm^2/s

[diff_Ne,diff_Ar,diff_Kr,diff_Xe,diff_O2,diff_N2,diff_He] = gas_diffusion(S,T);
[diff_CFC11,diff_CFC12,diff_SF6] = gas_diffusion_CFC(S,T);

Schmidt_Ne = visc./diff_Ne;
Schmidt_Ar = visc./diff_Ar;
Schmidt_N2 = visc./diff_N2;
Schmidt_O2 = visc./diff_O2;
Schmidt_He = visc./diff_He;
Schmidt_Kr = visc./diff_Kr;
Schmidt_Xe = visc./diff_Xe;
Schmidt_CFC11 = visc./diff_CFC11;
Schmidt_CFC12 = visc./diff_CFC12;
Schmidt_SF6 = visc./diff_SF6;