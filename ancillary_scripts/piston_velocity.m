function kw = piston_velocity( u10, sc, formulation, avg )
% function kw = piston_velocity( u10, sc, formulation, avg ) [m/d]
%
% compute gas transfer velocity as a function of wind speed
%
% u10 is the wind speed [m/s]
%
% sc is the schmidt number
%
% formulation (optional): 
%   n00 = Nightingale et al. 2005;
%   w92 = Wanninkhof 1992.
%   h06 = Ho et al., 2006
%   RI Additions:
%   w14 = Wanninkhof 2014 (updated / revisited)
%   bm16 = Butterworth & Miller, 2016 (eddy covariance / Southern Ocn.)
%   bm16b = Butterworth & Miller, 2016 forced through origin
%   lm86 = Liss & Merlivat 1986
%
% avg (optional):
%   0 = no averaging performed;
%   1 = mean kw returned based on linear combination of mean(u10) and
%       rms(u10).
%
% Matthew Long, 2007
% Edited by Robert Izett, 2020
%__________________________________________________________________________

% set default arguments
if nargin == 2, formulation = 'n00'; end
if nargin <= 3, avg = 0; end

% set coefficients
switch formulation
    case 'n00'
        a = [0; 0.3330; 0.2220; 0];
        sc_norm = 660;
        sc_exp = 1/2;
    case 'w92'
        a = [0; 0.0000; 0.3100; 0];
        sc_norm = 660;
        sc_exp = 1/2;
    case 'h06'
        a = [0; 0; 0.266; 0];
        sc_norm = 600;
        sc_exp = 1/2;
    case 'w14'
        a = [0; 0; 0.251; 0];
        sc_norm = 660;
        sc_exp = 1/2;
    case 'bm16'
        a = [1.3; 0; 0.245; 0];
        sc_norm = 660;
        sc_exp = 1/2;
    case 'bm16b' %intercept forced through origin
        a = [0; 0; 0.26; 0];
        sc_norm = 660;
        sc_exp = 1/2;
    case 'lm86'
        a = [0; 0.17; 0; 0]; %u10 < 3.6
        b = [-9.65; 2.85; 0; 0]; %u10 > 3.6 & u10 < 13
        c = [-49.3; 5.9; 0; 0]; %u10 > 13 
        sc_norm = 600;
        sc_exp = 2/3;
    otherwise
        error('Gas transfer formulation %s not supported.\n', formulation );
end

% perform averaging or not
if avg
    u10_sq = rms( u10 ) ^ 2;
    u10 = mean( u10 );
else
    u10_sq = u10 .^ 2;
end

% compute piston velocity (cm hr-1)
kw = (a(1) + a(2) .* u10 + a(3) .* u10_sq + a(4) .* u10.^3) .* ( ( sc ./ sc_norm ) .^ -(sc_exp) ); %units cm hr-1

if strcmp(formulation, 'lm86')
    if numel(sc) == 1
        sc = repmat(sc, size(u10));
    end
    ui = find(u10 >= 3.6 & u10 < 13);
        kw(ui) = (b(1) + b(2) .* u10(ui) + b(3) .* u10_sq(ui) + b(4) .* u10(ui).^3) .* ( ( sc(ui) ./ sc_norm ) .^ -0.5 ); %units cm hr-1
    ui = find(u10 >= 13);
        kw(ui) = (c(1) + c(2) .* u10(ui) + c(3) .* u10_sq(ui) + c(4) .* u10(ui).^3) .* ( ( sc(ui) ./ sc_norm ) .^ -0.5 ); %units cm hr-1
end
    
kw=kw*24/100;  %units m d-1




