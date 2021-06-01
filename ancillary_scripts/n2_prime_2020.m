%--------------------------------------------------------------------------
% This function performs calculates N2-prime from MLD budget calculations
% based on "historic" forcing datasets and observations of N2 saturation
% state in seawater. Also estimates the mixed layer O2 requilibration
% time-scale. 
% Calculations based on details in Izett & Tortell, in rev.
% 
% USAGE:
%   [n2pr,tro2] = n2_prime_2020(backdat,N,mix)
%
% INPUT:
%   backdat = structure containing "historic" forcing data during an
%       approximated O2 residence time, over which differences in Ar and N2
%       saturation state are estimates. Fields are:
%     .dt        = time-increment (days; 1x1) 
%     .time      = time arrary (days; 1xN) over which N2' / budget calculations are performed. 
%                 Incremented in dt time-steps.
%                 backdat.time(1) = 0 
%                 backdat.time(end) = upper-estimate for O2-requilibration time.
%     .mld_t     = mixed layer temperature array (deg-C; 1xN)
%     .u10       = wind speed at 10-m array (m/s; 1xN)
%     .slp       = sea level pressure array (mbar; 1xN)
%     .ice       = fraction of sea-ice coverage (%/100; 1xN); set to 0xN if
%                 non-polar .
%     .mld       = mixed layer depth array (m; 1xN); will often be assumed constant
%     .mld_s     = mixed layer salinity array (PSU; 1xN); will often be assumed constant
%       NOTE: for each of the arrays above, the FIRST value corresponds with
%         data at time 1 / the beginning of the O2 re-equilibration time-frame
%         (i.e. the most historic observation; at t-trO2). 
%         The LAST value corresponds with data at the time of in-situ 
%         observations (i.e. the most recent observation)
%     .n2sat     = N2 saturation state (%; e.g. 101; 1x1) at the time of observation
%     .param     = gas exchange parameterization. One of:
%                     'l13' = Liang et al., 2013
%                     'v10' = Vagle et al., 2010
%                     'w97' = Woolf 1997
%                     's09' = Stanley et al., 2009
%                     's07' = Sweeney et al., 2007
%                     'i11' = Ito et al., 2011
%                     'n16' = Nicholson et al., 2016   
%     .beta     = gas exchange scaling coefficient
% 
%   N = number of days of data overwhich N2 prime calculations are
%       performed; should be larger than O2 residence time (recommend 60-90
%       days)
% 
%   mix = structure containing mixing information to perform N2 prime
%       calculations. Fields are:
%     .kz        = mixing coefficients / dz (m/d; 1xN). NOTE: units m/d (not m2/s)
%     .ardeep    = subsurface Ar concentration (mmol/m3; 1xN)
%     .n2deep    = subsurface N2 concentration (mmol/m3; 1xN)
%   
% OUTPUT:
%   n2pr = estimated N2-prime (%)
%   tro2 = estimated O2 re-equilibration time (days)
%
% MODEL / SCRIPT AUTHOR
%   R. Izett
%   rizett@eoas.ubc.ca
%   UBC Oceanography
%   Last modified: Nov. 2020
% 
% REFERENCES:
% Izett, R. W. and Tortell, P. D. In review. delO2/N2' as a Tracer of 
%   Mixed Layer Net Community Production: Theoretical Considerations and 
%   Proof-of-Concept.
%--------------------------------------------------------------------------

function [n2pr,tro2] = n2_prime_2020(backdat,N,mix,err)

%--- Set model domain parameters
    dt = backdat.dt; %time-increment, days
    ti = N/dt; %number of time points
    
%--- Extract variables during specified / historic time interval
    new_t   = backdat.time(end-ti:end); %time, days
    mld     = backdat.mld(end-ti:end); %MLD, m
    mld_t   = backdat.mld_t(end-ti:end); %surface temp, deg-C
    mld_s   = backdat.mld_s(end-ti:end); %surface sal, PSU
    u10     = backdat.u10(end-ti:end); %wind speed, m/d
    slp     = backdat.slp(end-ti:end); %SLP, mbar
    n2sat   = backdat.n2sat; %observed/true n2 saturation, %
    kz      = mix.kz(end-ti:end); %mixing coefficient, m/d
    n2deep  = mix.n2deep(end-ti:end); %subsurface [];mmol/m3
    ardeep  = mix.ardeep(end-ti:end);   
    
    if ~isfield(backdat,'ice'); %ice cover, %/100
        ice = repmat(0,size(new_t));
    else
        ice = backdat.ice(end-ti:end);
    end
    ice(ice==1) = 0.99;
    
%--- Adjust parameters by error if input variable 'err' provided
    if exist('err','var')
        if isfield(err,'kz')
            e = rand(size(kz)) .* err.kz;
                kz = kz + (e - nanmean(e));
                clear e
        end
        if isfield(err,'deep')
            e = rand(size(n2deep)) .* err.deep;
                e = e-nanmean(e);
                n2deep = n2deep .*(1+e);
                clear e
        end
        if isfield(err,'mld')
            e = rand(size(mld)) .* err.mld;
                mld = mld + (e - nanmean(e));
                clear e
        end
        if isfield(err,'sst')
            e = rand(size(mld_t)) .* err.sst;
                mld_t = mld_t + (e - nanmean(e));
                clear e
        end
        if isfield(err,'u10')
            e = rand(size(u10)) .* err.u10;
                u10 = u10 + (e - nanmean(e));
                clear e
        end
        if isfield(err,'slp')
            e = rand(size(slp)) .* err.slp;
                slp = slp + (e - nanmean(e));
                clear e
        end
        if isfield(err,'beta')
            e = rand(1) .* 0.3;
            backdat.beta = backdat.beta + (e-.15+1);
        end
    end    
    
%--- Get N2 and Ar equilibrium concentrations over backdat time period, mmol/m3
    n2sol = N2sol(mld_s,mld_t).*sw_dens(mld_s,mld_t,0)./1000;
    arsol = Arsol(mld_s,mld_t).*sw_dens(mld_s,mld_t,0)./1000;
                    
%--- Calculate O2 re-equilibration time from weighted kT (combined
%--- diffusive and bubble-mediated piston velocity) and kz (mixing coefficient)
    %Piston velocities
        if strcmp(backdat.param,'l13') %Liang et al., 2013           
            [~, ~, ~, ~, ko2,~,ko2b] = fas_L13(0,u10,mld_s,mld_t,1,'O2',1);
            %Combine bubble and diffusive k
                ko2 = (ko2+ko2b.*backdat.beta).*3600.*24; %m/s --> m/d

        elseif strcmp(backdat.param,'v10') %Vagle et al., 2010
            [~, ko2, ko2b,o2_deq] = vagle_gas_exch(mld_s,mld_t,u10,1013.25,0,'o2',1.7);
            %Combine bubble and diffusive k
                ko2 = ko2 + ko2b.*backdat.beta; %m/d

        elseif strcmp(backdat.param,'w97') %Woolf, 1997
            [~, ko2, ko2b,o2_deq] = woolf_gas_exch(mld_s,mld_t,u10,1013.25,0,'o2');
            %Combine bubble and diffusive k
                ko2 = ko2 + ko2b.*backdat.beta; %m/d

        elseif strcmp(backdat.param,'s09') %Stanley et al, 2009
            [~, ~, ~, ~, ko2] = fas_S09(0,u10,mld_s,mld_t,1,'O2',1);            
            %Combine bubble and diffusive k
                ko2 = (ko2).*3600.*24; %m/s --> m/d

        elseif strcmp(backdat.param,'s07') %Sweeney 2007
            [~, ~, ~, ~, ko2] = fas_Sw07(0,u10,mld_s,mld_t,1,'O2',1);            
            %Combine bubble and diffusive k
                ko2 = (ko2).*3600.*24; %m/s --> m/d
            
        elseif strcmp(backdat.param,'i11') %Ito et al., 2011
            [~, ~, ~, ~, ko2] = fas_I11(0,u10,mld_s,mld_t,1,'O2',1);            
            %Combine bubble and diffusive k
                ko2 = (ko2).*3600.*24; %m/s --> m/d
            
        elseif strcmp(backdat.param,'n16') %Nicholson et al., 2016
            [~, ~, ~, ~, ko2] = fas_N16(0,u10,mld_s,mld_t,1,'O2',1);
            %Combine bubble and diffusive k
                ko2 = (ko2).*3600.*24; %m/s --> m/d
        else
            error('Please select a valid parameterization (l13, v10, s09, i11, n16 or w97) or modify the script')
        end

    %Adjust piston velocities if ice is present
        if any(ice ~= 0)
            [~,ko2i]=fas_Fd(0,u10,mld_s,mld_t,1,'O2','BM16',1);
                ko2(find(ice>0.01)) = ko2i(find(ice>0.01)) .* (1-ice(find(ice>0.01))) .*3600.*24;
        end

    %Weighted k and O2 re-equilibration time
        %Weighted O2 piston velocity
        [kwo2,wt] = kw_weighting(ko2(end-ti:end),dt,N,mld(end-ti:end));
        
        %Weighted MLD and mixing coefficient
        kwz     = nansum(wt.*mix.kz(end-ti:end))./nansum(wt);
        mldw    = nansum(wt.*mld(end-ti:end))./nansum(wt);
        
        %O2 residence time as MLD/k_weight
        tro2(1) = (mldw ./ kwo2);
        
        %Definition of O2 re-equilibration time from Izett & Tortell, 2020
        %use weighted mixing coefficient and mld (if time-variable values supplied)
        tro2(2) = -log(0.005)./(kwo2+kwz).*mldw; 
        
        if tro2(2) < 5; 
            tro2 = ceil(tro2(1)); 
        else
            tro2 = ceil(tro2(2));
        end
        if tro2 < 5; tro2 = 5; end
        
        ti = (tro2)/dt; %data indices over N2-prime time-scale (O2 re-equilibration time)
        
%--- N2 prime budget calculations over O2 re-equilibration time
        if tro2 > N
            ti = N/dt;
            tro2 = -tro2;
            warning('O2 re-equilibration time longer than supplied data record')
        end
        
    %Evaluate MLD model
        %Indices over which to evaluate Ar and N2 
            bi = length(new_t) - ti : length(new_t);    

        %Set initial conditions and respone variables
            respmat_n2 = nan(1,numel(bi)); 
            respmat_ar = nan(1,numel(bi));

            %starting concentrations to equilibrium (i.e. at t-tO2)
            respmat_n2(1) = n2sol(bi(1)) .* (slp(bi(1)))./1013.25;
            respmat_ar(1) = arsol(bi(1)) .* (slp(bi(1)))./1013.25;

        %Step through model
            %dC/dt * MLD = Fd + Fc + Fp + Fm
            for jj = 2:numel(bi)
                
                %Variables
                    here = bi(jj-1);
                    next = bi(jj);

                    h       = nanmean(mld([here]));
                    wind    = nanmean(u10([here]));
                    icec    = nanmean(ice([here]));
                    deep_n2 = nanmean(n2deep([here]));
                    deep_ar = nanmean(ardeep([here]));
                    pslp    = nanmean(slp([here]));
                    kappa   = nanmean(kz([here]));   
                    sal     = nanmean(mld_s([here]));
                    tem     = nanmean(mld_t([here]));
                
                    fn2     = zeros(1,4);
                    far     = zeros(1,4);
                    
                %Mixing flux, Fm (mmol/m2/d)
                    %recall, kappa is already in units of [m/d]
                    %         [m/d]            [mmol/m3]
                    far(1)  = kappa .* (deep_ar-respmat_ar(jj-1));
                    fn2(1)  = kappa .* (deep_n2-respmat_n2(jj-1));                    
                    
                %Air-sea exchange, Fd + Fc + Fp (mmol/m2/d)
                    if strcmp(backdat.param,'l13') %Liang et al., 2013                           
                        [far(2),far(3),far(4)] = fas_L13(respmat_ar(jj-1)/1000,wind,sal,tem,pslp./1013.25,'ar',1);
                        [fn2(2),fn2(3),fn2(4)] = fas_L13(respmat_n2(jj-1)/1000,wind,sal,tem,pslp./1013.25,'n2',1);

                        %Scale and convert units
                        far([3,4]) = far([3,4]) .* backdat.beta;
                            far([2:4]) = far([2:4]) *1000.*3600.*24; %mol/m2/s --> mmol/m2/d
                        fn2([3,4]) = fn2([3,4]) .* backdat.beta;
                            fn2([2:4]) = fn2([2:4]) *1000.*3600.*24; %mol/m2/s --> mmol/m2/d

                    elseif strcmp(backdat.param,'v10') %Vagle et al., 2010
                        [far(2)] = vagle_gas_exch(sal,tem,wind,pslp,respmat_ar(jj-1)/1000,'ar',1.7);
                        [fn2(2)] = vagle_gas_exch(sal,tem,wind,pslp,respmat_n2(jj-1)/1000,'n2',1.7);
                        
                        %Convert units
                        far(2) = -far(2).*1000; %mol/m2/d --> mmol/m2/d
                        fn2(2) = -fn2(2).*1000; %mol/m2/d --> mmol/m2/d

                    elseif strcmp(backdat.param,'w97') %Woolf, 1997
                        [far(2)] = woolf_gas_exch(sal,tem,wind,pslp,respmat_ar(jj-1)/1000,'ar');
                        [fn2(2)] = woolf_gas_exch(sal,tem,wind,pslp,respmat_n2(jj-1)/1000,'n2',1.7);
                        
                        %Convert units
                        far(2) = -far(2).*1000; %mol/m2/d --> mmol/m2/d
                        fn2(2) = -fn2(2).*1000; %mol/m2/d --> mmol/m2/d

                    elseif strcmp(backdat.param,'s09') %Stanley et al, 2009
                        [far(2),far(3),far(4)] = fas_S09(respmat_ar(jj-1)/1000,wind,sal,tem,pslp./1013.25,'ar',1);
                        [fn2(2),fn2(3),fn2(4)] = fas_S09(respmat_n2(jj-1)/1000,wind,sal,tem,pslp./1013.25,'n2',1);

                        %Scale and convert units
                        far([3,4]) = far([3,4]) .* backdat.beta;
                            far([2:4]) = far([2:4]) *1000.*3600.*24; %mol/m2/s --> mmol/m2/d
                        fn2([3,4]) = fn2([3,4]) .* backdat.beta;
                            fn2([2:4]) = fn2([2:4]) *1000.*3600.*24; %mol/m2/s --> mmol/m2/d
                        
                    elseif strcmp(backdat.param,'s07') %Sweeney 2007
                        [far(2),far(3),far(4)] = fas_Sw07(respmat_ar(jj-1)/1000,wind,sal,tem,pslp./1013.25,'ar',1);
                        [fn2(2),fn2(3),fn2(4)] = fas_Sw07(respmat_n2(jj-1)/1000,wind,sal,tem,pslp./1013.25,'n2',1);

                        %Scale and convert units
                        far([3,4]) = far([3,4]) .* backdat.beta;
                            far([2:4]) = far([2:4]) *1000.*3600.*24; %mol/m2/s --> mmol/m2/d
                        fn2([3,4]) = fn2([3,4]) .* backdat.beta;
                            fn2([2:4]) = fn2([2:4]) *1000.*3600.*24; %mol/m2/s --> mmol/m2/d

                    elseif strcmp(backdat.param,'i11') %Ito et al., 2011
                        [far(2),far(3),far(4)] = fas_I11(respmat_ar(jj-1)/1000,wind,sal,tem,pslp./1013.25,'ar',1);
                        [fn2(2),fn2(3),fn2(4)] = fas_I11(respmat_n2(jj-1)/1000,wind,sal,tem,pslp./1013.25,'n2',1);

                        %Scale and convert units
                        far([3,4]) = far([3,4]) .* backdat.beta;
                            far([2:4]) = far([2:4]) *1000.*3600.*24; %mol/m2/s --> mmol/m2/d
                        fn2([3,4]) = fn2([3,4]) .* backdat.beta;
                            fn2([2:4]) = fn2([2:4]) *1000.*3600.*24; %mol/m2/s --> mmol/m2/d
                            
                    elseif strcmp(backdat.param,'n16') %Nicholson et al., 2016
                        [far(2),far(3),far(4)] = fas_N16(respmat_ar(jj-1)/1000,wind,sal,tem,pslp./1013.25,'ar',1);
                        [fn2(2),fn2(3),fn2(4)] = fas_N16(respmat_n2(jj-1)/1000,wind,sal,tem,pslp./1013.25,'n2',1);

                        %Scale and convert units
                        far([3,4]) = far([3,4]) .* backdat.beta;
                            far([2:4]) = far([2:4]) *1000.*3600.*24; %mol/m2/s --> mmol/m2/d
                        fn2([3,4]) = fn2([3,4]) .* backdat.beta;
                            fn2([2:4]) = fn2([2:4]) *1000.*3600.*24; %mol/m2/s --> mmol/m2/d
                            
                    else
                        error('Please select a valid parameterization (l13, v10, s09, i11, n16 or w97) or modify the script')
                    end
                
                    %Re-calc. air-sea exch. if sea ice is present
                    if icec > 0.01
                        %Set bubble fluxes to 0 if there is ice present
                            far([3,4]) = [0 0];
                            fn2([3,4]) = [0 0];

                        %Calc. diffusive flux (mmol/m2/d)
                            far(2)= fas_Fd(respmat_ar(jj-1)/1000,wind,sal,tem,pslp./1013.25,'ar','BM16',1) .* (1-icec) *1000.*3600.*24; %mol/m2/s --> mmol/m2/d
                            fn2(2)= fas_Fd(respmat_n2(jj-1)/1000,wind,sal,tem,pslp./1013.25,'n2','BM16',1) .* (1-icec) *1000.*3600.*24; %mol/m2/s --> mmol/m2/d
                    end
                    
                %Apply errors to gas flux terms, if relevant
                    if exist('err','var')
                        if isfield(err,'fd')
                            e = rand(1) .* 0.3;
                                far(2) = far(2) .* (e-.15+1);
                                fn2(2) = fn2(2) .* (e-.15+1);
                                clear e
                        end
                        if isfield(err,'fp')
                            e = rand(1) .* 0.3;
                                far(3) = far(3) .* (e-.15+1);
                                fn2(3) = fn2(3) .* (e-.15+1);
                                clear e
                            clear e
                        end
                        if isfield(err,'fc')
                            e = rand(1) .* 0.3;
                                far(4) = far(4) .* (e-.15+1);
                                fn2(4) = fn2(4) .* (e-.15+1);
                                clear e
                        end
                    end
                    
                %Sum fluxes and determine new gas concentrations (mmol/m3)
                    respmat_ar(jj) = nansum(far) ./ h .* dt + respmat_ar(jj-1);
                    respmat_n2(jj) = nansum(fn2) ./ h .* dt + respmat_n2(jj-1);
                    
                clear far fn2 here next h wind icec deep_n2 deep_ar pslp kappa sal tem
            end
            
            %Convert to saturation state (%)
                respmat_ar = 100*respmat_ar ./ (arsol(bi).*slp(bi)./1013.25);
                respmat_n2 = 100*respmat_n2 ./ (n2sol(bi).*slp(bi)./1013.25);
        
            %Calculate the difference between modeled N2 and Ar on the last day of calculations
                dif_resp = respmat_n2(end) - respmat_ar(end);
                
            %N2-prime is observed N2sat minus modeled difference between N2 and Ar
                n2pr = n2sat - dif_resp; %[%]
                
        
return
