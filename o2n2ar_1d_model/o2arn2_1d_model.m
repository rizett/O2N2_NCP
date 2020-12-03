%--------------------------------------------------------------------------
% This function performs 1D box model calcualtions of O2, Ar and N2.
% Refer to Izett & Tortell, 2020 for details
%
% INPUT:
%   met = structure containing atmospheric and ocean surface forcing data
%   profile = structure containing starting profile / initial conditions data
%   dom = structure containing information about model / run domain
%   expt = (optional) structure containing experimental forcing conditions
% 
% OUTPUT:
%   model = structure containing time-variable model output
%
% MODEL / SCRIPT AUTHOR
%   R. Izett
%   rizett@eoas.ubc.ca
%   UBC Oceanography
%   Last modified: August 2020

% REFERENCES:
% Izett, R. W. and Tortell, P. D. YYYY. delO2/N2' as a Tracer of Mixed Layer
%   Net Community Production: Theoretical Considerations and 
%   Proof-of-Concept.
%--------------------------------------------------------------------------

function model = o2arn2_1d_model(met, profile, dom, expt)

%----------------------
%--- Setup model domain
%----------------------
    %time step; days
        dt = dom.dt; 
    %total number of days
        tot = dom.days;     
    %start day; corresponding with input forcing data
        t0 = datenum(0,dom.strt_mon,dom.strt_day);
    %time array
        tt = t0:dt:t0+tot; %time 
        
    %Eddy diffusivity, m2/s        
        kz = dom.rkz_surf;
        if numel(kz)==1
            kz = repmat(dom.rkz_surf,size(tt)); 
        end        
        
    %Mixed layer depth, m 
        if isfield(dom,'mld')
            %if specified by user
            mld = dom.mld;
            if numel(mld)==1;
                mld = repmat(mld,size(tt));
            end
        else 
            %otherwise, estimate from starting T/S profile
            mld = depML(sw_dens(profile.s,profile.t,profile.z),profile.z,[],'N');
            mld = repmat(mld.a125,size(tt));
        end
        
        %Scale MLD, if specified by user, to examine how changes in MLD
        %affect gases
        if isfield(dom,'mld_scale')
            mld = mld .* dom.mld_scale;
        end
        
    %Biological O2 production (NCP), mmol O2/m2/d
        N = dom.biol; %mmol o2/m2/d
        if numel(N) == 1
            N = repmat(N,size(tt)); 
        end
        
    %N2-fixation, mmol N2/m2/d
        if isfield(dom,'n2fix')
            n2fix = dom.n2fix;
        else
            n2fix = 0;
        end
        
    %Depth of deep layer / vertical gradient integration depth
        if isfield(dom,'dz')
            dz = dom.dz;
        else
            dz = 25;
        end
        if numel(dz) == 1
            dz = repmat(dz,size(tt));
        end
        
    %Bubble mediated gas exchange scaling, beta
    %Use same value for large and small bubbles
        if isfield(dom,'wind_beta')
            gas_beta = dom.wind_beta;
        else
            gas_beta = 1;
        end
    
    %Interpolate environmental forcing to model time grid   
    %MLD salinity, PSU
        if ~isempty(profile)
            %Constant salinity from initial profile
            s_cycle = repmat(nanmean(profile.s(profile.z<=mld(1))), size(tt));
        else
            %Salinity time-series
            s_cycle = interp1(met.time,met.sss,tt);
        end
        
    %Temperature, deg-C
    %Wind speed, m/s
    %Density, kg/m3
    %Sea level pressure, mbar
    %Upwelling velocity, m/s
    %Ice cover, %/100
    if ~exist('expt','var') %if no experimental condition is applied
        %Wind speed
            u10_cycle = dom.wind_scale.*interp1(met.time,met.u10,tt); %wind speed
        
        %Surface Temperature
            t_cycle = interp1(met.time,met.sst,tt);
            
        %Density
            den_cycle = sw_dens(s_cycle,t_cycle,0);
        
        %SLP
            slp_cycle = interp1(met.time,met.slp,tt);
        
        %Upwelling velocity
            uw = interp1(met.time,met.pw,tt) ./ den_cycle;
            if all(dom.rkz_surf == 0)
                uw(:) = 0;
            end
        
        %Subsurface hydrography
        if ~isempty(profile)
            %Constant values from initial profile
            %Temperature
                t_deep = nanmean(profile.t((profile.z>nanmean(mld) & profile.z<=nanmean(mld)+nanmean(dz))));
                t_mld = nanmean(profile.t((profile.z<=nanmean(mld))));
                t_dif = t_deep - t_mld;
                t_deep = t_cycle + t_dif;
                clear t_mld 
            %Salinity
                s_deep = nanmean(profile.s(profile.z>nanmean(mld) & profile.z<=nanmean(mld)+nanmean(dz)));
                s_deep = repmat(s_deep, size(t_deep));
            %Density
                den_deep = sw_dens(s_deep,t_deep,mld+dz);
        else
            %From time-series observations
            %Temperature
                t_deep = interp1(met.deep_time,met.deep_t,tt);
            %Salinity
                s_deep = interp1(met.deep_time,met.deep_s,tt);
            %Density
                den_deep = interp1(met.deep_time,met.deep_d,tt);
        end
        
        %Ice cover
            if isfield(met,'ice_t')
                ice_cycle = interp1(met.ice_t+tt(1),met.sea_ice,tt);
            else
                ice_cycle = zeros(size(tt));
            end
        
    else %if experimental condition is specified
        %Wind speed
            u10_cycle = expt.u10;
            
        %Surface temperature
            t_cycle = expt.sst;
        
        %Density
            den_cycle = sw_dens(s_cycle,t_cycle,0);
            
        %SLP
            if isfield(expt,'slp')
                slp_cycle = expt.slp;
            else
                slp_cycle = repmat(1013.25,size(tt));
            end
            
        %Upwelling velocity
            uw = repmat(nanmean(met.pw)./nanmean(den_cycle),size(tt));
        
        %Subsurface hydrography
            %Temperature
            if isfield(expt,'tem_deep')
                t_deep = expt.tem_deep;
                
            elseif ~isempty(profile)
                t_deep = t_cycle - (nanmean(profile.t((profile.z<=nanmean(mld)))) - nanmean(profile.t((profile.z>nanmean(mld) & profile.z<=nanmean(mld)+nanmean(dz)))));
                
            else
                t_deep = interp1(met.deep_time,met.deep_t,tt);
            end
            
            %Salinity 
            if isfield(expt,'sal_deep')
                s_deep = expt.sal_deep;
                
            elseif ~isempty(profile)
                s_deep = nanmean(profile.s((profile.z>nanmean(mld) & profile.z<=nanmean(mld)+nanmean(dz))));
                s_deep = repmat(s_deep, size(t_deep));
                
            else
                s_deep = interp1(met.deep_time,met.deep_s,tt);
            end
            
        %Subsurface density
            den_deep = sw_dens(s_deep,t_deep,mld+dz);
            
        %Adjust mixing (on or off)
            if isfield(expt,'mix')
                kz = kz.*expt.mix;
                uw = uw.*expt.mix;
            end
        
        %Adjust upwelling (on or off)
            if isfield(expt,'upw')
                uw = uw.*expt.upw;
            end
        
        %Ice cover
            if isfield(expt,'ice')
                ice_cycle = expt.ice;
            else
                ice_cycle = zeros(size(tt));
            end
    end    
    uw = abs(uw);    
    
    %Initialize gas Variables
    %MLD concentrations, mmol/m3
        o2 = nan(size(tt)); 
            o2(1) = O2sol(s_cycle(1),t_cycle(1)).*den_cycle(1)./1000 .* dom.o2_sat_start .* slp_cycle(1)./1013.25;
        ar = nan(size(tt)); 
            ar(1) = Arsol(s_cycle(1),t_cycle(1)).*den_cycle(1)./1000 .* dom.sat_start .* slp_cycle(1)./1013.25;
        n2 = nan(size(tt)); 
            n2(1) = N2sol(s_cycle(1),t_cycle(1)).*den_cycle(1)./1000 .* dom.sat_start .* slp_cycle(1)./1013.25;
    
    %Subsurface, mmol/m3
        o2_deepeq = (O2sol(s_deep,t_deep).*den_deep./1000);
        ar_deepeq = (Arsol(s_deep,t_deep).*den_deep./1000);
        n2_deepeq = (N2sol(s_deep,t_deep).*den_deep./1000);
    
    %Adjust deep Ar if specified
        if isfield(dom,'sar_deep')
            ar_deep = ar_deepeq .* dom.sar_deep;
        else
            ar_deep = ar_deepeq;
        end
        
    %Adjust deep N2 (if specified)
        if isfield(dom, 'dn2ar_deep')
            n2_deep = n2_deepeq .* (1+dom.dn2ar_deep) .* (ar_deep./ar_deepeq);
        else
            n2_deep = n2_deepeq;
        end
        
    %Adjust deep O2 (if specified)
        if isfield(met,'deep_o2')
            o2_deep = interp1(met.deep_time,met.deep_o2,tt);        
        else
            o2_deep = o2_deepeq;
        end
        
    %Gas saturation, %
        o2sat = nan(size(tt)); 
            o2sat(1) = dom.o2_sat_start*100;
        arsat = nan(size(tt)); 
            arsat(1) = dom.sat_start*100;
        n2sat = nan(size(tt)); 
            n2sat(1) = dom.sat_start*100;
    
    %Delta-O2/Ar and Delta-O2/N2, %
        do2ar = nan(size(tt));  
            do2ar(1) = (dom.o2_sat_start./dom.sat_start - 1)*100;
        do2n2 = nan(size(tt)); 
            do2n2(1) = (dom.o2_sat_start./dom.sat_start - 1)*100;
    
%------------------------------------
%--- Initialize output data strucutre
%------------------------------------
    kk=1;
        model.domain       = dom;                           %domain information
        model.domain.dz    = dz;                            %subsurface integration depth, m
        model.time         = tt(kk)-tt(1);                  %time array, days since start
        model.mld          = mld(kk);                       %mixed layer depth, m
        model.mld_t        = t_cycle(kk);                   %surface temp., deg-C
        model.deep_t       = t_deep;                        %subsurface temp., deg-C
        model.mld_s        = s_cycle(kk);                   %surface salinity, PSU 
        model.deep_s       = s_deep;                        %subsurface salinity, PSU
        model.mld_dens     = den_cycle;                     %surface density, kg/m3
        model.deep_dens    = den_deep;                      %subsurface density, kg/m3
        model.u10          = u10_cycle(kk);                 %wind speed, m/s
        model.ice          = ice_cycle(kk);                 %ice cover, %/100
        model.slp          = slp_cycle(kk);                 %sea level pressure, mbar
        model.kz           = kz;                            %Eddy diffusivity array, m2/s
        model.upw          = uw;                            %Upwelling velocity array, m/s
        model.do2ar        = do2ar(kk);                     %del-O2/Ar
        model.do2n2        = do2n2(kk);                     %del-O2/N2        
        model.mld_o2sat    = o2sat(kk);                     %surface O2 saturation, %
        model.mld_arsat    = arsat(kk);                     %surface Ar saturation, %
        model.mld_n2sat    = n2sat(kk);                     %surface N2 saturation, %
        model.ceq_o2       = o2(kk)./dom.o2_sat_start;      %surface O2 equil conc., mmol/m3
        model.ceq_ar       = ar(kk)./dom.sat_start;         %surface Ar equil conc., mmol/m3
        model.ceq_n2       = n2(kk)./dom.sat_start;         %surface N2 equil conc., mmol/m3
        model.o2_deep      = o2_deep;                       %subsurface O2 conc., mmol/m3
        model.ar_deep      = ar_deep;                       %subsurface Ar conc., mmol/m3
        model.n2_deep      = n2_deep;                       %subsurface n2 conc., mmol/m3
        model.o2_deq       = 0;                             %Bubble-induced supersaturation, %/100
        model.ar_deq       = 0;
        model.n2_deq       = 0;
        %Fluxes, mmol/m3/d
        %O2
        model.o2_fd        = 0;                             %Diffusive air-sea flux of O2
        model.o2_fc        = 0;                             %Small bubble mediated air-sea flux of O2
        model.o2_fp        = 0;                             %Large bubble mediated air-sea flux of O2
        model.o2_fk        = 0;                             %Eddy diffusivity flux of O2
        model.o2_fw        = 0;                             %Upwelling flux of O2
        model.o2_fe        = 0;                             %Entrainment flux of O2
        model.o2_fN        = N(kk)./mld(kk);                %O2 NCP
        %Ar
        model.ar_fd        = 0;                             
        model.ar_fc        = 0;
        model.ar_fp        = 0;
        model.ar_fk        = 0;
        model.ar_fw        = 0;
        model.ar_fe        = 0;
        %N2
        model.n2_fd        = 0;
        model.n2_fc        = 0;
        model.n2_fp        = 0;
        model.n2_fk        = 0;
        model.n2_fw        = 0;
        model.n2_fe        = 0;
        model.n2fix        = n2fix./mld(kk);                %N2 fixation       
        
        %Subsurface dN2/Ar, %
        model.dn2ar_deep   = 100.*((n2_deep./n2_deepeq)./(ar_deep./ar_deepeq)-1);
        model.dn2ar_subs   = 100.*((n2_deep./(n2(1)./dom.sat_start))./(ar_deep./(ar(1)./dom.sat_start))-1); %relative to surface N2/Ar
        
        model

%-----------------------
%--- STEP THROUGH MODEL
%-----------------------
    fo2 = zeros(numel(tt),6);
    far = zeros(numel(tt),6);
    fn2 = zeros(numel(tt),6);
    deq = zeros(numel(tt),3);
    
    for kk = 2:numel(tt)
        
        %Calculate equilibrium concentrations from T/S/dens at current step
        %umol/kg * kg/m3 * 1mmol/1000umol --> %mmol/m3
            o2eq = O2sol(s_cycle(kk-1),t_cycle(kk-1)).*den_cycle(kk)./1000.*slp_cycle(kk-1)./1013.25;
            areq = Arsol(s_cycle(kk-1),t_cycle(kk-1)).*den_cycle(kk)./1000.*slp_cycle(kk-1)./1013.25;
            n2eq = N2sol(s_cycle(kk-1),t_cycle(kk-1)).*den_cycle(kk)./1000.*slp_cycle(kk-1)./1013.25;
                
        %Advection flux on gases ([C]/d)
            fo2(kk,1) = adv(o2(kk-1),o2_deep(kk-1),uw(kk-1).*3600.*24,dz(kk-1)); 
            far(kk,1) = adv(ar(kk-1),ar_deep(kk-1),uw(kk-1).*3600.*24,dz(kk-1)); 
            fn2(kk,1) = adv(n2(kk-1),n2_deep(kk-1),uw(kk-1).*3600.*24,dz(kk-1));
            
        %Eddy diffusivity flux ([C]/d)
            fo2(kk,2)=eddy_diff(o2(kk-1),o2_deep(kk-1),kz(kk-1),dz(kk-1)) ./ mld(kk-1);
            far(kk,2)=eddy_diff(ar(kk-1),ar_deep(kk-1),kz(kk-1),dz(kk-1)) ./ mld(kk-1);
            fn2(kk,2)=eddy_diff(n2(kk-1),n2_deep(kk-1),kz(kk-1),dz(kk-1)) ./ mld(kk-1);
            
        %Entrainment flux ([C]/d)
            dmld = (mld(kk)-mld(kk-1))./dt;
            if all(kz ==0) %i.e., no mixing flux
                dmld = 0;
            end
            if dmld > 0
                %m/d * mmol/m3 / m --> mmol/m3/d
                fo2(kk,6) = dmld .* (o2_deep(kk-1) - o2(kk-1)) ./ mld(kk); 
                far(kk,6) = dmld .* (ar_deep(kk-1) - ar(kk-1)) ./ mld(kk);
                fn2(kk,6) = dmld .* (n2_deep(kk-1) - n2(kk-1)) ./ mld(kk);
            else
                fo2(kk,6) = 0;
                far(kk,6) = 0;
                fn2(kk,6) = 0;
            end                
            
        %Gas exchange flux 
            %NOTE: User may modify below to add additional parameterizations
            if strcmp(dom.param,'l13') %Liang et al., 2013
                [fo2(kk,3),fo2(kk,4),fo2(kk,5),deq(kk,1)] = fas_L13(o2(kk-1)/1000,u10_cycle(kk-1),s_cycle(kk-1),t_cycle(kk-1),slp_cycle(kk-1)./1013.25,'o2',1);
                [far(kk,3),far(kk,4),far(kk,5),deq(kk,2)] = fas_L13(ar(kk-1)/1000,u10_cycle(kk-1),s_cycle(kk-1),t_cycle(kk-1),slp_cycle(kk-1)./1013.25,'ar',1);
                [fn2(kk,3),fn2(kk,4),fn2(kk,5),deq(kk,3)] = fas_L13(n2(kk-1)/1000,u10_cycle(kk-1),s_cycle(kk-1),t_cycle(kk-1),slp_cycle(kk-1)./1013.25,'n2',1);
            
            elseif strcmp(dom.param,'s09') %Stanley et al., 2009
                [fo2(kk,3),fo2(kk,4),fo2(kk,5),deq(kk,1)] = fas_S09(o2(kk-1)/1000,u10_cycle(kk-1),s_cycle(kk-1),t_cycle(kk-1),slp_cycle(kk-1)./1013.25,'o2',1);
                [far(kk,3),far(kk,4),far(kk,5),deq(kk,2)] = fas_S09(ar(kk-1)/1000,u10_cycle(kk-1),s_cycle(kk-1),t_cycle(kk-1),slp_cycle(kk-1)./1013.25,'ar',1);
                [fn2(kk,3),fn2(kk,4),fn2(kk,5),deq(kk,3)] = fas_S09(n2(kk-1)/1000,u10_cycle(kk-1),s_cycle(kk-1),t_cycle(kk-1),slp_cycle(kk-1)./1013.25,'n2',1);
                                
            elseif strcmp(dom.param,'v10') %Vagle et al., 2010
                [f_tot, ~, ~, deq(kk,1)] = vagle_gas_exch(s_cycle(kk-1),t_cycle(kk-1),u10_cycle(kk-1),slp_cycle(kk-1),o2(kk-1)/1000,'o2',1.7);
                    fo2(kk,3) = -f_tot./24./3600; %mol/m2/d --> mol/m2/s 
                    fo2(kk,4:5) = 0;
                    
                [f_tot, ~, ~, deq(kk,2)] = vagle_gas_exch(s_cycle(kk-1),t_cycle(kk-1),u10_cycle(kk-1),slp_cycle(kk-1),ar(kk-1)/1000,'ar',1.7);
                    far(kk,3) = -f_tot./24./3600; %mol/m2/d --> mol/m2/s 
                    far(kk,4:5) = 0;
                    
                [f_tot, ~, ~, deq(kk,3)] = vagle_gas_exch(s_cycle(kk-1),t_cycle(kk-1),u10_cycle(kk-1),slp_cycle(kk-1),n2(kk-1)/1000,'n2',1.7);
                    fn2(kk,3) = -f_tot./24./3600; %mol/m2/d --> mol/m2/s 
                    fn2(kk,4:5) = 0;
                
                clear f_tot
            end
                
        %Re-calc. gas exchange fluxes if there is sea ice present
            %NOTE: User may modify below to add additional parameterizations
            if ice_cycle(kk-1) > 0.01
                %Set bubble fluxes to 0 if there is ice present
                    fo2(kk,[4:5]) = [0 0];
                    far(kk,[4:5]) = [0 0];
                    fn2(kk,[4:5]) = [0 0];

                %Calc. diffusive flux (mol m-2 s-1)
                    fo2(kk,3)= fas_Fd(o2(kk-1)/1000,u10_cycle(kk-1),s_cycle(kk-1),t_cycle(kk-1),slp_cycle(kk-1)./1013.25,'o2','BM16',1) .* (1-ice_cycle(kk-1));
                    far(kk,3)= fas_Fd(ar(kk-1)/1000,u10_cycle(kk-1),s_cycle(kk-1),t_cycle(kk-1),slp_cycle(kk-1)./1013.25,'ar','BM16',1) .* (1-ice_cycle(kk-1));
                    fn2(kk,3)= fas_Fd(n2(kk-1)/1000,u10_cycle(kk-1),s_cycle(kk-1),t_cycle(kk-1),slp_cycle(kk-1)./1013.25,'n2','BM16',1) .* (1-ice_cycle(kk-1));
            end
            
        %Scale bubble fluxes
            fo2(kk,[4:5]) = fo2(kk,[4:5]) .* gas_beta;
            fn2(kk,[4:5]) = fn2(kk,[4:5]) .* gas_beta;
            far(kk,[4:5]) = far(kk,[4:5]) .* gas_beta;
            
%             error('re-calc del-eq using new beta-scaled terms')
        
        %Convert air-sea flux units
            fo2(kk,[3:5]) = fo2(kk,[3:5]).*1000.*3600.*24./mld(kk-1); %mol/m2/s --> mmol/m3/d
            far(kk,[3:5]) = far(kk,[3:5]).*1000.*3600.*24./mld(kk-1);
            fn2(kk,[3:5]) = fn2(kk,[3:5]).*1000.*3600.*24./mld(kk-1);            
            
        %Sum fluxes and calculate new gas concentrations [mmol/m3/d]*[d] = [mmol/m3]
            o2(kk) = o2(kk-1) + nansum(fo2(kk,[1:6])) .* dt + N(kk-1)./mld(kk-1).*dt;
            ar(kk) = ar(kk-1) + nansum(far(kk,[1:6])) .* dt; 
            n2(kk) = n2(kk-1) + nansum(fn2(kk,[1:6])) .* dt - n2fix./mld(kk-1).*dt; 
            
        %Calculate new saturation states and deltas
            o2sat(kk) = 100*o2(kk)./o2eq;
            arsat(kk) = 100*ar(kk)./areq;
            n2sat(kk) = 100*n2(kk)./n2eq;
            
            do2ar(kk) = (o2sat(kk)./arsat(kk) - 1)*100;
            do2n2(kk) = (o2sat(kk)./n2sat(kk) - 1)*100;
            
        %Write data to structure
            model.time(:,end+1)         = tt(kk)-tt(1);
            model.mld_t(:,end+1)        = t_cycle(kk);
            model.mld_s(:,end+1)		= s_cycle(kk);
            model.u10(end+1)            = u10_cycle(kk);
            model.slp(end+1)            = slp_cycle(kk);
            model.ice(end+1)            = ice_cycle(kk);
            model.do2ar(end+1)          = do2ar(kk);
            model.do2n2(end+1)          = do2n2(kk);
            model.mld(end+1)            = mld(kk);
            model.mld_o2sat(end+1)      = o2sat(kk);
            model.mld_arsat(end+1)      = arsat(kk);
            model.mld_n2sat(end+1)      = n2sat(kk);
            %Ceq [mmol/m3]
            model.ceq_o2(end+1)         = o2eq;
            model.ceq_ar(end+1)         = areq;
            model.ceq_n2(end+1)         = n2eq;
            %Fluxes [mmol/m3/d]
            model.o2_fd(end+1)          = fo2(kk,3);
            model.o2_fc(end+1)          = fo2(kk,4);
            model.o2_fp(end+1)          = fo2(kk,5);
            model.o2_fk(end+1)          = fo2(kk,2);
            model.o2_fw(end+1)          = fo2(kk,1);
            model.o2_fN(end+1)          = N(kk)./mld(kk);
            model.o2_fe(end+1)          = fo2(kk,6);
            model.ar_fd(end+1)          = far(kk,3);
            model.ar_fc(end+1)          = far(kk,4);
            model.ar_fp(end+1)          = far(kk,5);
            model.ar_fk(end+1)          = far(kk,2);
            model.ar_fw(end+1)          = far(kk,1);
            model.ar_fe(end+1)          = far(kk,6);
            model.n2_fd(end+1)          = fn2(kk,3);
            model.n2_fc(end+1)          = fn2(kk,4);
            model.n2_fp(end+1)          = fn2(kk,5);
            model.n2_fk(end+1)          = fn2(kk,2);
            model.n2_fw(end+1)          = fn2(kk,1);
            model.n2_fe(end+1)          = fn2(kk,6);
            model.n2fix(end+1)          = n2fix./mld(kk);
            model.o2_deq(end+1)         = deq(kk,1);
            model.ar_deq(end+1)         = deq(kk,2);
            model.n2_deq(end+1)         = deq(kk,3);
            model.dn2ar_subs(end+1)     = 100.*((n2_deep(kk)./n2eq)./(ar_deep(kk)./areq)-1);
                
    end
    
end

%Advection flux
function afl = adv(c_up,c_lo,w,dz)
    % Calculate vertical advection flux [units c / d]
    % w is m/d profile
    % c_up and c_lo are concentrations in upper and lower boxes
    
    %Vertical concentration gradient
    c_der = (c_lo  - c_up)  ./ (dz); %units [c] / m
   
    %Flux
    afl = c_der .* w; %units [c]/m * m/d == [c]/d
end

%Eddy diffuvity flux
function efl=eddy_diff(c_up,c_lo,kz,dz)
    % Calculate eddy diffusivity flux [units c * m / d]
    % kz is m2/s profile
    % c_up and c_lo are concentrations in upper and lower boxes

	%Vertical concentration gradient
    c_der = (c_lo  - c_up)  ./ (dz); %units [c] / m
   
    %Flux
    efl = c_der .* (kz*3600*24); %units [c]/m * m2/d == [c]*m/d

end


