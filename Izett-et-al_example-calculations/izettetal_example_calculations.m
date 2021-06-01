%% Example calculations from Izett et al. 
%
% All calculations are performed on a subset of the Arctic CAA/Shelf dataset
%
% Sections 1 and 2 can also be performed using the function:
%   [ncp,kw_o2,tro2,wt_t] = calc_o2n2_ncp(do2n2,S,T,mld,u10mat,Ndays,dt,mod);
%   (found in the O2N2_NCP_toolbox\o2n2_ncp folder)
% 
% Refer to Izett et al. for futher details:
%     Izett, R. W., Hamme, R. C., McNeil, C., Manning, C. C. M.,
%       Bourbonnais, A., and Tortell, P. D. ΔO2/N2' as a new tracer of marine net 
%       community production: Application and evaluation in the Subarctic Northeast 
%       Pacific and Canadian Arctic Ocean. Submitted. May 2021.

%% 1) Calculate O2 gas transfer velocity (kO2) from wind speed

clear all; close all; clc

% Some indicators / dummy variables
    arctic = 1; %1 = true; 0 = false

% CD to directory containing data
    cd('data\directory')
    
% Add folder containing ancillary scripts to the Matlab path
    addpath(genpath('ancillary\script\directory'));
    
% Load relevant data: underway, profile and historic environmental
    %underway:
    load('Izett-et-al_Arctic-Summer_underway.mat');
        uw = data;
    %profile:
    load('Izett-et-al_Arctic-Summer_profile.mat')
        prof = data;
    %Historic:
    load('Izett-et-al_Arctic-Summer_historic.mat')
        hist = data;
    clear data
        
% Get data incides for a subset of CAA/shelf data (data indices: di)
    di = find(uw.reg_indx==2 & uw.lat >= 75.5);
   
% Make historic matrices of relevant data
    tmat = hist.sst(di,:); %temp; deg-C
    smat = repmat(uw.sal(di),size(tmat,2),1)'; %sal; PSU
    windmat = hist.wind_speed(di,:); %u10, m/s
    zMLmat = hist.MLD(di,:); %MLD; m
    if arctic
        icemat = hist.ice_conc(di,:); %ice_conc; %
    else
        icemat = zero(size(tmat));
    end
    Pmat = hist.slp(di,:); %P_SLP; mbar
        
%Calculations 
    %1) O2 gas transfer velocity, kO2
    if arctic
        %Using Butterworth & Miller, 2016
        [~, k_o2] = fas_Fd(0,windmat,smat,tmat,Pmat,'o2','BM16'); 
        k_o2 = k_o2 .*3600.*24;
        k_o2 = k_o2 .*(1-icemat);
    else
        %Using Liang et al. 2013
        [~,~,~,~, k_o2] = fas_L13(0,windmat,smat,tmat,Pmat,'o2');
        k_o2 = k_o2 .*3600.*24;
    end
    
    % 2) Weighted kO2
        kw_o2 = nan(size(di)); %dummy variable for weighted k
        wt_t = 30; wt_t = repmat(wt_t,size(di)); %weighting time period, days
        for kk = 1:length(di)
            [kw_o2(kk)] = kw_weighting(k_o2(kk,end-wt_t(kk)*4:end), .25, wt_t(kk), zMLmat(kk,end-wt_t(kk)*4:end));
        end; clear kk 
        
% Plot results
    figure(1);
    subplot(3,1,1); hold on
        plot(uw.time(di),kw_o2,'k.')
        set(gca,'box','on')
        ylabel({'weighted k_{O2}';'[m/d]'})
    subplot(3,1,2); hold on
        plot(uw.time(di),100*nanmean(icemat,2),'k.')
        set(gca,'box','on')
        ylabel({'mean ice conc.';'over last 60 days';'[%]'})
    subplot(3,1,3); hold on
        plot(uw.time(di),nanmean(windmat,2),'k.')
        set(gca,'box','on')
        ylabel({'mean wind speed';'over last 60 days';'[m/s]'}) 
        xlabel('Julian Day')
    
%% 2) N2' calculations

clearvars -except uw di arctic kw_o2 hist prof
clc

N = 60; %min. number of days of data needed for calculation  

del_deep = 0.8; %Subsurface N2/Ar (N2/Ardeep) [%] value applied for BB; lookup in table S1

% Collect input data
    % time, lat, long, and sst/sss arrays from cruise / underway data
    tt = uw.time(di);
    dd = uw.dist(di);
    la = uw.lat(di);
    lo = uw.long(di);
    cr_t = uw.intake_SST(di);
    cr_s = uw.sal(di);
    cr_mld = hist.MLD(di,:); %MLD; m
    
    %historic matrices of u10, sst, slp, ice, sss from reanalysis products
    u10 = hist.wind_speed(di,:);
    slp = hist.slp(di,:);
    sst = hist.sst(di,:);
    if arctic
        ice = hist.ice_conc(di,:); %ice_conc; %
    else
        ice = zero(size(tmat));
    end
    sss = cr_s;
        sss = interp1(tt(~isnan(sss)),sss(~isnan(sss)),tt);
        
    %Observed N2, Ar data
    n2sat = uw.n2sat(di)*100;
    arsat = uw.o2sat(di)./(uw.do2ar(di)./100+1)*100;
    
    %backwards time array
    back_time = -N:.25:0; %time stamp of historic data
    
    %subsruface S & T: interpolate between stations
    %take average over 10 m below MLD
    mld = cr_mld(:,end)';
    for zz = 1:11;
        deep_T(zz,:) = griddata(prof.long,prof.lat,prof.depth,prof.temp,lo,la,mld+(zz)-1,'nearest');
        deep_S(zz,:) = griddata(prof.long,prof.lat,prof.depth,prof.sal,lo,la,mld+(zz)-1,'nearest');
    end
    deep_T=nanmean(deep_T);
    deep_S=nanmean(deep_S);
    
    %N2 and Ar at equil at base of ML
    n2sol = N2sol(deep_S,deep_T) .* sw_dens(deep_S,deep_T,mld) ./ 1000;
    arsol = Arsol(deep_S,deep_T) .* sw_dens(deep_S,deep_T,mld) ./ 1000;
    
    %mixing terms
    dn2ar_deep = repmat(del_deep, size(di));
    kmix = uw.uw_kz(di);
    
% Perform N2-prime calculation
    %Make some dummy variables
    n2pr_noMix = nan(size(tt));
    n2pr = nan(size(tt));
        
    %Go through all cruise / underway observation time points
    for kk = 1:10:numel(tt);
        %Don't do calculation if measured N2 is nan
        if isnan(n2sat(kk))
            continue
        end

        %Make an array of data going backwards in time
        backdat.dt = 0.5;    

        time            = tt(kk) + back_time;

        backdat.time    = time(1):backdat.dt:time(end);
        backdat.mld     = repmat(cr_mld(kk),size(backdat.time));
        backdat.mld_t   = interp1(time(~isnan(sst(kk,:))),sst(kk,~isnan(sst(kk,:))),backdat.time);
            if isnan(backdat.mld_t(end)); backdat.mld_t(end) = backdat.mld_t(end-1); end
        backdat.mld_s   = repmat(sss(kk),size(backdat.time));
        backdat.u10     = interp1(time,u10(kk,:),backdat.time);

        backdat.n2sat   = n2sat(kk);
        backdat.arsat   = arsat(kk);

        backdat.slp     = interp1(time,slp(kk,:),backdat.time);
        backdat.ice     = interp1(time,ice(kk,:),backdat.time);

        backdat.beta = 0.5; %wind speed scaling coefficient

        backdat.param = 'l13'; %gas transfer parameterization

        dz = 10; %integration depth; set to 10 m here, but should be set to (pycnocline depth - MLD)

        %mixing terms
        mix.kz = repmat(kmix(kk),size(backdat.time)) .* 3600 .* 24 ./ dz; %convert to m/d
        mix.ardeep = repmat(arsol(kk),size(backdat.time)); %deep Ar concentration at equil.
        mix.n2deep = repmat((dn2ar_deep(kk)./100 + 1) .* n2sol(kk),size(backdat.time)); %subsurface N2

        %Perform N2-rpime calculation
        [prime] = n2_prime_2020(backdat,N,mix); clc
        n2pr(kk) = prime;
        
        %Repeat by setting kz to zero
        mix.kz = repmat(0,size(backdat.time)); 
        [prime] = n2_prime_2020(backdat,N,mix); clc
        n2pr_noMix(kk) = prime;
        
        clear time prime mix backdat 
    end

    n2pr(n2pr==0)=nan;
    n2pr_noMix(n2pr_noMix==0)=nan;
    
%Plot results
    figure(2); 
    subplot(2,1,1); hold on
        plot(tt,n2sat,'k.--');
        plot(tt,arsat,'r.--');
        plot(tt,n2pr,'c.--')
        plot(tt,n2pr_noMix,'y.--')
        set(gca,'box','on')
        ylabel('% saturation')
        legend('Meas. N2','Meas. Ar','N2''','N2'' (without mixing)','location','best');
    
%% 3) Calculate NCP from O2/Ar (apply the same equations for O2/N2 and O2/N2')

clearvars -except uw di arctic kw_o2 hist prof n2pr n2pr_noMix

% units mmol O2 / m2 / d
    o2ar_ncp = uw.do2ar(di) ./ 100  .* ...
        kw_o2 .* ...
        (O2sol(uw.sal(di), uw.intake_SST(di)) .* sw_dens(uw.sal(di), uw.intake_SST(di),5) ./1000);

%Plot results
    figure(2); 
    subplot(2,1,2); hold on
        plot(uw.time(di),nanmoving_average(o2ar_ncp,60),'k.');
        set(gca,'box','on')
        ylabel({'O2/Ar-NCP';'[mmol O_2/m^2/d]'})
        xlabel('Julian day')
    