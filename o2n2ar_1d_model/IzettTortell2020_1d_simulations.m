%% README

%THIS SCRIPT CONTAINS CODE TO PERFORM 1D NUMERICAL SIMULATIONS DESCRIBED IN
%IZETT & TORTELL, 2020 (GBC). TWO EXAMPLES (for simulations ExIF-1c and
%real-OSP c) ARE PROVIDED, WITH INPUT FORCING FILES. 

%SIMULATION RESULTS AND FIGURES ARE SAVED TO USER-DEFINED LOCATIONS.

%ADJUST THE CODE BELOW TO PERFORM ADDITIONAL SIMULATIONS.

%SCRIPT AUTHOR:
%   R. Izett
%   rizett@eoas.ubc.ca
%   UBC Oceanography
%   Last modified: July 2020

%REFERENCES:
% Izett, R. W. and Tortell, P. D. 2020. ?O2/N2' as a Tracer of Mixed Layer
%   Net Community Production: Theoretical Considerations and 
%   Proof-of-Concept. Global Biogeochemical Cycles.

%NOTE: 
% Many of the functions used herewithin were wirtten by others, but are
% replicated here for ease. All functions are referenced accordingly.
% Kindly ackowledge the contributions of other script authors if using 
% these functions in subsequent work.

%% 1) Perform experimental simulation (e.g., ExIF-1c)

clear all; close all; clc

%-------------------------------------------------------------------------
% INPUT INFORMATION

% Identify files containing forcing data
    %Note: input should be full path directory with .mat extension
    %Atmospheric and surface ocean focing
    force_env = 'directory\to\forcing\data\2017-NEP-offshore_met-forcing_v20191031.mat';
    %Initial profile
    force_pro = 'directory\to\forcing\data\NEP-Summer-offshore_profile-ini_v20191031.mat';
    
% Identify directory where output results and figure will be saved
    res_save_dir = 'output\directory';
    fig_save_dir = 'output\directory';
    
% Specify simulation name
    sim_name = 'ExIF_1c';
    %Note: the output results and figure file save names will contain this extension
    %e.g. 'o2arn2_1d_ExIF_1_vYYYYMMDD.mat' or 'o2arn2_1d_ExIF_1_vYYYYMMDD.jpg'
    
% Setup simulation domain and conditions
    % Domain variables
        dom.dt              = 3/24;     %simulation time-increment (days)
        dom.dz              = 10;       %depth of subsurface box (meters below MLD); subsurface properties are averaged over the depth range MLD to MLD + dom.dz 
        dom.strt_mon        = 6;        %starting month of model run corresponding with forcing/input data
        dom.strt_day        = 1;        %starting day of model run (day of month)
        dom.days            = 180;      %number of days to run
        dom.rkz_surf        = 1e-4;     %constant surface diffusivity (m2/s); leave empty, [], or set to 0 if no mixing
        dom.biol            = 0;        %mixed layer NCP [mmol o2/m2/d]
        dom.dn2ar_deep      = 1.5/100;  %subsurface dN2/Ar (%/100)
        dom.wind_scale      = 1;        %wind speed scaling factor
        dom.o2_sat_start    = 1;        %starting O2 saturation state
        dom.sat_start       = 1;        %starting Ar and N2 saturation state
        dom.param           = 'l13';    %Air-sea exchange parameterization; choose from l13, s09, v10
        dom.wind_beta       = 0.5;      %Bubble-flux scaling coefficient between 0 and 1
        %OPTIONAL (uncomment)
        %dom.n2fix = N2-fixation rate (mmol/m2/d_
        %dom.sar_deep = subsurface Ar supersaturation anomaly (%/100)

    % Experimental conditions
        %set to expt = []; if running realistic simulations
        expt.u10_1      = 5;        %wind speed at start and end of run (m/s)
        expt.u10_2      = 15;       %wind speed after step change (m/s); set equal to u10_l if constant
        expt.sst_1      = 10;       %SST at start and end of run (deg-C)
        expt.sst_2      = 10;       %SST after step change (deg-C); set equal to sst_1 if constant
        expt.ice_1      = 0;        %ice percent at start and end of run (%/100)
        expt.ice_2      = 0;        %ice percent after step change (%/100); set equal to ice_1 if constant
        expt.mix        = 1;        %mixing OFF (0) or ON (1)
        expt.t1         = 0.4;      %percent (%/100) of total model run time when first experimental step change occurs
        expt.t2         = 0.7;      %percent (%/100) of total model run time when first experimental step change occurs
        
%-------------------------------------------------------------------------

%-----------------------%
%--- RUN MAIN SCRIPT ---%
%-----------------------%

%--- Load the input forcing and initial condition data files
    load(force_env);
    load(force_pro); 
       
%--- Set experimental conditions (if applicable)  
    if ~isempty(expt)
        %time variable
            expt.time = datenum(0,dom.strt_mon,dom.strt_day):dom.dt:datenum(0,dom.strt_mon,dom.strt_day)+dom.days;
            
        %Timing of step changes
            ti1 = 1:floor(length(expt.time).*expt.t1); 
            ti2 = ti1(end)+2/dom.dt:floor(length(expt.time).*expt.t2);
            ti3 = ti2(end)+2/dom.dt:length(expt.time);
            
        %Adjust wind speed (u10) during experimental time periods
            expt.u10 = repmat(expt.u10_1,size(expt.time)); %nominal u10
            expt.u10(ti2) = expt.u10_2; %u10 during step change
            %interpolate between step change and nominal condition
            expt.u10(ti1(end):ti2(1)) = interp1([expt.time(ti1(end)) expt.time(ti2(1))],[expt.u10_1; expt.u10_2], expt.time(ti1(end):ti2(1)));
            expt.u10(ti2(end):ti3(1)) = interp1([expt.time(ti2(end)) expt.time(ti3(1))],[expt.u10_2; expt.u10_1], expt.time(ti2(end):ti3(1)));
            
        %Adjust temperature (sst) during experimental time periods
            expt.sst = repmat(expt.sst_1,size(expt.time)); %nominal sst        
            expt.sst(ti2) = expt.sst_2; %sst during step change
            %interpolate between step change and nominal condition
            expt.sst(ti1(end):ti2(1)) = interp1([expt.time(ti1(end)) expt.time(ti2(1))],[expt.sst_1; expt.sst_2], expt.time(ti1(end):ti2(1)));
            expt.sst(ti2(end):ti3(1)) = interp1([expt.time(ti2(end)) expt.time(ti3(1))],[expt.sst_2; expt.sst_1], expt.time(ti2(end):ti3(1)));     
        
        %Adjust ice cover (ice) during experimental time periods
            expt.ice = repmat(expt.ice_1,size(expt.time)); %nominal ice
            expt.ice(ti2) = expt.ice_2; %ice during step change
            %interpolate between step change and nominal condition
            expt.ice(ti1(end):ti2(1)) = interp1([expt.time(ti1(end)) expt.time(ti2(1))],[expt.ice_1; expt.ice_2], expt.time(ti1(end):ti2(1)));
            expt.ice(ti2(end):ti3(1)) = interp1([expt.time(ti2(end)) expt.time(ti3(1))],[expt.ice_2; expt.ice_1], expt.time(ti2(end):ti3(1)));
            
        %Cleanup expt variable
            expt = rmfield(expt,{'u10_1';'u10_2';'sst_1';'sst_2';'ice_1';'ice_2'});           
    end
    
%--- Run Model
    if ~isempty(expt)
        model = o2arn2_1d_model(met, profile, dom, expt);
    else
        model = o2arn2_1d_model(met, profile, dom);
    end
        
%--- Output figure
    set(groot, 'defaultFigurePosition',[10 10 700 600]) %figure position
    figure;
    %MLD Gas saturation
    ax(1)=subplot(5,1,1); hold on
        plot(model.time,model.mld_o2sat,'b.-','markersize',4)
        plot(model.time,model.mld_arsat,'r.-','markersize',4)
        plot(model.time,model.mld_n2sat,'g.-','markersize',4)
        legend('O2','Ar','N2')
        title('Gas saturation [%]')
        set(gca,'xticklabel',[],'box','on');
        
    %MLD dN2/Ar
    ax(2)=subplot(5,1,2); hold on
        plot(model.time, 100*(model.mld_n2sat./model.mld_arsat-1),'k.-','markersize',4)
        plot(model.time, model.dn2ar_deep,'r--')
        legend('Surface','Deep')
        title('\DeltaN_2/Ar [%]')
        set(gca,'xticklabel',[],'box','on');

    %Wind speed and ice
    subplot(5,1,3); hold on
        [ax([3,4]),p1,p2] = plotyy(model.time,model.u10,model.time,model.ice);
        set(p1,'color','k','linewidth',1.5); set(ax(3),'ycolor','k')
        set(p2,'color','c','linewidth',1.5); set(ax(4),'ycolor','c')
        title('Wind speed & ice')
        ylabel(ax(3),'Wind speed [m/s]');
        ylabel(ax(4),'Ice cover [%]');
        set(ax(3),'xticklabel',[]);
        set(ax(4),'xticklabel',[]);
        
    %Surface and deep temperature
    ax(5)=subplot(5,1,4); hold on
        plot(model.time, model.mld_t,'k-','linewidth',1.5);
        plot(model.time, model.deep_t,'-','color',[.5 .5 .5],'linewidth',1.5);
        legend('Surface','Deep')
        title('Temperature [deg-C]')
        set(gca,'xticklabel',[],'box','on');
        
    %MLD and mixing coef.
    subplot(5,1,5); hold on
        [ax([6,7]),p1,p2] = plotyy(model.time,model.mld,model.time,model.kz);
        set(p1,'color','k','linewidth',1.5); set(ax(6),'ycolor','k')
        set(p2,'color','m','linewidth',1.5); set(ax(7),'ycolor','m')
        title('MLD & Eddy Diffusivity')
        ylabel(ax(6),'MLD [m]');
        ylabel(ax(7),'\kappa_Z [m2/s]');
        xlabel(ax(6),'Model time')
        
    linkaxes(ax,'x')
        
%--- Save model output and figure
    sn = ['o2arn2_1d_',sim_name];

    sname=save_version(res_save_dir,sn,model,'model'); %save output data
    clc
    display('Model output:')
    model
    display(['Model output data saved as ' sname])
    
    sname = save_fig(fig_save_dir,sn,{'j'}) 
    display(' ')
    display(['Model output figure saved as ' sname])
    
clearvars -except dom model profile met

%% 2) Perform realistic simulation (e.g., Real-OSP c)

clear all; close all; clc

%-------------------------------------------------------------------------
% INPUT INFORMATION

% Identify files containing forcing data
    %Note: input should be full path directory with .mat extension
    %Atmospheric and surface ocean focing
    force_env = 'directory\to\forcing\data\NEP_smooth1011_met-forcing_v20200402.mat';
    %Initial profile
    force_pro = 'directory\to\profile\data\NEP_smooth1011_2010_profile-ini_v20200402.mat';
    
% Identify directory where output results and figure will be saved
    res_save_dir = 'output\directory';
    fig_save_dir = 'output\directory';
    
% Specify simulation name
    sim_name = 'Real-OSPc';
    %Note: the output results and figure file save names will contain this extension
    %e.g. 'o2arn2_1d_ExIF_1_vYYYYMMDD.mat' or 'o2arn2_1d_ExIF_1_vYYYYMMDD.jpg'
    
% Setup simulation domain and conditions
    % Domain variables
        dom.dt              = 0.5;      %simulation time-increment (days)
        dom.dz              = 25;       %depth of subsurface box (meters below MLD); subsurface properties are averaged over the depth range MLD to MLD + dom.dz 
        dom.strt_mon        = 1;        %starting month of model run corresponding with forcing/input data
        dom.strt_day        = 1;        %starting day of model run (day of month)
        dom.days            = 500;      %number of days to run
        %modify kz below for time variability
            dom.rkz_surf        = [];       %constant surface diffusivity (m2/s); leave empty, [], or set to 0 if no mixing
        %modify biol below for time variability
            dom.biol            = 0;        %mixed layer NCP [mmol o2/m2/d]
        %modify dN2/Ar_deep below for time variability
            dom.dn2ar_deep      = 0.075;    %subsurface dN2/Ar (%/100)
            dom.wind_scale      = 1;        %wind speed scaling factor
            dom.o2_sat_start    = 1;        %starting O2 saturation state
            dom.sat_start       = 1;        %starting Ar and N2 saturation state
            dom.param           = 'l13';    %Air-sea exchange parameterization; choose from l13, s09, v10
        %OPTIONAL (uncomment)
        %dom.n2fix = N2-fixation rate (mmol/m2/d)
        %dom.sar_deep = subsurface Ar supersaturation anomaly (%/100)

    % Experimental conditions
        %set to expt = []; if running realistic simulations
        expt = [];        
%-------------------------------------------------------------------------

%-----------------------%
%--- RUN MAIN SCRIPT ---%
%-----------------------%

%--- Load the input forcing and initial condition data files
    load(force_env);
    load(force_pro); 
       
%--- Set experimental conditions (if applicable)  
    if ~isempty(expt)
        %time variable
            expt.time = datenum(0,dom.strt_mon,dom.strt_day):dom.dt:datenum(0,dom.strt_mon,dom.strt_day)+dom.days;
            
        %Timing of step changes
            ti1 = 1:floor(length(expt.time).*expt.t1); 
            ti2 = ti1(end)+2/dom.dt:floor(length(expt.time).*expt.t2);
            ti3 = ti2(end)+2/dom.dt:length(expt.time);
            
        %Adjust wind speed (u10) during experimental time periods
            expt.u10 = repmat(expt.u10_1,size(expt.time)); %nominal u10
            expt.u10(ti2) = expt.u10_2; %u10 during step change
            %interpolate between step change and nominal condition
            expt.u10(ti1(end):ti2(1)) = interp1([expt.time(ti1(end)) expt.time(ti2(1))],[expt.u10_1; expt.u10_2], expt.time(ti1(end):ti2(1)));
            expt.u10(ti2(end):ti3(1)) = interp1([expt.time(ti2(end)) expt.time(ti3(1))],[expt.u10_2; expt.u10_1], expt.time(ti2(end):ti3(1)));
            
        %Adjust temperature (sst) during experimental time periods
            expt.sst = repmat(expt.sst_1,size(expt.time)); %nominal sst        
            expt.sst(ti2) = expt.sst_2; %sst during step change
            %interpolate between step change and nominal condition
            expt.sst(ti1(end):ti2(1)) = interp1([expt.time(ti1(end)) expt.time(ti2(1))],[expt.sst_1; expt.sst_2], expt.time(ti1(end):ti2(1)));
            expt.sst(ti2(end):ti3(1)) = interp1([expt.time(ti2(end)) expt.time(ti3(1))],[expt.sst_2; expt.sst_1], expt.time(ti2(end):ti3(1)));     
        
        %Adjust ice cover (ice) during experimental time periods
            expt.ice = repmat(expt.ice_1,size(expt.time)); %nominal ice
            expt.ice(ti2) = expt.ice_2; %ice during step change
            %interpolate between step change and nominal condition
            expt.ice(ti1(end):ti2(1)) = interp1([expt.time(ti1(end)) expt.time(ti2(1))],[expt.ice_1; expt.ice_2], expt.time(ti1(end):ti2(1)));
            expt.ice(ti2(end):ti3(1)) = interp1([expt.time(ti2(end)) expt.time(ti3(1))],[expt.ice_2; expt.ice_1], expt.time(ti2(end):ti3(1)));
            
        %Cleanup expt variable
            expt = rmfield(expt,{'u10_1';'u10_2';'sst_1';'sst_2';'ice_1';'ice_2'});           
    end
    
%--- Modify time-variabilie forcing components
    %kz, eddy diffusivity
        %Values from Cronin et al., 2015
        t_kz = datenum(0,[1:24],0); %months over which kz data are extracted
        kz = [[3.7 4.7 5.1 4.5 3.4 2.5 1.7 1.4 1.3 1.4 1.9 2.6]-.5] * 1e-4; %Cronin et al. kz values
        kz = [kz,kz];    
        t0 = datenum(0,dom.strt_mon,dom.strt_day);
        t_model = t0:dom.dt:t0+dom.days; %time array in model
        dom.rkz_surf = interp1(t_kz,kz,t_model,'spline');
        clear t_kz kz t0
        
    %MLD
        %Based on observations from OSP (see within 'met' structure)
        t_mld = datenum(0,met.mld_mon,1);
        dom.mld = interp1(t_mld,met.mld,t_model,'spline'); 
        clear t_mld
        
    %NCP (can set to 0 if you're not worried about modelling O2 production)
        %Values from Fassbender et al., 2016
        t_fas = datenum(0,[1:24],0); %months over which NCP data are extracted
        n = [-9,-2,4,18,25,8,12,14,10,0,-10,-15]*1.4; %NCP values from Fassbender et al.
        n = [n,n];    
        dom.biol = interp1(t_fas,n,t_model,'spline');
        clear t_fas n
        
    %Deep dN2/Ar and deep dAr
        %Values from Hamme et al., 2019 @ OSP
        %dn2/ar
        t_hamme = datenum(0,[2,6,8,14,18,20],0);  %months over which gas data are extracted
        deep = [.5 .25 0 .5 .25 0]/100; %gas data from Hamme et al., 2019
        dom.dn2ar_deep = interp1(t_hamme,deep,t_model,'pchip');
        dom.dn2ar_deep = nanmoving_average(dom.dn2ar_deep,60);
        
        %dar
        t_hamme = datenum(0,[1,2,6,8,13,14,18,20],0);
        deep = [0 0 1 1 0 0 1 1]/100+1;
        dom.sar_deep = interp1(t_hamme,deep,t_model,'pchip'); %deep ar saturation
        dom.sar_deep = nanmoving_average(dom.sar_deep,60);
        clear t_hamme deep t_model 
        
%--- Run Model
    if ~isempty(expt)
        model = o2arn2_1d_model(met, profile, dom, expt);
    else
        model = o2arn2_1d_model(met, profile, dom);
    end
        
%--- Output figure
    set(groot, 'defaultFigurePosition',[10 10 700 600]) %figure position
    figure;
    %MLD Gas saturation
    ax(1)=subplot(5,1,1); hold on
        plot(model.time,model.mld_o2sat,'b.-','markersize',4)
        plot(model.time,model.mld_arsat,'r.-','markersize',4)
        plot(model.time,model.mld_n2sat,'g.-','markersize',4)
        legend('O2','Ar','N2')
        title('Gas saturation [%]')
        set(gca,'xticklabel',[],'box','on');
        
    %MLD dN2/Ar
    ax(2)=subplot(5,1,2); hold on
        plot(model.time, 100*(model.mld_n2sat./model.mld_arsat-1),'k.-','markersize',4)
        plot(model.time, model.dn2ar_deep,'r--')
        legend('Surface','Deep')
        title('\DeltaN_2/Ar [%]')
        set(gca,'xticklabel',[],'box','on');

    %Wind speed and ice
    subplot(5,1,3); hold on
        [ax([3,4]),p1,p2] = plotyy(model.time,model.u10,model.time,model.ice);
        set(p1,'color','k','linewidth',1.5); set(ax(3),'ycolor','k')
        set(p2,'color','c','linewidth',1.5); set(ax(4),'ycolor','c')
        title('Wind speed & ice')
        ylabel(ax(3),'Wind speed [m/s]');
        ylabel(ax(4),'Ice cover [%]');
        set(ax(3),'xticklabel',[]);
        set(ax(4),'xticklabel',[]);
        
    %Surface and deep temperature
    ax(5)=subplot(5,1,4); hold on
        plot(model.time, model.mld_t,'k-','linewidth',1.5);
        plot(model.time, model.deep_t,'-','color',[.5 .5 .5],'linewidth',1.5);
        legend('Surface','Deep')
        title('Temperature [deg-C]')
        set(gca,'xticklabel',[],'box','on');
        
    %MLD and mixing coef.
    subplot(5,1,5); hold on
        [ax([6,7]),p1,p2] = plotyy(model.time,model.mld,model.time,model.kz);
        set(p1,'color','k','linewidth',1.5); set(ax(6),'ycolor','k')
        set(p2,'color','m','linewidth',1.5); set(ax(7),'ycolor','m')
        title('MLD & Eddy Diffusivity')
        ylabel(ax(6),'MLD [m]');
        ylabel(ax(7),'\kappa_Z [m2/s]');
        xlabel(ax(6),'Model time')
        
    linkaxes(ax,'x')
        
%--- Save model output and figure
    sn = [sim_name,'_o2arn2_1d-results'];

    sname=save_version(res_save_dir,sn,model,'model'); %save output data
    clc
    display('Model output:')
    model
    display(['Model output data saved as ' sname])
    
    sname = save_fig(fig_save_dir,sn,{'j'}) 
    display(' ')
    display(['Model output figure saved as ' sname])
    
clearvars -except dom model profile met

%% 3) Calculate N2-prime (N2') on simulations

clear all; close all; clc

%-------------------------------------------------------------------------
% INPUT INFORMATION

%Identify directory where output files are saved
    mod_dir = 'mode\data\dir';
    
% Identify directory where output figure will be saved
    fig_save_dir = 'output\directory';    
    
%Identify simulations for which N2' calculations will be performed
    %These should match the simulation names set above
    sim = {'ExIF_1c';'Real-OSPc'};

%Specify number of days of data needed for calculations
    %i.e. estimate an upper limit of O2-residence time
    %can be a single number, or one number per simulation run
    N = [60,100];
    
%Specify the time increment (in days) over which N2' calculations will be
%discretized
    dt = 0.5;
%-------------------------------------------------------------------------

%-----------------------%
%--- RUN MAIN SCRIPT ---%
%-----------------------%
    cd(mod_dir)
    if numel(N) == 1; N = repmat(N,size(sim)); end

%--- Go through each simulation
    for kk = 1:numel(sim)
        %load model output
        [data,fname] = load_newest(sim{kk}); clc
            model = data.model;
            clear data

        %N2 prime calculations
            %Initialize N2 prime variables
            model.n2pr = nan(size(model.time)); 
            model.n2pr_noMix = nan(size(model.time));
            model.o2_tres = nan(size(model.time));
            model.o2_tres_noMix = nan(size(model.time));

            %Go through all model observations 
            for tt = 1:numel(model.time)
                if tt <= N(kk)/model.domain.dt + 5
                    continue
                else

                %Standard N2-prime w/ variable mld, true kz, and L13 param.
                    ti = 1:tt;
                    
                    %'backdat' represents the conditions over the MLD
                    %residence time before the time of the current
                    %observation
                    time            = model.time(ti);
                    backdat.dt      = dt;
                    backdat.time    = model.time(ti(1)):backdat.dt:model.time(ti(end));

                    %THING WE KNOW (e.g., from reanalysis products)
                    backdat.mld_t   = interp1(time,model.mld_t(ti),backdat.time);
                    backdat.u10     = interp1(time,model.u10(ti),backdat.time);
                    backdat.slp     = interp1(time,model.slp(ti),backdat.time);
                    backdat.ice     = interp1(time,model.ice(ti),backdat.time);

                    %THINGS WE ESTIMATE / ASSUME / HOLD CONSTANT
                    backdat.mld     = repmat(model.mld(ti(end)),size(backdat.time));

                    backdat.mld_s   = repmat(model.mld_s(ti(end)),size(backdat.time));

                    %kz / upw [m/d]
                    mix.kz          = (model.kz + model.upw.*model.mld) .* 3600 .* 24./ (model.domain.dz); %/d
                    mix.kz          = repmat(nanmean(mix.kz(ti(end))),size(backdat.time));

                    %deep N2 and Ar [mmol/m3]
                    %Subsurface Ar is equilibrium concentration
                    mix.ardeep      = Arsol(model.deep_s(ti(end)),model.deep_t(ti(end))) .* model.deep_dens(ti(end)) ./ 1000;
                    mix.ardeep      = repmat(mix.ardeep,size(backdat.time));
                    %Subsurface N2 is based on dN2/Ar
                    mix.n2deep      = (model.dn2ar_deep(ti(end))/100+1) .* N2sol(model.deep_s(ti(end)),model.deep_t(ti(end))) .* model.deep_dens(ti(end)) ./ 1000;
                    mix.n2deep      = repmat(mix.n2deep,size(backdat.time));

                    %OTHER
                    %"true" N2 and Ar saturation, %
                    backdat.n2sat   = model.mld_n2sat(ti(end));
                    %Gas exchange parameterization
                    %choose from: l13, v10, s09, i11, n16 or w97
                    backdat.param = 'l13';
                    %Bubble flux scaling
                    backdat.beta = 0.5;

                    %RUN N2-prime
                    clear ardeep n2deep arsurf n2surf
                    [model.n2pr(tt),model.o2_tres(tt)] = n2_prime_2020(backdat,N(kk),mix);        
                    
                    %re-do N2-prime calculation by setting mixing component to 0
                    mix.kz(:) = 0;
                    [model.n2pr_noMix(tt),model.o2_tres_noMix(tt)] = n2_prime_2020(backdat,N(kk),mix);
                    
                end
            end  
            clear backdat
            
            %Plot results
                figure;
                ax(1)=subplot(4,1,1); hold on
                    plot(model.time,100*(model.mld_n2sat./model.mld_arsat-1),'k-','linewidth',1.5);
                    plot(model.time,100*(model.n2pr./model.mld_arsat-1),'r-','linewidth',1.5);
                    plot(model.time,100*(model.n2pr_noMix./model.mld_arsat-1),'b--','linewidth',1);
                    ylabel('\DeltaN_2/Ar [%]');
                    legend('Original','N2-prime','N2-prime (no mix.)','location','best')
                    title(strrep(sim{1},'_','-'))
                    set(gca,'box','on')
                subplot(4,1,2); hold on
                [ax([2,3]),p1,p2] = plotyy(model.time,[0,diff(model.mld)],model.time,[0,diff(model.mld_t)]);
                    set(p1,'color','k');set(p2,'color',[.5 .5 .5]);
                    set(ax(2),'ycolor','k');set(ax(3),'ycolor',[.5 .5 .5]);
                    ylabel(ax(2),'dMLD/dt [m/d]');
                    ylabel(ax(3),'dSST/dt [C/d]');
                    set(ax(2),'box','on')
                ax(4)=subplot(4,1,3); hold on
                    plot(model.time,model.o2_tres,'r-','linewidth',1.5)
                    plot(model.time,model.o2_tres_noMix,'b--','linewidth',1)
                    ylabel('O_2 re-equil. time [days]');
                    set(gca,'box','on')
                ax(5)=subplot(4,1,4); hold on
                    plot(model.time,-log10(model.kz),'k-','linewidth',1.5)
                    ylabel('-log_{10}(\kappa) [m2/s]');
                    xlabel('Model time [days]')
                    set(gca,'box','on')
                linkaxes(ax,'x')
                    
            %save output again
                save(fname,'model')
                
            %Save figure
                sn = [sim{kk},'_N2prime'];
                sname = save_fig(fig_save_dir,sn,{'j'}) 
                display(' ')
                display(['Model output figure saved as ' sname])

            clearvars -except kk sim N N1 mod_dir dt fig_save_dir
    end
    