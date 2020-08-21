%This is an example NCP calculation based on simulated delO2/N2 and
%hydrographic data. This scrip uses the "calc_o2n2_ncp.m" function.

clear all; close all; clc

%--- Set default figure position
    set(groot, 'defaultFigurePosition',[100 0 700 700]) %figure position

%--- Set domain parameters
    Ndays   = 60; %days over which piston velocity is weighted
    dt      = .25; %time-increment, days
    mod     = 'w14'; %gas exchange model

%--- Simulated data
    tt      = 1:100; %time array
    do2n2   = nanmoving_average(.2*rand(1,100),2); %delO2/N2; %/100
    T       = nanmoving_average(rand(1,100)*20,2); %SST; deg-C
    S       = nanmoving_average(35 + (1-2*rand(1,100)),2); %Sal; PSU
    u10mat  = nanmoving_average(10 + (1-5*rand(100,Ndays * 1/dt + 1)),2); %u10 matrix; m/d
    mld     = nanmoving_average(50 + (1-10*rand(1,100)),2); %MLD

%--- Use function to calculate NCP
    [ncp,ko2,tro2,wt_t] = calc_o2n2_ncp(do2n2,S,T,mld,u10mat,Ndays,dt,mod);

%--- Plot
    subplot(8,1,1); hold on
        plot(tt,do2n2*100,'k','linewidth',2)
        ylabel('\DeltaO2/N2 [%]')
    subplot(8,1,2); hold on
        plot(tt,ncp,'r','linewidth',2)
        ylabel('NCP [mmol O2/m2/d]')
        set(gca,'yaxisloc','right')
    subplot(8,1,3); hold on
        plot(tt,T,'k','linewidth',2)
        ylabel('SST')
    subplot(8,1,4); hold on
        plot(tt,S,'k','linewidth',2)
        ylabel('Sal.')
        set(gca,'yaxisloc','right')
    subplot(8,1,5); hold on
        plot(tt,u10mat(:,end),'k','linewidth',2)
        plot(tt,u10mat);
        plot(tt,u10mat(:,end),'k','linewidth',2)
        ylabel('u10 [m/s]')
        legend('at time of underway obs.','during weighting period','location','eastoutside')
    subplot(8,1,6); hold on
        plot(tt,mld,'k','linewidth',2)
        ylabel('MLD [m]')
        set(gca,'yaxisloc','right')
        axis ij
    subplot(8,1,7); hold on
        plot(tt,tro2,'k','linewidth',2)
        plot(tt,wt_t,'b','linewidth',2)
        ylabel('Time [days]')
        legend('O2 residence time','Adjusted weighting time','location','eastoutside')
    subplot(8,1,8); hold on
        plot(tt,ko2,'k','linewidth',2)
        ylabel('ko2 [m/d]')
        set(gca,'yaxisloc','right')
        xlabel('Time')

    format_subplot(8,.02,.1,.1,.1,.25);