%% Execution Script for the Bram Vanthoor Model
%clean 
%% Simulation Time
global Simulation
Simulation.Start  = 0; % [s]
Simulation.End    = 100*24*3600; % [s]
Simulation.SimRes = 5*60; % [s]
Simulation.TimeV  = [Simulation.Start:Simulation.SimRes:Simulation.End]';
%% Generate Input Trajectories
% Temperature Trajectory
AvgTemp  = 18.2;    % [degree Celcius]
DayPer   = 24*3600; % [s]
DIF      = 14;      % Day Temp - Night Temp [degree Celcius]
Tcan     = AvgTemp + DIF/2*sin((1/DayPer)*Simulation.TimeV*2*pi);      % [degree Celcius]
% Radiation Trajectory
AvgRadPerDay   = 12.2;    % [MJ m^{-2} day^{-1}]
AvgRadPerSec   = 1000000*AvgRadPerDay/(24*3600);    % [J m^{-2} s^{-1}]
PARcan = AvgRadPerSec + AvgRadPerSec*sin(1/DayPer*Simulation.TimeV*2*pi);  % [MJ m^{-2} s^{-1}]
% CO2 Trajectory
AvgCO2   = 677; % [mumol mol^{-1}]
AmpCO2   = 301; % [mumol mol^{-1}]
CO2air   = AvgCO2 + AmpCO2*sin(1/DayPer*Simulation.TimeV*2*pi);      % [mumol mol^{-1}]
%% Define parameters 
% PhotoSynthesis Buffer ---------------------------------------------------
global LAI_Max THETA eta_ppm_mgm3
LAI_Max      =    3;      % [m2 {leaf} m^{-2}] Maximal value of LAI
THETA        =  0.7; 
eta_ppm_mgm3 = 1.804; % Unit Transformation for CO2 Concentration
% CarboHydrate Buffer   ---------------------------------------------------
global cfruit_g cleaf_g cstemroot_g rg_stemroot rg_leaf rg_fruit Tsum_needed T24_S1 T24_b1 T24_S2 T24_b2 Tinst_S1 Tinst_b1 Tinst_S2 Tinst_b2 tau_Tcan %LAI_Max Simulation
% Growth Inhibition (24 hour mean)
T24_S1      = 1.1587;   % [-]
T24_b1      = 15;       % [degree Celcius]
T24_S2      = -1.3904;   % [-]
T24_b2      = 24.5;     % [degree Celcius]
% Growth Inhibition (instanteneous)
Tinst_S1    = 0.869;   % [-]
Tinst_b1    = 10;       % [degree Celcius]
Tinst_S2    = -0.5793;   % [-]
Tinst_b2    = 34;       % [degree Celcius]
tau_Tcan    = 86400;    % Time Constant of Averaging Filter (24 hour)
Tsum_needed = 1035;      % degreedays needed for full partioning to generative parts
% Potential Growth
rg_fruit    = 0.328;              % mg CH20.s-1.m-2 hereby the maintenace is included, based upon N:\Projecten\adaptieve kas\IntkamOutput\ T1TomFastStandard.CSV!!!!
rg_leaf     = 0.0950;             % mg CH20.s-1.m-2
rg_stemroot = 0.0742;         % mg CH20.s-1.m-2
% Growth respiration coefficients based upon Heuvelink Phd, page 238
cfruit_g    = 0.27;%           (1.37-1)/1.37;
cleaf_g     = 0.28;%            (1.39-1)/1.39;
cstemroot_g = 0.30;%        (1.42-1)/1.42;
% Maintenance Respiration  ------------------------------------------------
global NrBox cfruit_m cleaf_m cstemroot_m Q10m Resp_fac
NrBox       = 2;        % [-] Number of development stages
cfruit_m    = 1.16e-7;
cleaf_m     = 3.47e-7;
cstemroot_m = 1.47e-7;
Resp_fac    = 1; %(1 -exp(-f_RGR*RGR));
Q10m      = 2;        % [-] Temperature Effect on Maintenance Respiration
% Leaf Harvest             ------------------------------------------------
global SLA % LAI_Max
SLA       = 2.66e-5;      % [m2 {leaf} mg^{-1} {CH2O}] Specific Leaf Area Index
% Plant Buffer
% global Simulation NrBox
% Tomato Harvest
% global Simulation
%% Initial Values 
LAI_Init = 2;
% for the state vector
%   x = [Cbuf            Cstemroot             Cleaf 
%        Cfruit                    Nfruit                
%        Char            Tsum_1              Tcan24_1
%                        Tsum_2              Tcan24_2];
xinit = [10e3/(10^6)     (0.5*LAI_Init/SLA)/(10^6)      (LAI_Init/SLA)/(10^6)  ...
         (NrBox*ones(1,2)*0.1)/(10^6)       NrBox*ones(1,2)*2   ...
         0               1035                  AvgTemp     ...     
                         1035                  AvgTemp];
%% Create Input Vectors 
% These vectors are aligned with the time vector in the fourth column of u;
% this will enable the interpolation of data in the StateSpaceSystem-file.
Input(:,1) = Tcan;                          Input(:,2) = PARcan;
Input(:,3) = CO2air;                        
%% Run the System Model 
StatesStored = zeros(length(Simulation.TimeV),length(xinit));
OutputFlowsStored  = zeros(length(Simulation.TimeV),12);
for TimeIter = 1:1:length(Simulation.TimeV)-1
    disp(['Evaluating: ',num2str(TimeIter)]);
    if TimeIter == 1
        StatesIn = xinit;
        StatesStored(1,:) = xinit;
        LAI = LAI_Init;
    end
    % Run the Model
    [StateOut,OutputFlows,LAIOut] = StateSpaceSystem(StatesIn,Input(TimeIter,:),TimeIter,LAI);
    % Store the Output of the Model
    StatesStored(TimeIter+1,:)      = StateOut;
    OutputFlowsStored(TimeIter+1,:) = OutputFlows;
    StatesIn = StateOut;
    LAI      = LAIOut;
end
% Parse State Vector
Cbuf         = StatesStored(:,1)*(10^6);
Cstemroot    = StatesStored(:,2)*(10^6);
Cleaf        = StatesStored(:,3)*(10^6);
Cfruit       = StatesStored(:,4:4+(NrBox-1))*(10^6);
Nfruit       = StatesStored(:,4+(NrBox-1)+1:4+2*(NrBox-1)+1);
Charv        = StatesStored(:,4+2*(NrBox-1)+2)*(10^6);
Tsum_1       = StatesStored(:,4+2*(NrBox-1)+3);
Tcan24_1     = StatesStored(:,4+2*(NrBox-1)+4);
Tsum_2       = StatesStored(:,4+2*(NrBox-1)+3);
Tcan24_2     = StatesStored(:,4+2*(NrBox-1)+4);
% Parse the Output Flows
MCleafhar        = OutputFlowsStored(:,1)*(10^6);
MCairbuf         = OutputFlowsStored(:,2)*(10^6); 
MCbuffruit_tot   = OutputFlowsStored(:,3)*(10^6); 
MCbufleaf        = OutputFlowsStored(:,4)*(10^6); 
MCbuf_stemroot   = OutputFlowsStored(:,5)*(10^6); 
MCstemroot_air_m = OutputFlowsStored(:,6)*(10^6); 
MCleafair_m      = OutputFlowsStored(:,7)*(10^6); 
MCfruit          = OutputFlowsStored(:,8:9)*(10^6); 
MCfruitair_m     = OutputFlowsStored(:,10:11)*(10^6); 
MCbufair_g       = OutputFlowsStored(:,12)*(10^6); 
%% Validation Set
VALIDATION = 1;
if VALIDATION
    Vali = load('Validation.mat');
end
%% Present Inputs + Temperature Depend States
Figure_Handle = figure('Name','Inputs & Temperature','NumberTitle','off','Color','w');
% Plot: Tcan (Input), Tcan24 (Output)
Ax1 = subplot(4,1,1);
hold all
% Tcan
Fig1Plot1 = plot(Simulation.TimeV./(3600*24),Tcan,...
      'DisplayName','Tcan','LineStyle','-','LineWidth',2,...
      'Marker','none','MarkerEdgeColor','auto','MarkerFaceColor','none');
% Tcan24
Fig1Plot2 = plot(Simulation.TimeV./(3600*24),Tcan24_1,...
      'DisplayName','Tcan24','LineStyle','-','LineWidth',2,...
      'Marker','none','MarkerEdgeColor','auto','MarkerFaceColor','none');
Fig1Plot6 = plot(Simulation.TimeV./(3600*24),Tcan24_2,...
      'DisplayName','Tcan24','LineStyle','-','LineWidth',2,...
      'Marker','none','MarkerEdgeColor','auto','MarkerFaceColor','none');
set(gca,'Color','w',...
    'FontUnits','points','FontSize',12,'FontWeight','normal',...
    'GridLineStyle','-','LineWidth',0.5,'XGrid','on','YGrid','on',...
    'XLim',[Simulation.Start./(3600*24) Simulation.End./(3600*24)]);legend('show');
xlabel('Time [Days]');ylabel('Temperature [degrees Celcius]');
% Plot: Temperature Sum (Output)
Ax2 = subplot(4,1,2);
Fig1Plot3 = plot(Simulation.TimeV./(3600*24),Tsum_1,...
      'DisplayName','Tsum','LineStyle','-','LineWidth',2,...
      'Marker','none','MarkerEdgeColor','auto','MarkerFaceColor','none');
Fig1Plot6 = plot(Simulation.TimeV./(3600*24),Tsum_2,...
      'DisplayName','Tsum','LineStyle','-','LineWidth',2,...
      'Marker','none','MarkerEdgeColor','auto','MarkerFaceColor','none');
set(gca,'Color','w',...
    'FontUnits','points','FontSize',12,'FontWeight','normal',...
    'GridLineStyle','-','LineWidth',0.5,'XGrid','on','YGrid','on',...
    'XLim',[Simulation.Start./(3600*24) Simulation.End./(3600*24)]);legend('show');
xlabel('Time [Days]');ylabel('Temperature Sum [degree Celcius days]');
% Plot: Radiation
Ax3 = subplot(4,1,3);
Fig1Plot4 = plot(Simulation.TimeV./(3600*24),PARcan,...
      'DisplayName','PARcan','LineStyle','-','LineWidth',2,...
      'Marker','none','MarkerEdgeColor','auto','MarkerFaceColor','none');
set(gca,'Color','w',...
    'FontUnits','points','FontSize',12,'FontWeight','normal',...
    'GridLineStyle','-','LineWidth',0.5,'XGrid','on','YGrid','on',...
    'XLim',[Simulation.Start./(3600*24) Simulation.End./(3600*24)]);legend('show');
xlabel('Time [Days]');ylabel('Radiation [umol photons m^{-2} s^{-1}]');
% Plot: CO2 Concentration
Ax4 = subplot(4,1,4);
Fig1Plot4 = plot(Simulation.TimeV./(3600*24),CO2air,...
      'DisplayName','CO2air','LineStyle','-','LineWidth',2,...
      'Marker','none','MarkerEdgeColor','auto','MarkerFaceColor','none');
set(gca,'Color','w',...
    'FontUnits','points','FontSize',12,'FontWeight','normal',...
    'GridLineStyle','-','LineWidth',0.5,'XGrid','on','YGrid','on',...
    'XLim',[Simulation.Start./(3600*24) Simulation.End./(3600*24)]);legend('show');
xlabel('Time [Days]');ylabel('Radiation [umol CO2 mol^{-1} air]');
% Link axes
linkaxes([Ax1 Ax2 Ax3 Ax4],'x');
%% Present Mass Flows
Figure_Handle = figure('Name','Mass Flows','NumberTitle','off','Color','w');
% Plot: Cbuf, MCAirBuf, MCBufAir_g, MCBuf_Stem, MCBuf_Leaf, MCBuf_Fruit
Ax1 = subplot(1,1,1);
hold all
% MCAirBuf
Fig2Plot1 = plot(Simulation.TimeV./(3600*24),MCairbuf,...
      'DisplayName','MCAir,Buf','LineStyle','-','LineWidth',2,...
      'Marker','none','MarkerEdgeColor','auto','MarkerFaceColor','none');
if VALIDATION
    Fig1PlotV = plot(Vali.TimeOut./(3600*24),Vali.MCAirBuf,...
          'DisplayName','VAL,MCAir,Buf','LineStyle','-','LineWidth',2,...
          'Marker','none','MarkerEdgeColor','auto','MarkerFaceColor','none');
end
% MCBuf_Stem
Fig2Plot3 = plot(Simulation.TimeV./(3600*24),MCbuf_stemroot,...
      'DisplayName','MCBuf,Stem','LineStyle','-','LineWidth',2,...
      'Marker','none','MarkerEdgeColor','auto','MarkerFaceColor','none');
if VALIDATION
    Fig2PlotV = plot(Vali.TimeOut./(3600*24),Vali.MCBuf_Stem,...
          'DisplayName','VAL,MCBuf,Stem','LineStyle','-','LineWidth',2,...
          'Marker','none','MarkerEdgeColor','auto','MarkerFaceColor','none');
end
% MCBuf_Leaf
Fig2Plot4 = plot(Simulation.TimeV./(3600*24),MCbufleaf,...
      'DisplayName','MCBuf,Leaf','LineStyle','-','LineWidth',2,...
      'Marker','none','MarkerEdgeColor','auto','MarkerFaceColor','none');
if VALIDATION
    Fig3PlotV = plot(Vali.TimeOut./(3600*24),Vali.MCBuf_Leaf,...
          'DisplayName','VAL,MCBuf,Leaf','LineStyle','-','LineWidth',2,...
          'Marker','none','MarkerEdgeColor','auto','MarkerFaceColor','none');
end
% MCBuf_Fruit
Fig2Plot5 = plot(Simulation.TimeV./(3600*24),MCfruit,...
      'DisplayName','MCBuf,Fruit','LineStyle','-','LineWidth',2,...
      'Marker','none','MarkerEdgeColor','auto','MarkerFaceColor','none'); 
if VALIDATION
    Fig4PlotV = plot(Vali.TimeOut./(3600*24),Vali.MCBuf_Fruit,...
          'DisplayName','VAL,MCBuf,Fruit','LineStyle','-','LineWidth',2,...
          'Marker','none','MarkerEdgeColor','auto','MarkerFaceColor','none');
end
set(gca,'Color','w',...
    'FontUnits','points','FontSize',12,'FontWeight','normal',...
    'GridLineStyle','-','LineWidth',0.5,'XGrid','on','YGrid','on',...
    'XLim',[Simulation.Start./(3600*24) Simulation.End./(3600*24)]);legend('show');
xlabel('Time [Days]');ylabel('CarboHydrate Mass Flow [mg CH2O m^{-2} s^{-1}]');
%% Present Buffer Contents
Figure_Handle = figure('Name','Buffer Contents','NumberTitle','off','Color','w');
% Plot: Cbuf, Cleaf, Cstemroot, Cfruit
Ax1 = subplot(2,2,1);
hold all
% Cbuf
Fig3Plot1 = plot(Simulation.TimeV./(3600*24),Cbuf,...
      'DisplayName','Cbuf','LineStyle','-','LineWidth',2,...
      'Marker','none','MarkerEdgeColor','auto','MarkerFaceColor','none');
if VALIDATION
    Fig3PlotV = plot(Vali.t./(3600*24),Vali.Cbuf,...
          'DisplayName','VAL,Cbuf','LineStyle','-','LineWidth',2,...
          'Marker','none','MarkerEdgeColor','auto','MarkerFaceColor','none');
end
set(gca,'Color','w',...
    'FontUnits','points','FontSize',12,'FontWeight','normal',...
    'GridLineStyle','-','LineWidth',0.5,'XGrid','on','YGrid','on',...
    'XLim',[Simulation.Start./(3600*24) Simulation.End./(3600*24)]);legend('show');
xlabel('Time [Days]');ylabel('CarboHydrate Amount [mg CH2O m^{-2}]');
Ax2 = subplot(2,2,2);
hold all
% Cleaf
Fig3Plot2 = plot(Simulation.TimeV./(3600*24),Cleaf,...
      'DisplayName','Cleaf','LineStyle','-','LineWidth',2,...
      'Marker','none','MarkerEdgeColor','auto','MarkerFaceColor','none');
if VALIDATION
    Fig3PlotV = plot(Vali.t./(3600*24),Vali.Cleaf,...
          'DisplayName','VAL,Cleaf','LineStyle','-','LineWidth',2,...
          'Marker','none','MarkerEdgeColor','auto','MarkerFaceColor','none');
end
set(gca,'Color','w',...
    'FontUnits','points','FontSize',12,'FontWeight','normal',...
    'GridLineStyle','-','LineWidth',0.5,'XGrid','on','YGrid','on',...
    'XLim',[Simulation.Start./(3600*24) Simulation.End./(3600*24)]);legend('show');
xlabel('Time [Days]');ylabel('CarboHydrate Amount [mg CH2O m^{-2}]');
Ax3 = subplot(2,2,3);
hold all
% Cstemroot
Fig3Plot3 = plot(Simulation.TimeV./(3600*24),Cstemroot,...
      'DisplayName','Cstemroot','LineStyle','-','LineWidth',2,...
      'Marker','none','MarkerEdgeColor','auto','MarkerFaceColor','none');
if VALIDATION
    Fig3PlotV = plot(Vali.t./(3600*24),Vali.Cstemroot,...
          'DisplayName','VAL,Cstemroot','LineStyle','-','LineWidth',2,...
          'Marker','none','MarkerEdgeColor','auto','MarkerFaceColor','none');
end
set(gca,'Color','w',...
    'FontUnits','points','FontSize',12,'FontWeight','normal',...s
    'GridLineStyle','-','LineWidth',0.5,'XGrid','on','YGrid','on',...
    'XLim',[Simulation.Start./(3600*24) Simulation.End./(3600*24)]);legend('show');
xlabel('Time [Days]');ylabel('CarboHydrate Amount [mg CH2O m^{-2}]');
Ax4 = subplot(2,2,4);
hold all
% Cfruit
Fig3Plot4 = plot(Simulation.TimeV./(3600*24),Cfruit,...
      'DisplayName','Cfruit','LineStyle','-','LineWidth',2,...
      'Marker','none','MarkerEdgeColor','auto','MarkerFaceColor','none');
if VALIDATION
    Fig3PlotV = plot(Vali.t./(3600*24),Vali.Cfruit,...
          'DisplayName','VAL,Cfruit','LineStyle','-','LineWidth',2,...
          'Marker','none','MarkerEdgeColor','auto','MarkerFaceColor','none');
end
set(gca,'Color','w',...
    'FontUnits','points','FontSize',12,'FontWeight','normal',...
    'GridLineStyle','-','LineWidth',0.5,'XGrid','on','YGrid','on',...
    'XLim',[Simulation.Start./(3600*24) Simulation.End./(3600*24)]);legend('show');
xlabel('Time [Days]');ylabel('CarboHydrate Amount [mg CH2O m^{-2}]');
% Link axes
linkaxes([Ax1 Ax2 Ax3 Ax4],'x');
%% Analyse Leafs
Figure_Handle = figure('Name','Leaf Analysis','NumberTitle','off','Color','w');
% Plot: MCleafhar, MCBuf_Leaf,MCleafair_m
Ax1 = subplot(1,1,1);
hold all
% MCleafhar
Fig4Plot1 = plot(Simulation.TimeV./(3600*24),MCleafhar,...
      'DisplayName','MCleafhar','LineStyle','-','LineWidth',2,...
      'Marker','none','MarkerEdgeColor','auto','MarkerFaceColor','none');
if VALIDATION
    Fig3PlotV = plot(Vali.TimeOut./(3600*24),Vali.MCleafhar,...
          'DisplayName','VAL,MCleafhar','LineStyle','-','LineWidth',2,...
          'Marker','none','MarkerEdgeColor','auto','MarkerFaceColor','none');
end
Fig4Plot2 = plot(Simulation.TimeV./(3600*24),MCbufleaf,...
      'DisplayName','MCBuf_Leaf','LineStyle','-','LineWidth',2,...
      'Marker','none','MarkerEdgeColor','auto','MarkerFaceColor','none');
if VALIDATION
    Fig3PlotV = plot(Vali.TimeOut./(3600*24),Vali.MCBuf_Leaf,...
          'DisplayName','VAL,MCBuf,Leaf','LineStyle','-','LineWidth',2,...
          'Marker','none','MarkerEdgeColor','auto','MarkerFaceColor','none');
end
Fig4Plot3 = plot(Simulation.TimeV./(3600*24),MCleafair_m,...
      'DisplayName','MCleafair_m','LineStyle','-','LineWidth',2,...
      'Marker','none','MarkerEdgeColor','auto','MarkerFaceColor','none');
if VALIDATION
    Fig3PlotV = plot(Vali.TimeOut./(3600*24),Vali.MCleafair_m,...
          'DisplayName','VAL,MCleafair,m','LineStyle','-','LineWidth',2,...
          'Marker','none','MarkerEdgeColor','auto','MarkerFaceColor','none');
end
set(gca,'Color','w',...
    'FontUnits','points','FontSize',12,'FontWeight','normal',...
    'GridLineStyle','-','LineWidth',0.5,'XGrid','on','YGrid','on',...
    'XLim',[Simulation.Start./(3600*24) Simulation.End./(3600*24)]);legend('show');
xlabel('Time [Days]');ylabel('CarboHydrate Mass Flow [mg CH2O m^{-2} s^{-1}]');