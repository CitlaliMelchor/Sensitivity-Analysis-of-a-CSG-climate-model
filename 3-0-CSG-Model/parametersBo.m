%% simulate 1 m length of a parameterically defined semi-passive greenhouse
% Goal: model to simulate statevariables over time ...
% Input: geometry of CSG and outdoor climate (air temperature, global radiation, wind speed, RH...)
% Output: Air temperature and RH inside CSG
% Author: Bo Zhou, revision2 <2018.12.01>
%% format
clc
% clear all
format long;

%% global parameters
global r r_PAR Sigma dens_air Cap_air T_soilbound
r          = 0.8;                      % the ratio of direct light of global radiation
r_PAR      = 0.5;                      % the ratio of PAR among global radiation
Sigma      = 5.67e-08;                 % Stefan Boltzmann constant
dens_air   = 1.3;                      % Air density                                                          [kg/m3]
Cap_air    = 1.0064;                   % Air specific heat                                                    [kJ/kgK]
T_soilbound= 26;                       % boundry soil temperature                                             [oC]

% natural ventilation parameters
global Cd g Cw r_net Cleakage
Cd         = 0.75;                     %discharge coefficient which depends on the greenhouse shape.
g          = 9.8;                      %acceleration of gravity [m/s2]
Cw         = 0.09;                     %Global wind pressure coefficient which depends on the greenhouse shape
r_net      = 0.6;                      %the vent area ratio of pest net
%air_leach  = 1.25;                     %the air leaching/infiltration ration                [h-1]
Cleakage   = 5e-4;                     %the leakage coefficient which depends on the greenhouse type

% vapour transpiration
global Mwater R gam Lat rb rsmin Cevap1 Cevap2 Cevap3_day Cevap3_night Cevap4_day Cevap4_night Cevap5_day Cevap5_night T_min_day T_min_night r1 Water_con Cst
Mwater     = 18;                       % Molar mass of water    [kg/kmol]
R          = 8.314e3;                  % molar gas constant     [J/kmolK]
gam        = 65.8;                     % psychometric constant   [Pa/K]
Lat        = 2.45e6;                   % latent heat of evaporation   [J/kg water]
rb         = 275;                      % boundary layer resistance of the canopy for vapour transport   [s/m]
rsmin      = 82;                       % minimum canopy resistance [s/m]        **vptranBo.m line 33
Cevap1     = 4.3;                      % coefficient of the stomatal resistance model to account for radiation effect  [W/m2]
Cevap2     = 0.54;                     % Coefficient of the stomatal resistance model to account for radiation effect   [W/m2]
Cevap3_day = 6.1e-7;                   % Coefficient of the stomatal resistance model to account CO2 effect   [ppm-2]
Cevap3_night = 1.1e-11; 
Cevap4_day   = 4.3e-6;                 % Coefficient of the stomatal resistance model to account for vapor pressure difference [Pa-2]
Cevap4_night = 5.2e-6;
Cevap5_day   = 2.3e-2;
Cevap5_night = 0.5e-2;
T_min_day    = 24.5;
T_min_night  = 33.6;
r1           = 275;                    % aerodynamic resistance above the soil [s/m]
Water_con    = 0.1;                    % Water content of the soil
Cst          = 1;%0.4;                    % stress factor ???

% direct light transmission
global n C_abs
n          = 1.48;%1.5; %1.48                                                              % reflection index of plastic (optical property of the material)
C_abs      = 2112;%2000;%2112;    2578;%200;                                                                % power absorption coefficient [m-1]

% heat exchange coefficient of convective 
global cin cout1 cout2 cout3
cin = 1.86;                               % Convective heat exchange parameter between greenhouse elements and inside air    [W m-2 K-1]
cout1 = 2.8;                              % Convective heat exchange parameter between greenhouse elements and outside air   [W m-2 K-1]
cout2 = 1.2;                                                                                                               % [J m-3 K-1]
cout3 = 1;                                                                                                                 % []

% PhotoSynthesis Buffer ---------------------------------------------------
global LAI_Max THETA eta_ppm_mgm3
LAI_Max      =    2.5;      % [m2 {leaf} m^{-2}] Maximal value of LAI
THETA        =    0.7; 
eta_ppm_mgm3 =    1.804; % Unit Transformation for CO2 Concentration

% CarboHydrate Buffer   ---------------------------------------------------
global cfruit_g cleaf_g cstemroot_g rg_stemroot rg_leaf rg_fruit Tsum_needed T24_S1 T24_b1 T24_S2 T24_b2 Tinst_S1 Tinst_b1 Tinst_S2 Tinst_b2 %tau_Tcan %LAI_Max Simulation
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
%tau_Tcan    = 86400;    % Time Constant of Averaging Filter (24 hour)
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
NrBox       = 10;        % [-] Number of development stages
cfruit_m    = 1.16e-7;
cleaf_m     = 3.47e-7;
cstemroot_m = 1.47e-7;
Resp_fac    = 1; %(1 -exp(-f_RGR*RGR));
Q10m        = 2;        % [-] Temperature Effect on Maintenance Respiration
% Leaf Harvest             ------------------------------------------------
global SLA 
SLA       = 2.66e-5;      % [m2 {leaf} mg^{-1} {CH2O}] Specific Leaf Area Index

global track
track = zeros();

global C_CO2out
C_CO2out = 400;
%% the material of Greenhouse
spg = readtable('CSGmaterial.xlsx','ReadRowNames',true);
Pspg = [8.5, 0.03, 0.15, 0.14, 0.12, 0.16, 0.2, 0.2; 4.6, 0.26, 0.28, 0.15, 0.1, 0.12, 0.09, 0.65];

%% climate data
% outdoor data = [ 1t/t_sky, 2tem, 3rad, 4wind, 5rh, 6time]
Out_center = xlsread('20181224totaldata',1);
Out_station = xlsread('20181224totaldata',3);
dewpoint_tem_center  = dptem(Out_center(:,2),Out_center(:,5));               %outside dew point temperature 
dewpoint_tem_station = dptem(Out_station(:,2),Out_station(:,5));               %outside dew point temperature 
T_sky_center       = (0.74+0.006.*dewpoint_tem_center).^0.25.*(Out_center(:,2)+273)-273;                                          %sky temperature
T_sky_station      = (0.74+0.006.*dewpoint_tem_station).^0.25.*(Out_station(:,2)+273)-273;                                          %sky temperature
Out_center(:,1)      = T_sky_center;
Out_station(:,1)     = T_sky_station;

% real greenhouse data = [1tem, 2rh, 3soilTem, 4CO2, 5time]
xreal_6_fangshan =  xlsread('20181224totaldata',4);
xreal_7_fangshan =  xlsread('20181224totaldata',5);
% control vector of ventilation and blanket = [1DOY, 2top vent, 3bottom vent,
% 4top&bottom vent, 5ventilate or not, 6vent start, 7vent end, 8blanket morning,
% 9blanket afternoon]
U_6 = xlsread('20181224totaldata',6);
U_7 = xlsread('20181224totaldata',7);

%% the geometry of semi-passive greenhouse
%1- the coordinate of the panel start and end points 
c = [spg.startpointX,spg.startpointY,spg.endpointX,spg.endpointY];
l_spg = length(spg.element); 
%2- view factor between each panel
vf_mat = viewfactorBo(c,spg.outlayer);
%3- angle of each panel
[angle1,length1] = angleBo(c);
%4- the greenhouse orientation from north to south      [deg]
Angle_spg = 90;
%5- the volum of inside air [m3 m-1]
Vair = vairBo(Pspg);
%6- the initial transimission
tran = [spg.FIRtran,spg.FIRabs,spg.NIRtran,spg.NIRabs,spg.PARtran,spg.PARabs,spg.difNIRtran,spg.difPARtran];
%7- the window position
Atop_vent  = 0.07;          % Maximum top ventilation vertical area [m2]
Aside_vent = 1.2;           % Maximum side ventilation vertical area[m2]
htop       = 4.38;          % vertical height of middle point of top vent [m]
hside      = 0.5;          % vertical height of middle point of side vent [m]
%8- the location of each element
n_floor   = find(strcmp(spg.element,'floor'));
n_soil_last = find(spg.labelcond==1,1,'last');           % find the last layer of soil  --- 1
n_canopy  = find(strcmp(spg.element,'canopy'));
n_walli   = find(strcmp(spg.element,'walli'));
n_walle   = find(strcmp(spg.element,'walle'));
n_roofi   = find(strcmp(spg.element,'roofi'));
n_roofe   = find(strcmp(spg.element,'roofe'));
n_screen  = find(strcmp(spg.element,'thermalS'));
n_cover   = find(strcmp(spg.element,'cover'));
n_blanketi= find(strcmp(spg.element,'blanketi'));
n_blankete= find(strcmp(spg.element,'blankete'));
%9- the cover shape of greenhouse
shape_cover = coverBo(Pspg);
%10- label1 for radiation transmission between multipul covers
spg.label1             = zeros(l_spg,1);
spg.label1(n_floor)    = 1;
spg.label1(n_cover)    = 2;
spg.label1(n_blanketi) = 5;
spg.label1(n_blankete) = 5;
spg.label1(end)        = 5;  % sky
%11- label2 for inter and external layers
spg.label2             = zeros(l_spg,1);
spg.label2(n_floor)    = 1;
spg.label2(n_canopy)   = 1;
spg.label2(n_walli)    = 1;
spg.label2(n_roofi)    = 1;
spg.label2(n_screen)   = 1;
spg.label2(n_canopy)   = 1;
spg.label2(n_cover)    = 1;
spg.label2(n_blanketi) = 1;
spg.label2(end)        = 1;  %sky
spg.label2(n_walle)    = -1;
spg.label2(n_roofe)    = -1;
%12- label3 for shadow
spg.label3             = zeros(l_spg,1);
spg.label3(n_walli)    = 1;
spg.label3(n_roofe)    = -1;

%% the capacity * density
Capd_air     = Cap_air*dens_air*Vair*1000;                                      % [J/K]
spg.capacity(2:n_soil_last) = spg.capacity(2:n_soil_last).*(1-Water_con)+Water_con*4.2;
%% initial time setting %%
Clock_start = 18;                                         %simulation start time of 24 hours                                     [hour]
Month       = 9;                                          %simulation start month among 12                                       []
Date        = 1;                                          %simulation start date in a month                                      []
day_num_s = daynumber(0,[Month,Date],Clock_start);


