function [dM_vp,dRH, M, H_lat] = vptranBo(Vent, T_air, Tcover, Tcanopy, Tsoil2, T_out, RH_air, RH_out, Acover, Rad_canopy, C_CO2, Capd_air, V_air, LAI, Afloor, Water_con,n_cover,n_floor,n_canopy,l_spg)
%% vapour evaporation and transpiration from soil and canopy, as well as
%condensation and ventilation
% M    = [condensation from cover, canopy transpiration, soil evaporation, vapour exchange from ventilation]  [kg/s per meter length of greenhouse];
% dM_vp = vapour exchange [kg/s per meter length of greenhouse]
%% 
HEC_covairin = 7.2 * Acover / Afloor;

global Mwater R gam Lat rb rsmin Cevap1 Cevap2 Cevap3_day Cevap3_night Cevap4_day Cevap4_night Cevap5_day Cevap5_night T_min_day T_min_night r1 Cst

tem        = [T_air, Tcover, Tcanopy, Tsoil2, T_out];
VPsat      = vapourPsat(tem);                                       % vapour pressure at saturation  [pa]
VPsat_out  = VPsat(5);                                              % outside vapour pressure at saturation [pa]
VPair      = RH_air * VPsat(1);                                     % inside air vapour pressure  [pa]
VPout      = RH_out * VPsat_out;                                    % outside air vapour pressure [pa]
Cvap_out   = VPout*Mwater/(R*(273.15+T_out));                       % outside vapour concentration [kg/m3]
Cvap_in    = VPair*Mwater/(R*(273.15+T_air));                       % inside vapour concentration [kg/m3]

%% condensation from cover
%MV_cov=(switch01((VPair-VPsat(2)),-1))*6.4e-9*HEC_covairin*(VPair-VPsat(2));   %the vapour exchange coefficient between air and cover   [kg/m2s]
MV_cov=max(0,6.4e-9*HEC_covairin*(VPair-VPsat(2)));
%% transpiration from canopy
Srs    = switch03((Rad_canopy-5),-1);                %switch function for night and day time
Cevap3 = Cevap3_night*(1-Srs) + Cevap3_day*Srs;
Cevap4 = Cevap4_night*(1-Srs) + Cevap4_day*Srs;
Cevap5 = Cevap5_night*(1-Srs) + Cevap5_day*Srs;
T_min  = T_min_night *(1-Srs) + T_min_day *Srs;

rf_R   = (Rad_canopy+Cevap1)/(Rad_canopy+Cevap2);  %resistance factor for high radiation levels
rf_CO2 = (1+Srs*Cevap3*(C_CO2-200)^2)*switch01(C_CO2-1100,1)+1.5*(1-switch01(C_CO2-1100,1));                %resistance factor for high CO2 levels
rf_VP  = 3.8*switch02(1+Cevap4*(VPsat(3)-VPair)^2-3.8,-1)+(1+Cevap4*(VPsat(3)-VPair)^2)*switch01(1+Cevap4*(VPsat(3)-VPair)^2-3.8,1);           %resistance factor for large vapour pressure
rf_T   = 1+Cevap5*(tem(3)-T_min)^2;                    %resistance factor for temperature
rs     = rsmin*rf_R*rf_CO2*rf_VP*rf_T;                 %stomatal resistance of the canopy for vapour transport                [s/m]
VEC_canair = 2*Capd_air/V_air*LAI/(Lat*gam*(rb+rs));  %vapour exchange coefficient between the canopy and air  [kg Pa/s]
dvp1 = max(VPsat(3)-VPair,0);
MV_canair  = Cst*VEC_canair*dvp1;   %canopy transpiration rate                 %[kg/m2 s]

%% evaporation from soil
r2 = 3.5*Water_con^(-2.3);%30+3.5*(Water_max/Water_con)^2.3;%500;  %the surface resistance of the bare soil [s/m]
beta_soil = 1/(1+r2/r1);
VEC_soilair =Capd_air/V_air/(Lat*gam*(r1+r2));                       %vapour exchange coefficient between the soil and air  [kg Pa/s]
dvp2 = max(VPsat(4)-VPair,0);
MV_soil =beta_soil*VEC_soilair*dvp2;     % % soil evaporation   %[kg/m2 s]

%% ventilation
MV_vent = Cvap_out-Cvap_in;   %vapour remove rate according to ventilation [kg/m3]

%% 

M     = [MV_cov*Afloor, MV_canair*Afloor, MV_soil*Afloor, MV_vent*Vent*Afloor];
dM_vp = M*[-1;1;1;1];                           %change of vapour mass    [kg/s per meter length of greenhouse]

H_lat = zeros(l_spg,1);
H_lat(n_cover)  = M(1)*Lat;
H_lat(n_canopy) = -M(2)*Lat;
H_lat(n_floor+5)  = -M(3)*Lat;

M_vp = Cvap_in * V_air;                                   % vapour mass inside           [kg]

VPair1  = (M_vp+dM_vp)/V_air*R*(T_air+273.15)/Mwater;    %new vapour pressure [pa]
% dRH  =  (VPair1/VPsat(1))+(1-VPair1/VPsat(1))/(1+exp(10*(VPsat(1)-VPair1)))-RH_air;   %switch function, if VPair1>VPsat1, RHair=100% 
dRH  =  min(VPair1/VPsat(1),0.99)-RH_air;
% dRH(isnan(dRH)==1) = 0;

end

