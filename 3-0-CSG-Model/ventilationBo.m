function [Vent,H_airvent] = ventilationBo( Atop_vent, Aside_vent, htop, hside, Afloor, V_wind, T_air, T_out, u_vent, u_venttopside, u_ventside, u_venttop, u_blanket)
%The natural ventilation of semi-passive greenhouse and air
%leaching/infiltration [m3/m2 s]
% Ventilation rate through both the roof and side vents (Kittas et al. 1997)

global Cd g Cw r_net Cleakage

Atop_vent  = Atop_vent*r_net;            %Maximum top ventilation vertical area [m2]
Aside_vent = Aside_vent*r_net;           %Maximum side ventilation vertical area[m2]

%% air leaching 
Cleakage1 = Cleakage*(1-u_blanket)+0.2*Cleakage*u_blanket;
Leaching_air = max(Cleakage1*0.25, Cleakage1*V_wind);      

%% top and side window
vent_topside = Cd/Afloor*((Atop_vent*Aside_vent/(Atop_vent^2+Aside_vent^2)^0.5)^2*...
    (2*g*(htop-hside)*(max((T_air-T_out),0))/(T_air/2+T_out/2+273.15))+(Atop_vent/2+Aside_vent/2)^2*Cw*V_wind^2)^0.5;

%% top window
% vent_top = Atop_vent*Cd/2/Afloor*(g*htop/2*abs(T_air-T_out)/((T_air+T_out)/2+273.15)+Cw*V_wind^2)^0.5;
vent_top = Atop_vent*Cd/2/Afloor*(g*htop/2*(max((T_air-T_out),0))/((T_air+T_out)/2+273.15)+Cw*V_wind^2)^0.5;
%% side window
vent_side = Cd/2/Afloor*Aside_vent*V_wind*Cw^0.5;
      
%% 

Vent = Leaching_air + (1-u_blanket)*u_vent*(vent_topside*u_venttopside + vent_side*u_ventside + vent_top*u_venttop);   % [m3/ m2 s]
% Vent(isnan(Vent)==1) = 0;  %if vent = nan then vent = 0;
global Cap_air dens_air
H_airvent = Afloor*Vent*(T_air-T_out)*Cap_air*dens_air*1000;

end

