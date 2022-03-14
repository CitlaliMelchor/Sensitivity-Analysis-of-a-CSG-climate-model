function [u_blanket,u_vent, u_ventside, u_venttop, u_venttopside] = ControlventandblanketBo(U0,day_num,hour_clock)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Blanket_morning = interp1(U0(:,1),U0(:,8),day_num);           % time of opening in the morning
Blanket_afternoon = interp1(U0(:,1),U0(:,9),day_num);        % time of closing in the afternoon
Blanket = interp1(U0(:,1),U0(:,10),day_num);

%u_blanket=Blanket.*switch01((hour_clock-Blanket_morning)*(hour_clock-Blanket_afternoon),-1);    %switch01((day_num-daynumber(0,[3,1],1))*(day_num-daynumber(0,[11,6],1)),-1)
u_blanket  = round(Blanket/(1+exp(-20*(hour_clock-Blanket_morning)*(hour_clock-Blanket_afternoon))),1);

u_vent=interp1(U0(:,1),U0(:,5),day_num);

vent_start     =interp1(U0(:,1),U0(:,6),day_num);
vent_end       =interp1(U0(:,1),U0(:,7),day_num);
u_venttop      =interp1(U0(:,1),U0(:,2),day_num)/(1+exp(10*((hour_clock-vent_start)*(hour_clock-vent_end))));
u_ventside     =interp1(U0(:,1),U0(:,3),day_num)/(1+exp(10*((hour_clock-vent_start)*(hour_clock-vent_end))));
u_venttopside  =interp1(U0(:,1),U0(:,4),day_num)/(1+exp(1e16*((hour_clock-vent_start)*(hour_clock-vent_end))));

end

