function [climate, C_CO2] = climatedataBo( Out, indoorCO2, t )
%U interp the climate data
% 1-T_sky, 2- T_out, 3- Rad_out, 4- V_wind, 5- RH_out
% 
% Rad_out    =  interp1(Out(:,6),Out(:,3),t);            %outside global radiation
% T_out      =  interp1(Out(:,6),Out(:,2),t);            %outside air temperature
% V_wind     =  interp1(Out(:,6),Out(:,4),t);            %outside wind speed
% RH_out     =  interp1(Out(:,6),Out(:,5),t);            %outside relative humidity
% tic
climate = interp1(Out(:,6),Out(:,(1:5)),t);
% toc
% disp(['interp1:',num2str(toc)]);
C_CO2      =  interp1(indoorCO2(:,5),indoorCO2(:,4),t);%x(27+2*NrBox);                               %inside CO2 concentration [ppm]
% C_CO2out   = 450;                                     %outside CO2 concentration [ppm]



end

