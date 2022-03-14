function [ sens_element, sens_air ] = convectiveBo( Tair, Tout, area, tem, Vwind, l_spg, n_floor, n_canopy, n_walli, n_walle, n_roofi, n_roofe, ~, n_cover, n_blankete,u_blanket )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

global cin cout1 cout2 cout3

sens_element = zeros(l_spg,1);

sens_element(n_floor) = cin*abs(tem(n_floor)-Tair)^0.33*(tem(n_floor)-Tair)*area(n_floor);
sens_element(n_canopy)= 10*(tem(n_canopy)-Tair)*area(n_canopy);
sens_element(n_walli) = cin*abs(tem(n_walli)-Tair)^0.33*(tem(n_walli)-Tair)*area(n_walli);
sens_element(n_walle) = (cout1 +cout2*Vwind^cout3)*(tem(n_walle)-Tout)*area(n_walle);
sens_element(n_roofi) = cin*abs(tem(n_roofi)-Tair)^0.33*(tem(n_roofi)-Tair)*area(n_roofi);
sens_element(n_roofe) = (cout1 +cout2*Vwind^cout3)*(tem(n_roofe)-Tout)*area(n_roofe);
%sens_element(n_screen)= cin*abs(tem(n_screen)-Tair)^0.33*(tem(n_screen)-Tair)*area(n_screen)*(zerodivided0(abs(tem(n_screen)-Tair),tem(n_screen)-Tair));
sens_element(n_cover) = 7.2*(tem(n_cover)-Tair)*area(n_cover)+(cout1 +cout2*Vwind^cout3)*(tem(n_cover)-Tout)*area(n_cover)*(1-u_blanket);
sens_element(n_blankete)=(cout1 +cout2*Vwind^cout3)*(tem(n_blankete)-Tout)*area(n_blankete)*u_blanket;

% control?

sens_air = sens_element(n_floor)+sens_element(n_canopy)+sens_element(n_walli)+sens_element(n_roofi)+7.2*(tem(n_cover)-Tair)*area(n_cover);


sens_element = - sens_element;

end

