function tran = directtran( angle_solar,angle_csg,cover,Material )
%directtran: cover transmission of direct solar light.
%Parameters: 
%angle_solar -[solar altitude,time angle, solar declination],
%angle_csg   -[greenhouse orientation(deg)], 
%cover       -[slop(rad); length(m)]
%Material    -[reflection index, thickness(m),power absorption coefficient(/m)],

solar_azimuth = asin(cos(angle_solar(:,3)).*sin(angle_solar(:,2)).*zerodivided0(1,cos(angle_solar(1))));                      % solar azimuth relative to N-S [rad]
angle_g       = angle_csg*pi/180;                                                                                             % greenhouse orientation to N-S [rad]
angle_inci    = acos(cos(angle_solar(1)).*sin(solar_azimuth+angle_g).*sin(cover(:,1))+sin(angle_solar(1)).*cos(cover(:,1)));  % angle of incidence between normal(perpendicular) line [rad]
c             = (Material(1)^2-sin(angle_inci).^2).^0.5;                                                                      % Material(1) is reflection index of plastic (optical property of the material)
Rpa_single    = ((cos(angle_inci)-c)./(cos(angle_inci)+c)).^2;                                                                % single reflection of parallel light
Rpe_single    = ((Material(1)^2*cos(angle_inci)-c)./(Material(1)^2*cos(angle_inci)+c)).^2;                                    % single reflection of perpendicular light
Pathway       = Material(2)./(1-(sin(angle_inci)/Material(1)).^2).^0.5;                                                       % the optical pathway [m]
absorp        = exp(-2*Material(3)*Pathway);                                                                                  % the absorption of the material
Rpa_mul       = Rpa_single+Rpa_single.*(1-Rpa_single).^2.*absorp./(1-Rpa_single.^2.*absorp);       % for multiple reflection with transmission through the pane-parallel
Rpe_mul       = Rpe_single+Rpe_single.*(1-Rpe_single).^2.*absorp./(1-Rpe_single.^2.*absorp);       % for multiple reflection with transmission through the pane-perpendicular
Tpa           = (1-Rpa_mul).^2.*absorp.^0.5./(1-Rpa_mul.*absorp);                                  % transmission of the parallel light
Tpe           = (1-Rpe_mul).^2.*absorp.^0.5./(1-Rpe_mul.*absorp);                                  % transmission of the perpendicular light

Tpa(Tpa<0)    = 0;               % transfor the negative value to 0
Tpe(Tpe<0)    = 0;

tran_tot      = (Tpa+Tpe)/2.*switch01(angle_solar(1),-1);                                              % total transmission of light (assuming non polarzed light)


%% total direct light transmission of cover

tran = sum(tran_tot.*cover(:,2).*sin(pi-cover(:,1)-angle_solar(1)))/(sum(cover(:,2).*sin(pi-cover(:,1)-angle_solar(1))));    


end

