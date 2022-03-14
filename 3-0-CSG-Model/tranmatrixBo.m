function [tran,spg,angle_solar] = tranmatrixBo(t,day_num,ini_time,ini_location,Angle_spg,shape_cover,spg,u_blanket,LAI,n_cover,n_canopy,n_floor)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
global n C_abs

angle_solar = solaraltBo(t,day_num,ini_time(1,1),ini_location);                                  % solar position

tran_cover = directtran(angle_solar,Angle_spg,shape_cover,[n,spg.thickness(n_cover,1),C_abs]);    % transmission of cover
tran_cover = 0.95*tran_cover;
% cover

spg.NIRtran(n_cover,1)  = tran_cover.*(1-u_blanket);
spg.PARtran(n_cover,1)  = tran_cover.*(1-u_blanket);
spg.NIRabs(n_cover,1)   = 0.1.*(1-u_blanket);    %1-tran_cover;
spg.PARabs(n_cover,1)   = 0.1.*(1-u_blanket);    %1-tran_cover;
spg.difNIRabs(n_cover,1) = 0.1.*(1-u_blanket); 
spg.difPARabs(n_cover,1) = 0.1.*(1-u_blanket); 
% canopy
spg.FIRtran(n_canopy,1)    = exp(-0.94*LAI);
spg.FIRabs(n_canopy,1)     = 1-exp(-0.94*LAI);
spg.PARtran(n_canopy,1)    = exp(-0.7*LAI);%exp(-(0.9+0.83*exp(-0.12*angle_solar(1)))*LAI);%exp(-0.7*LAI);
spg.PARabs(n_canopy,1)     = 1-exp(-0.7*LAI);%0.94-0.95*exp(-(0.88+2.6*exp(-0.18*angle_solar(1)))*LAI);
spg.NIRtran(n_canopy,1)    = exp(-0.27*LAI);%(0.05+0.06*exp(-0.08*angle_solar(1)))+(0.92-0.53*exp(-0.18*angle_solar(1)))*exp(-(0.48+0.54*exp(-0.13*angle_solar(1)))*LAI);%exp(-0.27*LAI);
spg.NIRabs(n_canopy,1)     = 1-exp(-0.27*LAI);%(0.67-0.06*exp(-0.08*angle_solar(1)))-(0.68-0.50*exp(-0.11*angle_solar(1)))*exp(-(0.25+0.38*exp(-0.12*angle_solar(1)))*LAI);
spg.difPARtran(n_canopy,1) = exp(-0.7*LAI);%exp(-0.92*LAI);%exp(-0.7*LAI);
spg.difPARabs(n_canopy,1)  = 1-exp(-0.7*LAI);%0.95-0.9*exp(-0.85*LAI);
spg.difNIRtran(n_canopy,1) = exp(-0.27*LAI);%0.05+0.91*exp(-0.5*LAI);%exp(-0.27*LAI);
spg.difNIRabs(n_canopy,1)  = 1-exp(-0.27*LAI);%0.65-0.65*exp(-0.27*LAI);

spg.area(n_canopy,1)       = LAI*spg.area(n_floor,1);

tran = [spg.FIRtran,spg.FIRabs,spg.NIRtran,spg.NIRabs,spg.PARtran,spg.PARabs,spg.difNIRtran,spg.difPARtran,spg.difNIRabs,spg.difPARabs];

end

