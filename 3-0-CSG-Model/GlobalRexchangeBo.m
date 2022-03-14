function [Glob_R, Rad_canopy, PAR_canopy] = GlobalRexchangeBo( R, alt_solar,tran, angle1, length1,area,label1,label2,label3, n_cover, n_canopy, n_screen, n_blanketi, n_blankete, u_blanket, l_spg,outlayer ,LAI)
%The Global radiation absorption 
%   R -- the outside surface global solar radiation [W m-2]
%   alt_solar--solar altitude.
%   tran -- [FIRtran,FIRabs,NIRtran,NIRabs,PARtran,PARabs,difNIRtran,difPARtran,difNIRabs,difPARabs] transmission of direct light and diffuse light
%   angle1 -- the angle of each panel of greenhouse
%   length1 -- the length of each panel
%   label2 -- the internal or middel layer  [1-inner face, 0-middle, -1-outer face]
%   label3 -- the internal or middel layer  [1-inner face, 0-middle, -1-outer face, -1-cover, -1-thermal screen]
%   label4 -- the shadow layer is 1(wall_i) from the layer -1(north roof_e)
%% solar radiation
global r r_PAR 

s1    = 1-exp(-20*alt_solar);                  % smooth function    20 is the slope of the differentiable function
R_dir = r*s1*R;                                % direct solar radiation  [W m-2]
R_dif = R-R_dir;                               % diffuse solar radiation [W m-2]
Rad_para=zerodivided0(R_dir,sin(alt_solar));   % Solar radiation intensity parallel to the solar ray [W/m2]

%% direct light absorption
tran1       = [tran(n_canopy,:);tran(n_cover,:);tran(n_screen,:)];
trannir     = RtranBo( label1, tran1(:,3), label2 , n_blanketi, n_blankete, u_blanket,outlayer);
tranpar     = RtranBo( label1, tran1(:,5), label2 , n_blanketi, n_blankete, u_blanket,outlayer);
trandifnir  = RtranBo( label1, tran1(:,7), label2 , n_blanketi, n_blankete, u_blanket,outlayer);
trandifpar  = RtranBo( label1, tran1(:,8), label2 , n_blanketi, n_blankete, u_blanket,outlayer);

length1(n_canopy)    = area(n_canopy);

L1 = length1.*sin(angle1-alt_solar);
label2(n_blanketi) = 0;
label2(n_blankete) = 1;

L1 = L1.*label2;
L1(L1<=0)=0;
%shadow ?
[i,~]=find(label3==1);      % find the wall inside layer
L1(i)=sum(L1.*label3);      % the shadow from the roof external panel
if L1(L1<0)
L1(1)=L1(1)+L1(L1<0);                      % the shadow from the roof to floor
L1(n_canopy) = L1(n_canopy)+L1(L1<0)*LAI;  % the shadow to canopy

L1(L1<0)=0;                                % transfer the wall to 0, if the wall is shadowed by roof
end

Glob_R_dir = L1.*Rad_para.*trannir(:,l_spg).*(1-r_PAR).*tran(:,4)+L1.*Rad_para.*tranpar(:,l_spg).*r_PAR.*tran(:,6);
Glob_R_dir(n_canopy)=L1(n_canopy).*Rad_para.*trannir(n_canopy,l_spg).*(1-r_PAR).*tran(n_canopy,4);

%% diffuse light absorption

Glob_R_dif = length1.*R_dif.*(1-r_PAR).*trandifnir(:,l_spg).*tran(:,9)+length1.*R_dif.*r_PAR.*trandifpar(:,l_spg).*tran(:,10);
% Glob_R_dif(n_canopy) = length1(n_canopy).*R_dif.*(1-r_PAR).*trandifnir(n_canopy,l_spg).*tran(n_canopy,9);

%% 
Glob_R = Glob_R_dir+Glob_R_dif;                 %[W per meter length of greenhouse]
% Glob_R(n_canopy) = L1(n_canopy).*Rad_para.*trannir(n_canopy,l_spg).*(1-r_PAR).*tran(n_canopy,4) + length1(n_canopy).*R_dif.*(1-r_PAR).*trandifnir(n_canopy,l_spg).*tran(n_canopy,9);

Rad_canopy = Glob_R_dir(n_canopy)/L1(n_canopy).*sin(alt_solar)+Glob_R_dif(n_canopy)/length1(n_canopy);                        % the global radiation above the canopy  [W m-2]  ----evaporation model
% Rad_canopy = R_dir.*u_blanket.*tran(n_cover,3).*tran(n_screen,3)+R_dif.*u_blanket.*tran(n_cover,7).*tran(n_screen,7);    % u_screen?
PAR_canopy = Rad_para.*tranpar(n_canopy,l_spg).*r_PAR.*sin(alt_solar)*tran(n_canopy,6)+R_dif.*r_PAR.*trandifpar(n_canopy,l_spg)*tran(n_canopy,10);  % the PAR radiation above the canopy  [W m-2] ---plant model


end

