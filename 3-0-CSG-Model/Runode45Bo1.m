%% ode45 function [I CHANGED ODE45 FOR ODE_RK4G]
[t,x]=ode_rk4G(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
    n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
    Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);

%% plot %%%

%% air temperature
figure;
plot(t,x(:,1),'-r')          %inside air temperature
hold on
plot(xreal(:,5),xreal(:,1),'-b');
box on,legend('T_{air}-simulation','T_{air}-measurement'),title('Temperature');

%% vapour pressure
vpsat_s      = vapourPsat(x(:,1));
vp_s         = vpsat_s.*x(:,2);
vpsat_m      = vapourPsat(xreal(:,1));
vp_m         = vpsat_m.*xreal(:,2);

figure
plot(t,vp_s,'-r')
hold on
plot(xreal(:,5),vp_m,'-b')
box on
legend('VP_{simulation}','VP_{measurement}')

%% soil temperature
figure
plot(t,x(:,(n_floor+5)),'-r')
hold on
plot(xreal(:,5),xreal(:,3),'ob')
box on,legend('T_{soil}-simulation','T_{soil}-measurement'),title('Soil Temperature');

%% distribution
% simulation = [x(:,1),x(:,6),vp_s,t];
% tem1    =  interp1(simulation(:,4),simulation(:,1),xreal(:,5)); 
% soitem1    =  interp1(simulation(:,4),simulation(:,2),xreal(:,5)); 
% vp1    =  interp1(simulation(:,4),simulation(:,3),xreal(:,5)); 
% s_data= [xreal(:,5),tem1,soitem1,vp1];
% m_data= [xreal(:,5),xreal(:,1),xreal(:,3),vp_m];
% a= s_data-m_data;
% a(:,4)=a(:,4)./1000;
% %%%%  normal kernel function
% [f,x2]=ksdensity(a(:,2));
% figure
% plot(x2,f)
% hold on
% [f,x3]=ksdensity(a(:,3));
% figure
% plot(x3,f)
% [f,x4]=ksdensity(a(:,4));
% figure
% plot(x4,f)
% ave = mean(a,1)   
% sd  = std (a,0,1) 
