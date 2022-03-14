%% Linear One-at-a-time sensitivity analysis
% MSc Thesis Sensitivity analysis of a Chinese Solar Greenhouse
% Author: Citlali Melchor Ramírez 
% December 2021

clear all
clc
%% parameter and data
tic
 parametersBo
toc
 
% parameters_ShunyiBo

%Cevap4_day_SA = [1].*Cevap4_day;     %%%%%%%%**********************************************************************************
%for i=1:length(Cevap4_day_SA)
%Cevap4_day = Cevap4_day_SA(i);
%% initial setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------Plant----------%%%%
LAI0        = 2;         % initial LAI  0.2
Tsum        = 1000;     % initial temperature sum -1000;

%------Time-----------%%%%    The time is from 18:00 Sep. 1st to 18:00 10th
%Jan.
set_month   = 11;  % Month   From 9 - 12 
set_day     = 7;  % Day
set_hour    = 0;  % hour
set_daynum  = 1;   % day number from 1 to 132.

%------Data-----------%%%%
Out       = Out_center;        % Out_center         % Out_station
xreal     = xreal_7_fangshan;  % xreal_6_fangshan   % xreal_7_fangshan   %xreal_shunyi
U = U_7;         % U_6 is the No.6 greenhouse; 
                 % U_7 is the No.7 greenhouse.
                 % U_shunyi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initial input value
Tair0   = 28;
RH0     = 0.5;                           % initial inside air RH
Tcan24       = 25;
Cbuf         = 1e4;
Cstemroot    = LAI0/SLA;     %  1.5038e+04;
Cleaf        = LAI0/SLA;     %  1.5038e+04;
Char         = 0;
Cfruit       = ones(1,NrBox)*0.1;
Nfruit       = ones(1,NrBox);
x_plant0 = [Tsum, Tcan24, Cbuf, Cstemroot, Cleaf, Char, Cfruit, Nfruit];
l_plant = length (x_plant0);

%% location
latitude    = 39.62;                                        %local latitude -Beijing is 40 oN                                      [deg]
longitude   = 115.97;                                       %local longitude -Beijing is 117 oN                                    [deg]
Hourcircle  = 120;                                          %hour circle applied in the experiment area -Beijing is 120 oN         [deg] 
ini_location = [latitude, longitude, Hourcircle];

%% %%%%%%%%%%%%%Run the model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------------------------------%%%%%%%%%%%%%
% simulation time
day_num_s1 = daynumber(0,[set_month,set_day],set_hour);
ini_time = [set_hour, set_month, set_day];
%input data
xreal(:,5)= round(xreal(:,5)./100).*100;
% time distance
dis_time = 890;%8640;%600;
n1=set_daynum*24*3600/dis_time;

%%%%  %%%%%%%%
if set_daynum < 100
    
   dtime = (day_num_s1 - day_num_s)*24*3600+(set_hour-Clock_start)*3600;
   a1 = find(Out(:,6)==dtime);
   a2 = find(xreal(:,5)==dtime);
   a3 = Out(3,6)-Out(2,6);
   a4 = xreal(3,5)-xreal(2,5);
   b1 = a1+ceil(n1*dis_time/a3);
   b2 = a2+ceil(n1*dis_time/a4);

   Out = Out((a1:b1),:);
   xreal  = xreal((a2:b2),:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Out(:,6) = Out(:,6) - dtime;
   xreal(:,5) = xreal(:,5) - dtime;
end

%%%%%     %%%%%%%%%%%

%% ode 45
x0 = [Tair0,RH0,LAI0,spg.iniTem(1:(l_spg-1))',x_plant0]';
tspan=[0:dis_time:(n1-1)*dis_time]' ;

%% Nominal values
[t_nom,x_nom]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
    n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
    Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
save xnom.mat x_nom
% Delta
load parametersList.mat
delta = [0.8 0.9 1.1 1.2];
index = 1;
%% Densities sensitivity
spg_nom = spg;
for i = 1:22
    for j = 1:4
        spg.density(i) = spg.density(i)*delta(j);
        
        [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
            n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
            Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
        
        s_density(:,:,j) = (x_delta-x_nom)./(x_nom.*(delta(j)-1));
        
        figure (1)
        fig1 = scatter(t_nom, s_density(:,1,j),'.'); hold on
        %fig1.AlphaData=-s_density(:,1,j);
        %fig1.MarkerFaceAlpha = 'flat';

        spg = spg_nom;
    end
    S.(parametersList.Properties.RowNames{index}) = s_density;
    index = index+1
end

%% Capacities sensitivity
for i = 1:22
    for j = 1:4
        
        spg.capacity(i) = spg.capacity(i)*delta(j);
     
        [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
            n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
            Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
        s_capacity(:,:,j) = (x_delta-x_nom)./(x_nom.*(delta(j)-1));
        
        figure (1)
        fig1 = scatter(t_nom, s_capacity(:,1,j),'.'); hold on
        %fig1.AlphaData = -s_capacity(:,1,j);
        %fig1.MarkerFaceAlpha = 'flat';

        spg = spg_nom;
    end
    S.(parametersList.Properties.RowNames{index}) = s_capacity;
    index = index+1
end

%% Thickness sensitivity
for i = 1:22
    for j = 1:4
        
        spg.thickness(i) = spg.thickness(i)*delta(j);
     
        [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
            n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
            Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
        s_thickness(:,:,j) = (x_delta-x_nom)./(x_nom.*(delta(j)-1));
        
        figure (1)
        fig1 = scatter(t_nom, s_thickness(:,1,j),'.'); hold on
        %fig1.AlphaData = -s_thickness(:,1,j);
        %fig1.MarkerFaceAlpha = 'flat';

        spg = spg_nom;
    end
    S.(parametersList.Properties.RowNames{index}) = s_thickness;
    index = index+1
end

%% Area sensitivity
for i = 1:22
    for j = 1:4
        
        spg.area(i) = spg.area(i)*delta(j);
     
        [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
            n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
            Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
        s_area(:,:,j) = (x_delta-x_nom)./(x_nom.*(delta(j)-1));
        
        figure (1)
        fig1 = scatter(t_nom, s_area(:,1,j),'.'); hold on
        %fig1.AlphaData = -s_area(:,1,j);
        %fig1.MarkerFaceAlpha = 'flat';

        spg = spg_nom;
    end
    S.(parametersList.Properties.RowNames{index}) = s_area;
    index = index+1
end
%% FIRtran sensitivity
for i = 1:22
    for j = 1:4
        
        spg.FIRtran(i) = spg.FIRtran(i)*delta(j);
     
        [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
            n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
            Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
        s_FIRtran(:,:,j) = (x_delta-x_nom)./(x_nom.*(delta(j)-1));
        
        figure (1)
        fig1 = scatter(t_nom, s_FIRtran(:,1,j),'.'); hold on
        %fig1.AlphaData = -s_FIRtran(:,1,j);
        %fig1.MarkerFaceAlpha = 'flat';

        spg = spg_nom;
    end
    S.(parametersList.Properties.RowNames{index}) = s_FIRtran;
    index = index+1
end

%% FIRabs sensitivity
for i = 1:22
    for j = 1:4
        
        spg.FIRabs(i) = spg.FIRabs(i)*delta(j);
     
        [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
            n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
            Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
        s_FIRabs(:,:,j) = (x_delta-x_nom)./(x_nom.*(delta(j)-1));
        
        figure (1)
        fig1 = scatter(t_nom, s_FIRabs(:,1,j),'.'); hold on
        %fig1.AlphaData = -s_FIRabs(:,1,j);
        %fig1.MarkerFaceAlpha = 'flat';

        spg = spg_nom;
    end
    S.(parametersList.Properties.RowNames{index}) = s_FIRabs;
    index = index+1
end

%% NIRtran sensitivity
for i = 1:22
    for j = 1:4
        
        spg.NIRtran(i) = spg.NIRtran(i)*delta(j);
     
        [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
            n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
            Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
        s_NIRtran(:,:,j) = (x_delta-x_nom)./(x_nom.*(delta(j)-1));
        
        figure (1)
        fig1 = scatter(t_nom, s_NIRtran(:,1,j),'.'); hold on
        %fig1.AlphaData = -s_NIRtran(:,1,j);
        %fig1.MarkerFaceAlpha = 'flat';
        
        spg = spg_nom;
    end
    S.(parametersList.Properties.RowNames{index}) = s_NIRtran;
    index = index+1
end

%% NIRabs sensitivity
for i = 1:22
    for j = 1:4
        
        spg.NIRabs(i) = spg.NIRabs(i)*delta(j);
     
        [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
            n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
            Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
        s_NIRabs(:,:,j) = (x_delta-x_nom)./(x_nom.*(delta(j)-1));
        
        figure (1)
        fig1 = scatter(t_nom, s_NIRabs(:,1,j),'.'); hold on
        %fig1.AlphaData = -s_NIRabs(:,1,j);
        %fig1.MarkerFaceAlpha = 'flat';

        spg = spg_nom;
    end
    S.(parametersList.Properties.RowNames{index}) = s_NIRabs;
    index = index+1
end

%% PARtran sensitivity
for i = 1:22
    for j = 1:4
        
        spg.PARtran(i) = spg.PARtran(i)*delta(j);
     
        [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
            n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
            Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
        s_PARtran(:,:,j) = (x_delta-x_nom)./(x_nom.*(delta(j)-1));
        
        figure (1)
        fig1 = scatter(t_nom, s_PARtran(:,1,j),'.'); hold on
        %fig1.AlphaData = -s_PARtran(:,1,j);
        %fig1.MarkerFaceAlpha = 'flat';

        spg = spg_nom;
    end
    S.(parametersList.Properties.RowNames{index}) = s_PARtran;
    index = index+1
end

%% PARabs sensitivity
for i = 1:22
    for j = 1:4
        
        spg.PARabs(i) = spg.PARabs(i)*delta(j);
     
        [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
            n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
            Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
        s_PARabs(:,:,j) = (x_delta-x_nom)./(x_nom.*(delta(j)-1));
        
        figure (1)
        fig1 = scatter(t_nom, s_PARabs(:,1,j),'.'); hold on
        %fig1.AlphaData = -s_PARabs(:,1,j);
        %fig1.MarkerFaceAlpha = 'flat';
        
        spg = spg_nom;
    end
    S.(parametersList.Properties.RowNames{index}) = s_PARabs;
    index = index+1
end

%% difNIRtran sensitivity
for i = 1:22
    for j = 1:4
        
        spg.difNIRtran(i) = spg.difNIRtran(i)*delta(j);
     
        [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
            n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
            Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
        s_difNIRtran(:,:,j) = (x_delta-x_nom)./(x_nom.*(delta(j)-1));
        
        figure (1)
        fig1 = scatter(t_nom, s_difNIRtran(:,1,j),'.'); hold on
        %fig1.AlphaData = -s_difNIRtran(:,1,j);
        %fig1.MarkerFaceAlpha = 'flat';

        spg = spg_nom;
    end
    S.(parametersList.Properties.RowNames{index}) = s_difNIRtran;
    index = index+1
end

%% difNIRabs sensitivity
for i = 1:22
    for j = 1:4
        
        spg.difNIRabs(i) = spg.difNIRabs(i)*delta(j);
     
        [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
            n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
            Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
        s_difNIRabs(:,:,j) = (x_delta-x_nom)./(x_nom.*(delta(j)-1));
        
        figure (1)
        fig1 = scatter(t_nom, s_difNIRabs(:,1,j),'.'); hold on
        %fig1.AlphaData = -s_difNIRabs(:,1,j);
        %fig1.MarkerFaceAlpha = 'flat';

        spg = spg_nom;
    end
    S.(parametersList.Properties.RowNames{index}) = s_difNIRabs;
    index = index+1
end

%% difPARtran sensitivity
for i = 1:22
    for j = 1:4
        
        spg.difPARtran(i) = spg.difPARtran(i)*delta(j);
     
        [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
            n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
            Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
        s_difPARtran(:,:,j) = (x_delta-x_nom)./(x_nom.*(delta(j)-1));
        
        figure (1)
        fig1 = scatter(t_nom, s_difPARtran(:,1,j),'.'); hold on
        %fig1.AlphaData = -s_difPARtran(:,1,j);
        %fig1.MarkerFaceAlpha = 'flat';

        spg = spg_nom;
    end
    S.(parametersList.Properties.RowNames{index}) = s_difPARtran;
    index = index+1
end

%% difPARabs sensitivity
for i = 1:22
    for j = 1:4
        
        spg.difPARabs(i) = spg.difPARabs(i)*delta(j);
     
        [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
            n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
            Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
        s_difPARabs(:,:,j) = (x_delta-x_nom)./(x_nom.*(delta(j)-1));
        
        figure (1)
        fig1 = scatter(t_nom, s_difPARabs(:,1,j),'.'); hold on
        %fig1.AlphaData = -s_difPARabs(:,1,j);
        %fig1.MarkerFaceAlpha = 'flat';
        
        spg = spg_nom;
    end
    S.(parametersList.Properties.RowNames{index}) = s_difPARabs;
    index = index+1
end

%% heatcond sensitivity
for i = 1:22
    for j = 1:4
        
        spg.heatcond(i) = spg.heatcond(i)*delta(j);
     
        [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
            n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
            Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
        s_heatcond(:,:,j) = (x_delta-x_nom)./(x_nom.*(delta(j)-1));
        
        figure (1)
        fig1 = scatter(t_nom, s_heatcond(:,1,j),'.'); hold on
        %fig1.AlphaData = -s_heatcond(:,1,j);
        %fig1.MarkerFaceAlpha = 'flat';

        spg = spg_nom;
    end
    S.(parametersList.Properties.RowNames{index}) = s_heatcond;
    index = index+1
end

%% Convective heat exchange parameters sensitivity

cin_nom = 1.86;                               % Convective heat exchange parameter between greenhouse elements and inside air    [W m-2 K-1]
cout1_nom = 2.8;                              % Convective heat exchange parameter between greenhouse elements and outside air   [W m-2 K-1]
cout2_nom = 1.2;                                                                                                               % [J m-3 K-1]
cout3_nom = 1;

for i = 1:4  
    cin= cin*delta(i);
    
    [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
        n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
        Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
    s_cin(:,:,i) = (x_delta-x_nom)./(x_nom.*(delta(i)-1));
        
    figure (1)
    fig1 = scatter(t_nom, s_cin(:,1,i),'.'); hold on
    %fig1.AlphaData = -s_cin(:,1,i);
    %fig1.MarkerFaceAlpha = 'flat';

    cin = cin_nom; 
end
S.(parametersList.Properties.RowNames{index}) = s_cin;
index = index+1

for i = 1:4   
    cout1= cout1*delta(i);
     %Evaluate the model 
    [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
        n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
        Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
    s_cout1(:,:,i) = (x_delta-x_nom)./(x_nom.*(delta(i)-1));
        
    figure (1)
    fig1 = scatter(t_nom, s_cout1(:,1,i),'.'); hold on
    %fig1.AlphaData = -s_cout1(:,1,i);
    %fig1.MarkerFaceAlpha = 'flat';

    cout1 = cout1_nom;
end
S.(parametersList.Properties.RowNames{index}) = s_cout1;
index = index+1

for i = 1:4   
    cout2= cout2*delta(i);
     %Evaluate the model 
    [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
        n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
        Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
    s_cout2(:,:,i) = (x_delta-x_nom)./(x_nom.*(delta(i)-1));
        
    figure (1)
    fig1 = scatter(t_nom, s_cout2(:,1,i),'.'); hold on
    %fig1.AlphaData = -s_cout2(:,1,i);
    %fig1.MarkerFaceAlpha = 'flat';

    cout2 = cout2_nom;
end
S.(parametersList.Properties.RowNames{index}) = s_cout2;
index = index+1

for i = 1:4   
    cout3= cout3*delta(i);
     %Evaluate the model 
    [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
        n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
        Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
    s_cout3(:,:,i) = (x_delta-x_nom)./(x_nom.*(delta(i)-1));
        
    figure (1)
    fig1 = scatter(t_nom, s_cout3(:,1,i),'.'); hold on
    %fig1.AlphaData = -s_cout3(:,1,i);
    %fig1.MarkerFaceAlpha = 'flat';

    cout3 = cout3_nom;
end
S.(parametersList.Properties.RowNames{index}) = s_cout3;
index = index+1

%% Discharge coefficient sensitivity (which depends on the greenhouse shape).
Cd_nom = 0.75;

for i = 1:4   
    Cd= Cd*delta(i);

    [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
        n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
        Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
    s_Cd(:,:,i) = (x_delta-x_nom)./(x_nom.*(delta(i)-1));
        
    figure (1)
    fig1 = scatter(t_nom, s_Cd(:,1,i),'.'); hold on
    %fig1.AlphaData = -s_Cd(:,1,i);
    %fig1.MarkerFaceAlpha = 'flat';
    
    Cd = Cd_nom; 
end
S.(parametersList.Properties.RowNames{index}) = s_Cd;
index = index+1

%% Global wind pressure coefficient sensitivity(which depends on the greenhouse shape)
Cw_nom = 0.09;

for i = 1:4   
    Cw= Cw*delta(i);

    [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
        n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
        Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
    s_Cw(:,:,i) = (x_delta-x_nom)./(x_nom.*(delta(i)-1));
        
    figure (1)
    fig1 = scatter(t_nom, s_Cw(:,1,i),'.'); hold on
    %fig1.AlphaData = -s_Cw(:,1,i);
    %fig1.MarkerFaceAlpha = 'flat';
    
    Cw = Cw_nom; 
end
S.(parametersList.Properties.RowNames{index}) = s_Cw;
index = index+1

%% Vent area ratio of pest net sensitivity
r_net_nom = 0.6; 

for i = 1:4   
    r_net= r_net_nom*delta(i);

    [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
        n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
        Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
    s_rnet(:,:,i) = (x_delta-x_nom)./(x_nom.*(delta(i)-1));
        
    figure (1)
    fig1 = scatter(t_nom, s_rnet(:,1,i),'.'); hold on
    %fig1.AlphaData = -s_rnet(:,1,i);
    %fig1.MarkerFaceAlpha = 'flat';
    
    r_net = r_net_nom; 
end
S.(parametersList.Properties.RowNames{index}) = s_rnet;
index = index+1

%% Leakage coefficient which depends on the greenhouse type
Cleakage_nom = 5e-4;

for i = 1:4   
    Cleakage= Cleakage_nom*delta(i);

    [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
        n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
        Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
    s_Cleakage(:,:,i) = (x_delta-x_nom)./(x_nom.*(delta(i)-1));
        
    figure (1)
    fig1 = scatter(t_nom, s_Cleakage(:,1,i),'.'); hold on
    %fig1.AlphaData = -s_Cleakage(:,1,i);
    %fig1.MarkerFaceAlpha = 'flat';
    
    Cleakage = Cleakage_nom; 
end

S.(parametersList.Properties.RowNames{index}) = s_Cleakage;
index = index+1

%% Direct light transmission sensitivity (?)
n_nom = 1.48; % reflection index of plastic (optical property of the material)
C_abs_nom = 2112;

for i = 1:4   
    n= n*delta(i);
     %Evaluate the model 
    [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
        n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
        Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
    s_n(:,:,i) = (x_delta-x_nom)./(x_nom.*(delta(i)-1));
        
    figure (1)
    fig1 = scatter(t_nom, s_n(:,1,i),'.'); hold on
    %fig1.AlphaData = -s_n(:,1,i);
    %fig1.MarkerFaceAlpha = 'flat';
    
    n = n_nom; 
end
S.(parametersList.Properties.RowNames{index}) = s_n;
index = index+1

for i = 1:4   
    C_abs= C_abs*delta(i);
     %Evaluate the model 
    [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
        n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
        Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
    s_C_abs(:,:,i) = (x_delta-x_nom)./(x_nom.*(delta(i)-1));
        
    figure (1)
    fig1 = scatter(t_nom, s_C_abs(:,1,i),'.'); hold on
    %fig1.AlphaData = -s_C_abs(:,1,i);
    %fig1.MarkerFaceAlpha = 'flat';
    
    C_abs = C_abs_nom;
end
S.(parametersList.Properties.RowNames{index}) = s_C_abs;
index = index+1

%% Sky temperature sensitivity
% T_sky_center_nom = T_sky_center;
% T_sky_station_nom = T_sky_station;
% 
% for i = 1:4   
%     T_sky_center= T_sky_center*delta(i);
%     Out_center(:,1)= T_sky_center;
%      %Evaluate the model 
%     [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
%         n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
%         Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
%     
%     s_TskyCenter(1,i,:) = (x_delta(3,:)-x_nom(3,:))/(T_sky_center(5)-T_sky_center_nom(1));
%     s_TskyCenter_norm(1,i,:) = ((x_delta(3,:)-x_nom(3,:))/(T_sky_center(5)-T_sky_center_nom(5))).*(T_sky_center_nom(5)./x_nom(3,:));
%     T_sky_center= T_sky_center_nom;
% 
% end
% 
% for i = 1:4   
%     T_sky_station= T_sky_station*delta(i);
%     Out_station(:,1)= T_sky_station;
%      %Evaluate the model 
%     [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
%         n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
%         Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
%     
%     s_TskyStation(1,i,:) = (x_delta(3,:)-x_nom(3,:))/(T_sky_station(5)-T_sky_station_nom(5));
%     s_TskyStation_norm(1,i,:) = ((x_delta(3,:)-x_nom(3,:))/(T_sky_station(5)-T_sky_station_nom(5))).*(T_sky_station_nom(5)./x_nom(3,:));
%     T_sky_station= T_sky_station_nom;
% end

%% Geometry of the GH sensitivities
% Window positions
Atop_vent_nom = 0.07;       % Maximum top ventilation vertical area [m2]
Aside_vent_nom = 1.2;       % Maximum side ventilation vertical area[m2]
htop_nom   = 4.38;          % vertical height of middle point of top vent [m]
hside_nom  = 0.5;

for i=1:4
    Atop_vent = Atop_vent*delta(i);

    [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
        n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
        Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
    s_Atop_vent(:,:,i) = (x_delta-x_nom)./(x_nom.*(delta(i)-1));
        
    figure (1)
    fig1 = scatter(t_nom, s_Atop_vent(:,1,i),'.'); hold on
    %fig1.AlphaData = -s_Atop_vent(:,1,i);
    %fig1.MarkerFaceAlpha = 'flat';
    
    Atop_vent = Atop_vent_nom;
end
S.(parametersList.Properties.RowNames{index}) = s_Atop_vent;
index = index+1
        
for i=1:4
    Aside_vent = Aside_vent*delta(i);

    [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
        n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
        Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
    s_Aside_vent(:,:,i) = (x_delta-x_nom)./(x_nom.*(delta(i)-1));
        
    figure (1)
    fig1 = scatter(t_nom, s_Aside_vent(:,1,i),'.'); hold on
    %fig1.AlphaData = -s_Aside_vent(:,1,i);
    %fig1.MarkerFaceAlpha = 'flat';
    
    Aside_vent = Aside_vent_nom;     
end
S.(parametersList.Properties.RowNames{index}) = s_Aside_vent;
index = index+1

for i=1:4
    htop = htop*delta(i);

    [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
        n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
        Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
    s_htop(:,:,i) = (x_delta-x_nom)./(x_nom.*(delta(i)-1));
        
    figure (1)
    fig1 = scatter(t_nom, s_htop(:,1,i),'.'); hold on
    %fig1.AlphaData = -s_htop(:,1,i);
    %fig1.MarkerFaceAlpha = 'flat';
    
    htop = htop_nom;      
end
S.(parametersList.Properties.RowNames{index}) = s_htop;
index = index+1

for i=1:4
    hside = hside*delta(i);
     %Evaluate the model 
    [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
        n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
        Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
    s_hside(:,:,i) = (x_delta-x_nom)./(x_nom.*(delta(i)-1));
        
    figure (1)
    fig1 = scatter(t_nom, s_hside(:,1,i),'.'); hold on
    %fig1.AlphaData = -s_hside(:,1,i);
    %fig1.MarkerFaceAlpha = 'flat';
    
    hside = hside_nom;     
end
S.(parametersList.Properties.RowNames{index}) = s_hside;
index = index+1

%% Aerodynamic resistance above the soil
r1_nom = 275;

for i=1:4
    r1 = r1*delta(i);

    [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
        n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
        Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
    s_r1(:,:,i) = (x_delta-x_nom)./(x_nom.*(delta(i)-1));
        
    figure (1)
    fig1 = scatter(t_nom, s_r1(:,1,i),'.'); hold on
    %fig1.AlphaData = -s_r1(:,1,i);
    %fig1.MarkerFaceAlpha = 'flat';
    
    r1 = r1_nom;  
end
S.(parametersList.Properties.RowNames{index}) = s_r1;
index = index+1

%% Greenhouse orientation sensitivity
Angle_spg_nom = 90;

for i=1:4
   Angle_spg = Angle_spg*delta(i);
     %Evaluate the model 
    [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
        n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
        Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
    s_anglespg(:,:,i) = (x_delta-x_nom)./(x_nom.*(delta(i)-1));
        
    figure (1)
    fig1 = scatter(t_nom, s_anglespg(:,1,i),'.'); hold on
    %fig1.AlphaData = -s_anglespg(:,1,i);
    %fig1.MarkerFaceAlpha = 'flat';
    
    Angle_spg = Angle_spg_nom;  
end
S.(parametersList.Properties.RowNames{index}) = s_anglespg;
index = index+1

%% Boundary soil temperature sensitivity 
T_soilbound_nom = 26;

for i=1:4
    
    T_soilbound = T_soilbound*delta(i);

    [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
        n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
        Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
    s_T_soilbound(:,:,i) = (x_delta-x_nom)./(x_nom.*(delta(i)-1));
        
    figure (1)
    fig1 = scatter(t_nom, s_T_soilbound(:,1,i),'.'); hold on
    %fig1.AlphaData = -s_T_soilbound(:,1,i);
    %fig1.MarkerFaceAlpha = 'flat';
    
    T_soilbound = T_soilbound_nom;  
end
S.(parametersList.Properties.RowNames{index}) = s_T_soilbound;
index = index+1

%% Boundary layer resistance of the canopy sensitivity 
rb_nom = 275;

for i=1:4
    rb = rb*delta(i);

    [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
        n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
        Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
    s_rb(:,:,i) = (x_delta-x_nom)./(x_nom.*(delta(i)-1));
        
    figure (1)
    fig1 = scatter(t_nom, s_rb(:,1,i),'.'); hold on
    %fig1.AlphaData = -s_rb(:,1,i);
    %fig1.MarkerFaceAlpha = 'flat';
    
    rb = rb_nom;  
end
S.(parametersList.Properties.RowNames{index}) = s_rb;
index = index+1

%% Minimum canopy resistance sensitivity 
rsmin_nom = 275;

for i=1:4
    rsmin = rsmin*delta(i);

    [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
        n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
        Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
    s_rsmin(:,:,i) = (x_delta-x_nom)./(x_nom.*(delta(i)-1));
        
    figure (1)
    fig1 = scatter(t_nom, s_rsmin(:,1,i),'.'); hold on
    %fig1.AlphaData = -s_rsmin(:,1,i);
    %fig1.MarkerFaceAlpha = 'flat';
    
    rsmin = rsmin_nom;  
end
S.(parametersList.Properties.RowNames{index}) = s_rsmin;
index = index+1

%% Coefficients of the stomatal resistance model 
Cevap1_nom = Cevap1;
Cevap2_nom = Cevap2;
Cevap3_day_nom = Cevap3_day;
Cevap3_night_nom = Cevap3_night;
Cevap4_day_nom = Cevap4_day;
Cevap4_night_nom = Cevap4_night;
Cevap5_day_nom = Cevap5_day;
Cevap5_night_nom = Cevap5_night;

for i=1:4
    Cevap1 = Cevap1*delta(i);
    
    [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
    n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
    Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
    s_Cevap1(:,:,i) = (x_delta-x_nom)./(x_nom.*(delta(i)-1));
        
    Cevap1 = Cevap1_nom; 
end
S.(parametersList.Properties.RowNames{index}) = s_Cevap1;
index = index+1;

for i=1:4
    Cevap2 = Cevap2*delta(i);
    
    [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
    n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
    Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
    s_Cevap2(:,:,i) = (x_delta-x_nom)./(x_nom.*(delta(i)-1));
        
    Cevap2 = Cevap2_nom; 
end
S.(parametersList.Properties.RowNames{index}) = s_Cevap2;
index = index+1;

for i=1:4
    Cevap3_day = Cevap3_day*delta(i);
    
    [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
    n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
    Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
    s_Cevap3_day(:,:,i) = (x_delta-x_nom)./(x_nom.*(delta(i)-1));
        
    Cevap3_day = Cevap3_day_nom; 
end
S.(parametersList.Properties.RowNames{index}) = s_Cevap3_day;
index = index+1;

for i=1:4
    Cevap3_night = Cevap3_night*delta(i);
    
    [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
    n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
    Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
    s_Cevap3_night(:,:,i) = (x_delta-x_nom)./(x_nom.*(delta(i)-1));
        
    Cevap3_night = Cevap3_night_nom; 
end
S.(parametersList.Properties.RowNames{index}) = s_Cevap3_night;
index = index+1;

for i=1:4
    Cevap4_day = Cevap4_day*delta(i);
    
    [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
    n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
    Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
    s_Cevap4_day(:,:,i) = (x_delta-x_nom)./(x_nom.*(delta(i)-1));
        
    Cevap4_day = Cevap4_day_nom; 
end
S.(parametersList.Properties.RowNames{index}) = s_Cevap4_day;
index = index+1;

for i=1:4
    Cevap4_night = Cevap4_night*delta(i);
    
    [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
    n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
    Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
    s_Cevap4_night(:,:,i) = (x_delta-x_nom)./(x_nom.*(delta(i)-1));
        
    Cevap4_night = Cevap4_night_nom; 
end
S.(parametersList.Properties.RowNames{index}) = s_Cevap4_night;
index = index+1;

for i=1:4
    Cevap5_day = Cevap5_day*delta(i);
    
    [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
    n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
    Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
    s_Cevap5_day(:,:,i) = (x_delta-x_nom)./(x_nom.*(delta(i)-1));
        
    Cevap5_day = Cevap5_day_nom; 
end
S.(parametersList.Properties.RowNames{index}) = s_Cevap5_day;
index = index+1;

for i=1:4
    Cevap5_night = Cevap5_night*delta(i);
    
    [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
    n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
    Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
    s_Cevap5_night(:,:,i) = (x_delta-x_nom)./(x_nom.*(delta(i)-1));
        
    Cevap5_night = Cevap5_night_nom; 
end
S.(parametersList.Properties.RowNames{index}) = s_Cevap5_night;
index = index+1;
%% Minimum temperatures
T_min_day_nom = T_min_day;
T_min_night_nom = T_min_night;

for i=1:4
    T_min_day = T_min_day*delta(i);
    
    [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
    n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
    Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
    s_T_min_day(:,:,i) = (x_delta-x_nom)./(x_nom.*(delta(i)-1));
        
    T_min_day = T_min_day_nom; 
end
S.(parametersList.Properties.RowNames{index}) = s_T_min_day;
index = index+1;

for i=1:4
    T_min_night = T_min_night*delta(i);
    
    [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
    n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
    Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
    s_T_min_night(:,:,i) = (x_delta-x_nom)./(x_nom.*(delta(i)-1));
        
    T_min_night = T_min_night_nom; 
end
S.(parametersList.Properties.RowNames{index}) = s_T_min_night;
index = index+1;
%% Water content and stress factor (?)
Water_con_nom = Water_con;
Cst_nom = Cst;

for i=1:4
    Water_con = Water_con*delta(i);
    
    [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
    n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
    Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
    s_Water_con(:,:,i) = (x_delta-x_nom)./(x_nom.*(delta(i)-1));
        
    Water_con = Water_con_nom; 
end
S.(parametersList.Properties.RowNames{index}) = s_Water_con;
index = index+1;

for i=1:4
    Cst = Cst*delta(i);
    
    [t_delta,x_delta]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
    n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
    Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    
    s_Cst(:,:,i) = (x_delta-x_nom)./(x_nom.*(delta(i)-1));
        
    Cst = Cst_nom; 
end
S.(parametersList.Properties.RowNames{index}) = s_Cst;
index = index+1;
%% Save the sensitivities data
 save S20.mat S
%% Load the sensitivities data (uncomment to run from this point)
load S20.mat



%% Calculate the average sensitivities (per day)
fn = fieldnames(S);
clear Smean S_temperature S_rh S_tsoil
for i=1:numel(fn)
    meanS = squeeze(mean(S.(fn{i}),1));
    Smean.(fn{i}) = meanS;

    S_temperature(i,:) =  meanS(1,:);
    S_rh(i,:) = meanS(2,:);
    S_tsoil(i,:) = meanS(5,:);
end

%% fix the files names
newFn = strrep(fn,'_','_{');
for i=1:length(fn)
    if contains(newFn(i),'_')
     newFn(i)=strcat(newFn(i),'}');
    end
end
%%
% %% Merge all sensitivities
% rowNames = {};
% counter=1;
% for i=[2:3 5:15]
%     for j=1:22
%         rowNames(counter,1)={sprintf('%s_{%s}', spg.Properties.RowNames{j},spg.Properties.VariableNames{i})};
%         counter=counter+1;
%     end
% end
% rowNames(end+1:end+18)={'cin';'cout1';'cout2';'cout3';'Cd';'Cw';'r_{net}';...
%     'Cleakage';'n'; 'C_{abs}';'T_{sky, center}'; 'T_{sky, station}';...
%     'Aside_{vent}';'Atop_{vent}';'htop';'hside';'r1';'Angle_{spg}'};
% 
% S = [s_cap_norm; s_thickness_norm; s_FIRtran_norm; s_FIRabs_norm;...
%      s_NIRtran_norm; s_NIRabs_norm; s_PARtran_norm; s_PARabs_norm;...
%      s_difNIRtran_norm; s_difNIRabs_norm; s_difPARtran_norm;...
%      s_difPARabs_norm; s_heatcond_norm; s_cHEC_norm; s_Cd_norm;...
%      s_Cw_norm; s_rnet_norm; s_Cleakage_norm; s_dir_norm;...
%      s_TskyCenter_norm; s_TskyStation_norm;...
%      s_Asidevent_norm; s_Atopvent_norm;s_htop_norm;s_hside_norm;...
%      s_r1_norm;s_anglespg_norm];
%  
S_temperature=array2table(S_temperature);
S_temperature.Properties.VariableNames={'0.1','1','1.1','2'};
S_temperature.Properties.RowNames=fn;%newFn;
S_temperature=sortrows(S_temperature,4,'descend','ComparisonMethod','abs');
S_temperature = S_temperature(~any(ismissing(S_temperature),2),:);
 save S_temperature.mat S_temperature
% 
S_rh=array2table(S_rh);
S_rh.Properties.VariableNames={'0.1','1','1.1','2'};
S_rh.Properties.RowNames=newFn;
S_rh=sortrows(S_rh,4,'descend','ComparisonMethod','abs');
S_rh = S_rh(~any(ismissing(S_rh),2),:);
 save S_rh.mat S_rh
% 
S_soilT=array2table(S_tsoil);
S_soilT.Properties.VariableNames={'0.1','1','1.1','2'};
S_soilT.Properties.RowNames=newFn;
S_soilT=sortrows(S_soilT,4,'descend','ComparisonMethod','abs');
S_soilT = S_soilT(~any(ismissing(S_soilT),2),:);
 save S_soilT.mat S_soilT
% 
%% Plot
negativeColor = ([[242 191 150];[234 150  81];[210 109 25];[140 72 16]])./255;
positiveColor = ([[209 218 180];[179 194 131];[143 163 80];[95 108 53]])./255;
% 
figure
% subplot (3,1,1)
b = bar(table2array(S_temperature(1:20,:)),'FaceColor','flat');
for i=1:20
    for ii=1:4
        if table2array(S_temperature(i,ii))<0
            newColor(ii,:) = negativeColor(ii,:);
        else
            newColor(ii,:) = positiveColor(ii,:);
        end
    end
    for k=1:4
        b(k).CData(i,:) = newColor(k,:);
    end
end
title('Sensitivity Analysis: Air Temperature')
xticks([1:1:290])
xticklabels(S_temperature.Properties.RowNames)
xtickangle(45)
ylim([-0.15 0.15])
ylabel('$\bar{S},[-]$','Interpreter','Latex')
grid on

figure
% subplot (3,1,2)
b = bar(table2array(S_rh(1:20,:)),'FaceColor','flat');
for i=1:20
    for ii=1:4
        if table2array(S_rh(i,ii))<0
            newColor(ii,:) = negativeColor(ii,:);
        else
            newColor(ii,:) = positiveColor(ii,:);
        end
    end
    for k=1:4
        b(k).CData(i,:) = newColor(k,:);
    end
end
title('Sensitivity Analysis: Relative Humidity')
xticks([1:1:290])
xticklabels(S_rh.Properties.RowNames)
xtickangle(45)
ylim([-0.15 0.15])
ylabel('$\bar{S},[-]$','Interpreter','Latex')
grid on

figure
% subplot (3,1,3)
b = bar(table2array(S_soilT(1:20,:)),'FaceColor','flat');
for i=1:20 
    for ii=1:4
        if table2array(S_soilT(i,ii))<0
            newColor(ii,:) = negativeColor(ii,:);
        else
            newColor(ii,:) = positiveColor(ii,:);
        end
    end
    for k=1:4
        b(k).CData(i,:) = newColor(k,:);
    end
end
title('Sensitivity Analysis: Soil Temperature')
xticks([1:1:290])
xticklabels(S_soilT.Properties.RowNames)
xtickangle(45)
ylim([-0.15 0.15])
ylabel('$\bar{S},[-]$','Interpreter','Latex')
grid on

%% Plot the sensitivity change over a day for each parameter
figure
set(gcf, 'units','centimeters','Position',  [0, 0, 18.3, 11.9])

names = S_temperature.Properties.RowNames;
ParamNames=readtable('ParamNames.xlsx');
greyColor =[203 221 235]./255;
% for i=1:50
% S_temperature.Properties.RowNames{i} = char(ParamNames.ParamName(i));
% end
yyaxis left
for i=1:361 %number of parameters (half)
%     subplot(5,10,i)
    temp =S.(names{i}); % temporal variable with sensitivities from parameter i
    for j=1:4 % relative changes
        plot(TimeData, temp(:,1,j),'-','Color',greyColor, 'LineWidth',0.8)%,...
%             [1:97], ones(97,1).*table2array(S_temperature(i,j)),'-k') %sensitivities for temperature
        hold on
%         xlabel(char(ParamNames.ParamName(i)))
    end
end
ylim([-0.5 0.5])
ylabel('$\bar{S},[-]$','Interpreter','Latex','FontSize', 14)
datetick('x','HH','keepticks')

yyaxis right
plot(TimeData,x_nom(:,1),'-','Color', orangeColor,'LineWidth',1)
ylim([15 32])
ylabel('T_{Air}, [°C]')
datetick('x','HH','keepticks')
xlabel('Time, [h]')
% ylim([-0.5 0.5])
%% Example state, sensitivity over a day
startDate = datenum( '7-Oct-2019 00:00:00');
endDate = datenum( '7-Oct-2019 23:59:00');
TimeData = linspace(startDate,endDate,97);
orangeColor = [210 109 25]./255;
blueColor =[62 117 166]./255;
vpsat_s      = vapourPsat(x_nom(:,1));
vp_s         = vpsat_s.*x_nom(:,2);
vpsat_m      = vapourPsat(xreal(:,1));
vp_m         = vpsat_m.*xreal(:,2);

figure
set(gcf, 'units','centimeters','Position',  [0, 0, 18.3, 8.9])

subplot(2,1,1)%Temperature state over a day
plot(TimeData,xreal(:,1),':k','LineWidth',1.5); hold on
plot(TimeData,x_nom(:,1),'-','Color', orangeColor,'LineWidth',1)
legend('Real', 'Simulated', 'Location','eastoutside')
xlabel('Time, [h]')
ylabel('T_{air}, [°C]')
% xlim([0 97])
datetick('x','HH','keepticks')
title('Temperature state')

% subplot(2,2,2)%RH state over a day
% plot(TimeData,vp_m,'-k',...
%     TimeData,vp_s,'--r')
% legend('real', 'simulated')
% xlabel('Time, [h]')
% ylabel('VP_{air}, [kPa]')
% datetick('x','HH','keepticks')
% title('VP state')
%-------------------------------------------------------------------
subplot(2,1,2)%SA_temperature for Cd state over a day
temp = S.(names{1});
plot(TimeData,temp(:,1,3),'-','Color', blueColor,'LineWidth',1); hold on% Temperature, 20% increase in parameter
temp = S.(names{2});
plot(TimeData,temp(:,1,3),'--','Color', blueColor,'LineWidth',1); hold on% Temperature, 20% increase in parameter
temp = S.(names{5});
plot(TimeData,temp(:,1,3),'-.','Color', blueColor,'LineWidth',1); hold on% Temperature, 20% increase in parameter
temp = S.(names{7});
plot(TimeData,temp(:,1,3),'.','Color', blueColor,'LineWidth',1); hold on% Temperature, 20% increase in parameter

legend([char(ParamNames.ParamName([1,2,5,7]))], 'Location','eastoutside')

% for j=1:4 % relative changes
%         plot(TimeData, temp(:,1,j),'-')
%         hold on
% end
xlabel('Time, [h]')
ylabel('$\bar{S},[-]$','Interpreter','Latex','FontSize', 14)
ylim([-0.5 0.5])
datetick('x','HH','keepticks')
title('Sensitivities')

% subplot(2,2,4)%SA_VP for Cd state over a day
% temp =S.(names{2});
% for j=1:4 % relative changes
%         plot(TimeData, temp(:,2,j),'-')
%         hold on
% end
% xlabel('Time, [h]')
% ylabel('$\bar{S_{Cd}},[-]$','Interpreter','Latex','FontSize', 14)
% ylim([-0.35 0.35])
% datetick('x','HH','keepticks')
% title('Sensitivity to Cd: VP')

