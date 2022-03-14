%% Variance-based Global SA - latin hypercube sampling
% MSc Thesis Sensitivity analysis of a Chinese Solar Greenhouse
% Author: Citlali Melchor Ramírez
% December 2021

clear 
clc
%% Read the bounds table, with the lower and upper bounds of each parameter
boundsTable = readtable('ParamBoundsFive.xlsx');
bounds = table2array(boundsTable(:,3:4));
%% initial setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 parametersBo
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
% initial input value
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

% location
latitude    = 39.62;                                        %local latitude -Beijing is 40 oN                                      [deg]
longitude   = 115.97;                                       %local longitude -Beijing is 117 oN                                    [deg]
Hourcircle  = 120;                                          %hour circle applied in the experiment area -Beijing is 120 oN         [deg] 
ini_location = [latitude, longitude, Hourcircle];

% %%%%%%%%%%%%%Run the model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------------------------------%%%%%%%%%%%%%
%% simulation time
day_num_s1 = daynumber(0,[set_month,set_day],set_hour);
ini_time = [set_hour, set_month, set_day];
%input data
xreal(:,5)= round(xreal(:,5)./100).*100;
% time distance
dis_time = 800;
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

% ode 45
x0 = [Tair0,RH0,LAI0,spg.iniTem(1:(l_spg-1))',x_plant0]';
tspan=[0:dis_time:(n1-1)*dis_time]' ;
%%
[t_Bo, x_Bo]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
    n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
    Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
%  

%% Variance-based method
% For this method, only N = 100 iterations per matrix are considered, for
% the k = 50 parameters
N = 1000;
k = 5;
%%
w = size(x_Bo,1);
d = size(x_Bo,2);
y_A = zeros(N,w,d);
y_B = zeros(N,w,d);
y_C = zeros(N,w,d,k);
y_D = zeros(N,w,d,k);

%%
% 1. Define a matrix A with size (N,k) and evaluate the model y(A)
A = SA_LatHyp(k,bounds,N);
%%
for i=1:N
    Cd              = A(i,1);
    r_net           = A(i,2);
    spg.NIRtran(18) = A(i,3);
    n               = A(i,4);
    T_min_day       = A(i,5);

    [t_A,x_A]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
            n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
            Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    y_A(i,:,:) = x_A;%(48,1); % noon, temperature
%     y_A_night(i,:) = x_A(20,1); % night, temperature
%     yT_A0(i,:) = x_A(end,1); % midnight, temperature
%     yRH_A12(i,:) = x_A(48,2); % noon, temperature
%     yRH_A0(i,:) = x_A(end,2); % midnight, temperature

end
fprintf('y_A is done')

save GlobalSA_varMethod_latin5v2_A.mat A y_A 
%%
% 2. Define a matrix B with size (N,k) and evaluate the model y(B)
B = SA_LatHyp(k,bounds,N);

for i=1:N
    Cd              = B(i,1);
    r_net           = B(i,2);
    spg.NIRtran(18) = B(i,3);
    n               = B(i,4);
    T_min_day       = B(i,5);
    
    [t_B,x_B]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
            n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
            Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
    y_B(i,:,:) = x_B;%(48,1); % noon, temperature
%     y_B_night(i,:) = x_B(20,1); % dusk, temperature
end
fprintf('y_B is done')
%%
save GlobalSA_varMethod_latin5v2_B.mat B y_B 
%%
load GlobalSA_varMethod_latin5v2_A.mat
load GlobalSA_varMethod_latin5v2_B.mat
% 3. Define a C_i matrix with the columns of B and one i column of A  and
% evaluate the model y(C_i)
for i = 1:k
    C = B;
    C(:,i) = A(:,i);
    % 4. Evaluate the model 
    for j=1:N
        Cd              = C(j,1);
        r_net           = C(j,2);
        spg.NIRtran(18) = C(j,3);
        n               = C(j,4);
        T_min_day       = C(j,5);
        [t_C,x_C]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
                n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
                Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
        y_C(j,:,:,i) = x_C;%(48,1); % noon, temperature
%         y_C_night(j,i) = x_C(20,1); % noon, temperature
    end
end
fprintf('y_C is done')
save GlobalSA_varMethod_latin5v2_C.mat  y_C %S ST 
%% 
load GlobalSA_varMethod_latin5v2_A.mat
load GlobalSA_varMethod_latin5v2_B.mat
%4. Define a D_i matrix with the columns of A and one i column of B
for i = 1:k
    D = A;
    D(:,i) = B(:,i);
    % 4. Evaluate the model 
    for j=1:N
        Cd              = D(j,1);
        r_net           = D(j,2);
        spg.NIRtran(18) = D(j,3);
        n               = D(j,4);
        T_min_day       = D(j,5);
        [t_C,x_D]=ode45(@(t,x) SemiPGBo( t, x, U, Out, xreal, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
                n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
                Atop_vent, Aside_vent, htop, hside, Vair, shape_cover),tspan,x0);
        y_D(j,:,:,i) = x_D;%(48,1); % noon, temperature
%         y_D_night(j,i) = x_D(20,1);
        
    end
end
fprintf('y_D is done')
save GlobalSA_varMethod_latin5v2_D.mat y_D 
%%
% %% First order sensitivities, and total effects
% 
% % 5. Compute the first order sensitivities
% f_0sq = ((1/N)*sum(y_A))^2;
% FFF
% for i = 1:k 
%     S(i) = ((1/N)*(y_A'*y_C(:,i))- f_0sq)/((1/N)*sum(y_A.^2) - f_0sq);
% end
% % fprintf('S is done')
% 
% % 6. Compute the total effect indices 
% for i = 1:k 
%     ST(i) = 1-(((1/N)*(y_B'*y_C(:,i)) - f_0sq)/((1/N)*sum(y_A.^2) - f_0sq));
% end
% % fprintf('ST is done')

%%
%  save GlobalSA_varMethod_v5full_D.mat y_D y_D_night%S ST 
%% Save the results
% save GlobalSA_varMethod_v6D.mat A B y_A y_B y_C y_D y_A_night y_B_night y_C_night y_D_night%S ST 