%% parameter and data
tic
 parametersBo
%  spg.capacity(19) = 1.9;
toc
 tic
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
set_day     = 6;  % Day
set_hour    = 0;  % hour
set_daynum  = 5;   % day number from 1 to 132.

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
dis_time = 600;
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

%%
Runode45Bo1;

toc
%end
% distribution
