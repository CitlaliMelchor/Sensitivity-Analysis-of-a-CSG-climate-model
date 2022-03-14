function dx = SemiPGBo( t, x, U0,Out, indoorCO2, spg, l_spg, ini_time, ini_location, Angle_spg, angle1,length1,...
    n_cover, n_canopy, n_floor, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_blanketi, n_blankete,n_soil_last, vf_mat, Capd_air, ...
    Atop_vent, Aside_vent, htop, hside, Vair,shape_cover)
%The semi-passive greenhouse model and plant model

%% Run the greenhouse model

global T_soilbound Water_con 
%% outside conditation %%
[climate, C_CO2] = climatedataBo( Out, indoorCO2, t );    % climate = [1-T_sky, 2- T_out, 3- Rad_out, 4- V_wind, 5- RH_out]

%% temperature state
Tair  = x(1);
RHair = x(2);
LAI   = x(3);
x_spg = x(4:(2+l_spg));
x_plant = x((3+l_spg):end);         % x_plant = [Tsum, Tcan24, Cbuf, Cstemroot, Cleaf, Char, Cfruit, Nfruit];
x_spg(l_spg) = climate(1);

%% time 
hour_clock=ini_time(1)+t/3600-floor((ini_time(1)+t/3600)/24)*24; %indicate the time of one day in 24 hours  [hour]
day_num = daynumber(t,[ini_time(2),ini_time(3)],ini_time(1));

%% control 
[u_blanket,u_vent, u_ventside, u_venttop, u_venttopside] = ControlventandblanketBo(U0,day_num,hour_clock);

%% radiation transmission
[tran,spg,angle_solar] = tranmatrixBo(t,day_num,ini_time,ini_location,Angle_spg,shape_cover,spg,u_blanket,LAI,n_cover,n_canopy,n_floor);

%% Radiative heat flux
%1. global light absorption
[Glob_R , Rad_canopy,PAR_canopy] = GlobalRexchangeBo( climate(3), angle_solar(1),tran, angle1,length1, spg.area, spg.label1,spg.label2,spg.label3, n_cover, n_canopy, n_screen, n_blanketi, n_blankete, u_blanket, l_spg,spg.outlayer, LAI);
%2. longwave radiation exchange
FIR = FIRexchangeBo( vf_mat, x_spg,n_screen,n_canopy,n_cover,n_blanketi,n_blankete,n_floor,u_blanket,spg );

%% heat convective and conductive and ventilation
sens_cond = conductiveBo( x_spg, spg.thickness, spg.area, spg.heatcond, spg.labelcond, T_soilbound ,n_cover, n_blanketi, n_blankete, u_blanket, n_soil_last);     %heat conductive inside layers [W per meter length greenhouse]         
[sens_conv , sens_air ] = convectiveBo( Tair, climate(2), spg.area, x_spg, climate(4), l_spg, n_floor, n_canopy, n_walli, n_walle, n_roofi, n_roofe, n_screen, n_cover, n_blankete, u_blanket );
[Vent,H_airvent] = ventilationBo( Atop_vent, Aside_vent, htop, hside, spg.area(n_floor), climate(4), Tair, climate(2), u_vent, u_venttopside, u_ventside, u_venttop, u_blanket);

%% latent heat exchange
[~,dRH,~,H_lat] = vptranBo(Vent, Tair, x_spg(n_cover), x_spg(n_canopy), x_spg(n_floor+2), climate(2), RHair, climate(5), spg.area(n_cover), Rad_canopy, C_CO2, Capd_air, Vair, LAI, spg.area(n_floor), Water_con,n_cover,n_floor,n_canopy,l_spg);

%% the capacity * density * area
cap = [spg.capacity, spg.density, spg.thickness, spg.area,1000*ones(l_spg,1)];
Capd_spg = prod(cap,2);
% Capd_spg1 = spg.capacity.*spg.density.*spg.thickness.*spg.area.*1000;             % [J/K]

%% Temperature change
Heat_airdot = sens_air-H_airvent + 0.7*0.05*climate(3);     % + heating system
Heat_spgdot = Glob_R + FIR + sens_cond + sens_conv + H_lat;
T_airdot = Heat_airdot./Capd_air;
T_spgdot = zerodivided0(Heat_spgdot,Capd_spg);

%% Run plant model -------Bram Vanthoor
global NrBox
    % x_plant = [Tsum; Tcan24; Cbuf; Cstemroot; Cleaf; Char; Cfruit; Nfruit];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PhotoSynthesis Model
    [MCairbuf] = BramVanthoorPhotoSynthesis(x_spg(n_canopy), PAR_canopy, C_CO2, LAI, x_plant(3));
    % Assimilate Partitioning
    [MCbuffruit_tot,MCbufleaf,MCbuf_stemroot] = BramVanthoorPartioning(LAI, x_spg(n_canopy), x_plant(3), x_plant(2), x_plant(1));
    % Growth Respiration
    [MCbufair_g] = BramVanthoorGrowthRespiration(MCbuffruit_tot,MCbufleaf,MCbuf_stemroot,LAI,x_plant(3));
    % CarboHydrate Buffer 
    [Cbufdot] = BramVanthoorCarboBuffer(MCairbuf, MCbuffruit_tot, MCbufleaf, MCbuf_stemroot, MCbufair_g);
    % Maintenance Respiration
    [MCstemroot_air_m, MCleafair_m, MCfruitair_m] = BramVanthoorMaintenanceRespiration(MCbufleaf,x_spg(n_canopy),x_plant(4),x_plant(5),x_plant(7:6+NrBox));
    % Leaf Harvest
    [MCleafhar] = BramVanthoorLeafHarvest(x_plant(5));
    % Plant Buffers
    [Chardot,LAIdot,Cstemrootdot,Cleafdot,Cfruitdot,Nfruitdot] = BramVanthoorPlantBuffer(MCbuffruit_tot,MCbufleaf,MCbuf_stemroot,MCstemroot_air_m, MCleafair_m, MCfruitair_m,MCleafhar,x_spg(n_canopy),x_plant(1),x_plant(2),x_plant(7:(6+NrBox)),x_plant((7+NrBox):(6+2*NrBox)));
    % Tsum & Tcan24
    Tsumdot    = x_spg(n_canopy)/86400;
    Tcan24dot  = (x_spg(n_canopy)-x_plant(2))/86400;

%% dot------------------

dx = [T_airdot; dRH; LAIdot; T_spgdot(1:(l_spg-1));Tsumdot; Tcan24dot; Cbufdot;Cstemrootdot;Cleafdot;Chardot;Cfruitdot';Nfruitdot'];

end

