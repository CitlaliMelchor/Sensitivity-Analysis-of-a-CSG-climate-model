%% Doxygen Documentation
%> @file BramVanthoorPlantBuffer.m
%======================================================================
%> @brief Bram Vanthoor Plant Buffer Model
%>
%> This block describes the Bram Vanthoor *Plamt Buffer Model*.
%>
%> @param MCbuffruit_tot, MCbufleaf, MCbuf_stemroot, MCstemroot_air_m, MCleafair_m, MCfruitair_m, MCleafhar, Tsum, Tcan, Tcan24, Cstemroot, Cleaf, Cfruit, Nfruit <br><br>
%> <table>
%> <tr><th> Input Name <th> Category <th> Type <th> Unit
%> <tr><td> \b MCbuffruit_tot   <td> ConnectionCategory::DMFlow           <td> ConnectionType::Fruit_Flow     <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\;s^{-1}\right]\f$
%> <tr><td> \b MCbufleaf        <td> ConnectionCategory::DMFlow           <td> ConnectionType::Leaf_Flow      <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\;s^{-1}\right]\f$
%> <tr><td> \b MCbuf_stemroot   <td> ConnectionCategory::DMFlow           <td> ConnectionType::StemRoots_Flow <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\;s^{-1}\right]\f$
%> <tr><td> \b MCstemroot_air_m <td> ConnectionCategory::MaintenanceResp  <td> ConnectionType::Maint_Resp_StemRoot_kg_CH2O_per_m2_per_s <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\;s^{-1}\right]\f$
%> <tr><td> \b MCleafair_m      <td> ConnectionCategory::MaintenanceResp  <td> ConnectionType::Maint_Resp_Leaf_kg_CH2O_per_m2_per_s <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\;s^{-1}\right]\f$
%> <tr><td> \b MCfruitair_m     <td> ConnectionCategory::MaintenanceResp  <td> ConnectionType::Maint_Resp_Fruit_kg_CH2O_per_m2_per_s <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\;s^{-1}\right]\f$
%> <tr><td> \b MCleafhar        <td> ConnectionCategory::HarvestedLeaf    <td> ConnectionType::Har_Leaf_kg_CH2O_per_m2_per_s           <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\;s^{-1}\right]\f$
%> <tr><td> \b Tcan             <td> ConnectionCategory::Climate          <td> ConnectionType::Temperature_degrees <td> \f$\left[ ^\circ C \right]\f$
%> <tr><td> \b Tsum             <td> <td> StateType::TemperatureSum       <td> \f$\left[ ^\circ C \right]\f$
%> <tr><td> \b Tcan24           <td> <td> StateType::Temperature24h       <td> \f$\left[ ^\circ C \right]\f$
%> <tr><td> \b Cstemroot        <td> <td> StateType::StemRootBiomass      <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\right]\f$
%> <tr><td> \b Cleaf            <td> <td> StateType::LeafBiomass          <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\right]\f$
%> <tr><td> \b Cfruit           <td> <td> StateType::FruitBiomass_per_Dev <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\right]\f$
%> <tr><td> \b Nfruit           <td> <td> StateType::FruitNumbers_per_Dev <td> \f$\left[-\right]\f$
%> </table>
%>
%> @retval MCfruit, LAI, Tsum, Tcan24, Cstemroot, Cleaf, Cfruit, Nfruit <br><br>
%> <table>
%> <tr><th> Output Name <th> Category <th> Type <th> Unit
%> <tr><td> \b MCfruit          <td> ConnectionCategory::DMFlow           <td> ConnectionType::TemperatureSum       <td> \f$\left[ ^\circ C \right]\f$
%> <tr><td> \b LAI              <td> ConnectionCategory::LeafAreaIndex    <td> ConnectionType::LeafAreaIndex_m2_per_m2       <td> \f$\left[ ^\circ C \right]\f$
%> <tr><td> \b Tsum             <td> <td> StateType::TemperatureSum       <td> \f$\left[ ^\circ C \right]\f$
%> <tr><td> \b Tcan24           <td> <td> StateType::Temperature24h       <td> \f$\left[ ^\circ C \right]\f$
%> <tr><td> \b Cstemroot        <td> <td> StateType::StemRootBiomass      <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\right]\f$
%> <tr><td> \b Cleaf            <td> <td> StateType::LeafBiomass          <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\right]\f$
%> <tr><td> \b Cfruit           <td> <td> StateType::FruitBiomass_per_Dev <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\right]\f$
%> <tr><td> \b Nfruit           <td> <td> StateType::FruitNumbers_per_Dev <td> \f$\left[-\right]\f$
%> </table>
%>
%> ### Unit Conversions
%> The Framework works with \f$ kg \f$ instead of \f$ mg \f$.
%======================================================================
function [Chardot,LAIdot,Cstemrootdot,Cleafdot,Cfruitdot,Nfruitdot] = BramVanthoorPlantBuffer(MCbuffruit_tot,MCbufleaf,MCbuf_stemroot,MCstemroot_air_m, MCleafair_m, MCfruitair_m,MCleafhar,Tcan,Tsum,Tcan24,Cfruit,Nfruit)
	global NrBox SLA 
    %% Calculations        
    % Fruit development rate based upon the Koning
    r_dev = (-0.066+0.1*Tcan24)./(100*86400);        % Development rate, s-1 De Koning 56, linear approach because others don't work outside the range 17 - 27°C
    FGP = (1/(r_dev*86400));        % Fruit growth period
    
    % ------------------------------Maintenace respiraton of plant organs--------------------
    % Total maintenace respiration depends on "relative egrowth rate" see Heuvelink 239
%     f_RGR = 33*86400;
%     RGR = (MCbuffruit_tot + MCbufleaf + MCbuf_stemroot)/(sum(Cfruit) + Cleaf + Cstemroot);
      
    % Determine if a specific development stage contains fruits to judge if they receive assimilates from buffer. Fruits available, h_Cfruit_MCbuffruit == 1
    t_Kon=zeros(1,NrBox);
    for i = 1:NrBox   
        % Determine de Koning time in order to determine Gompertz equation (page 105 PhD De Koning)
        t_Kon(i) = (2*(i-1)+1)/(2*NrBox) * FGP;
    end
     
    % When assimlates flow to the fruit compartment (after a specific development stage) than the first development stage is always open to receive assimilates
    h_Tcansum_Cfruit = ones(NrBox,1);           % This inhibition is not needed anymore!
    
    % Determine relative fruit growth rate based upon Gompertz function page 105, 107 de Koning
    M =  -4.93 + 0.548*FGP;     % Inflexition point, asymptotical around this value, close to half of the FGP
    C = 10;                      % Determines maximum fruit yield
    B = 1/(2.44+0.403*M);
    
    % Growth rate based upon Gompertz, GR is a vector               
    GR = C*exp(-exp(-B*(t_Kon-M))).*B.*exp(-B*(t_Kon-M));           % g{DM}.day-1.fruit-1   equation 4.4.5 equation PhD de Koning 
    
    % Potential fruit weight based upon Gompertz
    Wfruit_pot = C.*exp(-exp(-B*(t_Kon-M)))*1000;        % mg{CH20}.N{fruit}-1
      
    % Determine the fruit number flow, also inhibited by existence of fruit development stages

    % Last part later introduced, needed if Tsum is low!
    MNfruit = r_dev.*Nfruit.*NrBox.*SmoothIfElse(Tsum, 0, -1/20);        % Fruit flow from dev i-1 to dev i
    
    % Fruit flow to development stage, carbohydrates comes from buffer
     MNfruit1max = (-1.708e-7 + 7.3125e-7*Tcan)*2.5;               % Maximum fruit set  N{fruits}.s-1, based upon De Koning page 36, with a 9 fruits/truss page 45
%      MCpot = 0.1;                                            % Above this assimilates supply the maximum fruit set is reached, under this value a linear decrease!
%      MNfruit1pot = MCbuffruit_tot* MNfruit1max/MCpot;                                 % Potential fruit set based upon available assimilates  
%      MNbuffruit = min(MNfruit1pot , MNfruit1max);                                % Potential fruit set can not exceed maximum fruit set
     
     % Introduced 8 January 2009 a smooth version
     MNbuffruit = SmoothIfElse(MCbuffruit_tot, 0.05, -58.9)*MNfruit1max;
     
     %Smooth Switch for 
     MNbuffruit = SmoothIfElse(Tsum, 0, -1/20)*MNbuffruit;
     
     % Determine potential carbohydrate flow when enough carbohydrates are available
     MCbuffruit(1) = MNbuffruit * Wfruit_pot(1);
    
    % Check if all carbohydrates needed for fruit set can be delivered
%     if (MNfruit1max/MCpot) > 1/Wfruit_pot(1)
%        display('Not Enough carbohydrates available for fruit set!!!') 
%     end
        
    % Determine relative growth factor to assure that all assimalates available for fruit growth flow to the fruits
%     GR_tot = sum(GR(2:end));
    logic1 = sum(h_Tcansum_Cfruit(2:end)) == 0 || sum(Nfruit(2:end)) == 0 || sum(GR(2:end))==0;
    rel_factor = 1/(sum( GR(2:end).*Nfruit(2:end)'))*switch02(logic1,1);

    % Determine all buffer flows, Assimilates a) from dev "i" to dev "i + 1", b) from Buf to dev "i", c) maintenace loss from dev "i"
    MCfruit=zeros(1,NrBox);
    % Increase Photosynthesis
    for i = 1: NrBox      
        MCfruit(i) = r_dev * NrBox*Cfruit(i); 
    end
 
    
    % Determine the fraction of carbohydrates that flow from the buffer to a specific development stage
    eta_Buf_fruit = zeros(1,NrBox);
    MCbuffruit(1,2:NrBox) = zeros(1,NrBox-1);
    for i = 2: NrBox   
        % Fraction of available assimilates that flow to a specific fruit developments stage
        eta_Buf_fruit(i) = GR(i)*Nfruit(i)*rel_factor;
        MCbuffruit(i) = eta_Buf_fruit(i)* (MCbuffruit_tot -MCbuffruit(1)) ;
    end
    
    % Stem and Root rate                
    Cstemrootdot = MCbuf_stemroot - MCstemroot_air_m;
    
    % Leaf rate
    Cleafdot     = MCbufleaf - MCleafair_m - MCleafhar;          % Leaf dry weight mg CH20 /m2
   
    % Calculate Leaf Area Index rate
    LAIdot          = SLA*Cleafdot;
    
    % Fruit rate
    Cfruitdot    = zeros(1,NrBox);
    Nfruitdot    = zeros(1,NrBox);
    Cfruitdot(1) = MCbuffruit(1) - MCfruit(1) - MCfruitair_m(1); % Fruit dry weight mg CH20 /m2
    Nfruitdot(1) = MNbuffruit - MNfruit(1);                      % Change of fruit number in fruit development stage 1 N{fruits}.s-1
    
    % For each fruit development stage the number of fruits and the carbon content is determined
    for i = 2:NrBox        
        Cfruitdot(i) = MCbuffruit(i) + MCfruit(i-1) - MCfruit(i) - MCfruitair_m(i);    % Carbon content of all fruit stages
        Nfruitdot(i) = MNfruit(i-1) - MNfruit(i);       % N{fruits}.s-1
    end
    Chardot = MCfruit(end);  % Final crop harvest, mg{DM} m-2

end

