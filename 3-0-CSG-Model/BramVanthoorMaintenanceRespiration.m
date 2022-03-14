%% Doxygen Documentation
%> @file BramVanthoorMaintenanceRespiration.m
%======================================================================
%> @brief Bram Vanthoor Maintenance Respiration Model
%>
%> This block describes the Bram Vanthoor *Maintenance Respiration Model*.
%>
%> @param MCbufleaf, Tcan, Cstemroot, Cleaf, Cfruit <br><br>
%> <table>
%> <tr><th> Input Name <th> Category <th> Type <th> Unit
%> <tr><td> \b MCbufleaf <td> ConnectionCategory::DMFlow <td> ConnectionType::Leaf_Flow <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\;s^{-1}\right]\f$
%> <tr><td> \b Tcan <td> ConnectionCategory::Climate <td> ConnectionType::Temperature_degrees <td> \f$\left[ ^\circ C \right]\f$
%> <tr><td> \b Cstemroot <td> <td> StateType::StemRootBiomass <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\right]\f$
%> <tr><td> \b Cleaf <td> <td> StateType::LeafBiomass <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\right]\f$
%> <tr><td> \b Cfruit <td> <td> StateType::FruitBiomass <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\right]\f$
%> </table>
%>
%> @retval MCstemroot_air_m, MCleafair_m, MCfruitair_m <br><br>
%> <table>
%> <tr><th> Output Name <th> Category <th> Type <th> Unit
%> <tr><td> \b MCstemroot_air_m <td> ConnectionCategory::MaintenanceResp <td> ConnectionType::Maint_Resp_StemRoot_kg_CH2O_per_m2_per_s <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\;s^{-1}\right]\f$
%> <tr><td> \b MCleafair_m <td> ConnectionCategory::MaintenanceResp <td> ConnectionType::Maint_Resp_Leaf_kg_CH2O_per_m2_per_s <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\;s^{-1}\right]\f$
%> <tr><td> \b MCfruitair_m <td> ConnectionCategory::MaintenanceResp <td> ConnectionType::Maint_Resp_Fruit_kg_CH2O_per_m2_per_s <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\;s^{-1}\right]\f$
%> </table>
%>
%> ### Unit Conversions
%> The Framework works with \f$ kg \f$ instead of \f$ mg \f$.
%======================================================================
function [MCstemroot_air_m, MCleafair_m, MCfruitair_m] = BramVanthoorMaintenanceRespiration(MCbufleaf,Tcan,Cstemroot,Cleaf,Cfruit)
    global NrBox cfruit_m cleaf_m cstemroot_m Q10m Resp_fac

    %% Calculations    
    MCfruitair_m=zeros(1,NrBox);

    for i = 1:NrBox      
        MCfruitair_m_without = Respiration(Cfruit(i), Tcan, 0, cfruit_m, Q10m, 0);
        MCfruitair_m(i) = Resp_fac*MCfruitair_m_without;
    end
   
    % ---------------Total fruit maintenance respiration -------
    MCfruitair_m_tot = sum(MCfruitair_m);
    
    %MCleafair_m = Resp_fac*Respiration(Cleaf, Tcan, 0, cleaf_m, Q10m, 0);
    %Ilias Tsafaras 17/12/2013 (on air....)
    MCleafair_m = min(Resp_fac*Respiration(Cleaf, Tcan, 0, cleaf_m, Q10m, 0),MCbufleaf);
    MCstemroot_air_m = Resp_fac*Respiration(Cstemroot, Tcan, 0, cstemroot_m, Q10m, 0);
    
    % Total maintenance respiration
    %MCorgair_m = MCfruitair_m_tot + MCleafair_m + MCstemroot_air_m;
    
    % -----------------End of maintenace----------------%

end