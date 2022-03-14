%% Doxygen Documentation
%> @file BramVanthoorGrowthRespiration.m
%======================================================================
%> @brief Bram Vanthoor Growth Respiration
%>
%> This block describes the Bram Vanthoor *Growth Respiration Model*.
%>
%> @param MCbuffruit_tot, MCbufleaf, MCbuf_stemroot, LAI, Cbuf <br><br>
%> <table>
%> <tr><th> Input Name <th> Category <th> Type <th> Unit
%> <tr><td> \b MCbuffruit_tot <td> ConnectionCategory::DMFlow <td> ConnectionType::Fruit_Flow <td> \f$\left[m^{2}\;[Leaf]\;m^{-2}\;[ground]\right]\f$
%> <tr><td> \b MCbufleaf <td> ConnectionCategory::DMFlow <td> ConnectionType::Leaf_Flow <td> \f$\left[m^{2}\;[Leaf]\;m^{-2}\;[ground]\right]\f$
%> <tr><td> \b MCbuf_stemroot <td> ConnectionCategory::DMFlow <td> ConnectionType::StemRoots_Flow <td> \f$\left[m^{2}\;[Leaf]\;m^{-2}\;[ground]\right]\f$
%> <tr><td> \b LAI <td> ConnectionCategory::LeafAreaIndex <td> ConnectionType::LeafAreaIndex_m2_per_m2 <td> \f$\left[m^{2}\;[Leaf]\;m^{-2}\;[ground]\right]\f$
%> <tr><td> \b Cbuf <td> <td> StateType::AssimilateBuffer <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\right]\f$
%> </table>
%>
%> @retval MCbufair_g <br><br>
%> <table>
%> <tr><th> Output Name <th> Category <th> Type <th> Unit
%> <tr><td> \b MCbufair_g <td> ConnectionCategory::GrowthResp <td> ConnectionType::Growth_Resp_kg_CH2O_per_m2_per_s <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\;s^{-1}\right]\f$
%> </table>
%>
%> ### Unit Conversions
%> The Framework works with \f$ kg \f$ instead of \f$ mg \f$.
%======================================================================
function [MCbufair_g] = BramVanthoorGrowthRespiration(MCbuffruit_tot,MCbufleaf,MCbuf_stemroot,LAI,Cbuf)
	global cfruit_g cleaf_g cstemroot_g LAI_Max
    %% Calculations
    Cbuf_max = min(20e3,20e3*LAI/LAI_Max);   %IliasTsafaras, 11/12/2013

    % Updated to avoid model that Cbuf becomes negetaive, Sept 24 October
    h_Buf_MCBufOrg = SmoothIfElse(Cbuf, 0.05*Cbuf_max, -20e-3);     % Buffer empty
    
    % Growth respiration based upon Heuvelink Phd, page 238
    MCfruitair_g     = cfruit_g    * MCbuffruit_tot;
    MCleafair_g      = cleaf_g     * MCbufleaf;
    MCstemroot_air_g = cstemroot_g * MCbuf_stemroot;
    
    % Total growth respiration
    MCbufair_g = MCfruitair_g + MCleafair_g + MCstemroot_air_g;
    MCbufair_g = h_Buf_MCBufOrg*MCbufair_g;         % When buffer is empty, no growth respiration
    
end

