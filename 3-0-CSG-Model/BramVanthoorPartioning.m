%% Doxygen Documentation
%> @file BramVanthoorPartitioning.m
%======================================================================
%> @brief Bram Vanthoor Assimilate Partitioning Model
%>
%> This block describes the Bram Vanthoor *Assimilate Partitioning Model*.
%>
%> @param LAI, Tcan, Cbuf, Tcan24, Tsum <br><br>
%> <table>
%> <tr><th> Input Name <th> Category <th> Type <th> Unit
%> <tr><td> \b LAI     <td> ConnectionCategory::LeafAreaIndex <td> ConnectionType::LeafAreaIndex_m2_per_m2 <td> \f$\left[m^{2}\;[Leaf]\;m^{-2}\;[ground]\right]\f$
%> <tr><td> \b Tcan    <td> ConnectionCategory::Climate <td> ConnectionType::Temperature_degrees <td> \f$\left[ ^\circ C \right]\f$
%> <tr><td> \b Cbuf    <td> <td> StateType::AssimilateBuffer <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\right]\f$
%> <tr><td> \b Tcan24  <td> <td> StateType::Temperature24h <td> \f$\left[ ^\circ C \right]\f$
%> <tr><td> \b Tsum    <td> <td> StateType::TemperatureSum <td> \f$\left[ ^\circ C \right]\f$
%> </table>
%>
%> @retval MCbuffruit_tot, MCbufleaf, MCbuf_stemroot, Tcan24, Tsum <br><br>
%> <table>
%> <tr><th> Output Name <th> Category <th> Type <th> Unit
%> <tr><td> \b MCbuffruit_tot <td> ConnectionCategory::DMFlow <td> ConnectionType::Fruit_Flow     <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\;s^{-1}\right]\f$
%> <tr><td> \b MCbufleaf      <td> ConnectionCategory::DMFlow <td> ConnectionType::Leaf_Flow      <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\;s^{-1}\right]\f$
%> <tr><td> \b MCbuf_stemroot <td> ConnectionCategory::DMFlow <td> ConnectionType::StemRoots_Flow <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\;s^{-1}\right]\f$
%> <tr><td> \b Tcan24  <td> <td> StateType::Temperature24h <td> \f$\left[ ^\circ C \right]\f$
%> <tr><td> \b Tsum    <td> <td> StateType::TemperatureSum <td> \f$\left[ ^\circ C \right]\f$
%> </table>
%>
%> ### Unit Conversions
%> The Framework works with \f$ kg \f$ instead of \f$ mg \f$.
%======================================================================
function [MCbuffruit_tot,MCbufleaf,MCbuf_stemroot,Tcan24,Tsum] = BramVanthoorPartioning(LAI, Tcan, Cbuf, Tcan24, Tsum)
	global rg_stemroot rg_leaf rg_fruit Tsum_needed T24_S1 T24_b1 T24_S2 T24_b2 Tinst_S1 Tinst_b1 Tinst_S2 Tinst_b2 LAI_Max 
    
    %% Calculations
    Cbuf_max = min(20e3,20e3*LAI/LAI_Max);   %IliasTsafaras, 11/12/2013
    % Updated to avoid model that Cbuf becomes negetaive, Sept 24 October
    h_Buf_MCBufOrg = SmoothIfElse(Cbuf, 0.05*Cbuf_max, -20e-3);     % Buffer empty
    % Temperature inhibition values by temperature filter
    % Describe smoothed, see page 11 PhD van Oothegem by adding 2 phi smoothed functions of Van Oothegem
    % smoothed version Oothegem page 12 and test_smth_func.m
    f1 = 1/Tsum_needed*Tsum;
    f2 = 0;
    f3 = 1/Tsum_needed.*(Tsum - Tsum_needed);

    % Smoothed versions
    a1_s = 0.5*(f2 + f1 + sqrt((abs(f2 - f1)).^2 +1e-4));
    a2_s = 0.5*(f2 + f3 + sqrt((abs(f2 - f3)).^2 +1e-4));
    g_MCBufFruit_Tsum = a1_s - a2_s;
    
    [h_Tcan_Char,h_Tcan24_Char,~,~,~,~,~,~] = GrowthInhibition2(Tcan,Tcan24,T24_S1,T24_b1,T24_S2,T24_b2,Tinst_S1,Tinst_b1,Tinst_S2,Tinst_b2);
    % Temperature dependent growth for all organs
    g_MCBufOrg_Tcan24 = 0.047*Tcan24 + 0.060;                       % - Based upon de Koning page 30, truss appereance is related to fruit growth
      
    MCbuffruit_tot =   h_Buf_MCBufOrg * h_Tcan24_Char * h_Tcan_Char * g_MCBufOrg_Tcan24 * g_MCBufFruit_Tsum * rg_fruit;
    MCbufleaf =        h_Buf_MCBufOrg * h_Tcan24_Char *               g_MCBufOrg_Tcan24 *                     rg_leaf ;
    MCbuf_stemroot =   h_Buf_MCBufOrg * h_Tcan24_Char *               g_MCBufOrg_Tcan24 *                     rg_stemroot;    
    
end

