%% Doxygen Documentation
%> @file BramVanthoorPhotoSynthesis.m
%======================================================================
%> @brief Bram Vanthoor Photosynthesis Model
%>
%> This block describes the Bram Vanthoor *Photosynthesis Model*.
%>
%> @param Tcan, R_PARcan, CO2air, LAI, Cbuf <br><br>
%> <table>
%> <tr><th> Input Name <th> Category <th> Type <th> Unit
%> <tr><td> \b Tcan <td> ConnectionCategory::Climate <td> ConnectionType::Temperature_degrees <td> \f$\left[ ^\circ C \right]\f$
%> <tr><td> \b R_PARcan <td> ConnectionCategory::Climate <td> ConnectionType::GlobalRadiation_J_per_m2 <td> \f$\left[J\;m^{-2}\;s^{-1}\right]\f$
%> <tr><td> \b CO2air <td> ConnectionCategory::Climate <td> ConnectionType::CO2_ppm <td> CO_{2} Concentration \f$\left[ ppm \right]\f$
%> <tr><td> \b LAI <td> ConnectionCategory::LeafAreaIndex <td> ConnectionType::LeafAreaIndex_m2_per_m2 <td> \f$\left[m^{2}\;[Leaf]\;m^{-2}\;[ground]\right]\f$
%> <tr><td> \b Cbuf <td> <td> StateType::AssimilateBuffer <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\right]\f$
%> </table>
%>
%> @retval MCAirBuf <br><br>
%> <table>
%> <tr><th> Output Name <th> Category <th> Type <th> Unit
%> <tr><td> \b MCAirBuf <td> ConnectionCategory::Photosynthesis <td> ConnectionType::PhotoSynthesis_kg_CH2O_per_m2_per_s <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\;s^{-1}\right]\f$
%> </table>
%>
%> ### Unit Conversions
%> The Framework works with \f$ kg \f$ instead of \f$ mg \f$.
%======================================================================
function [MCAirBuf] = BramVanthoorPhotoSynthesis(Tcan, R_PARcan, CO2air, LAI, Cbuf)
	global LAI_Max THETA eta_ppm_mgm3
    %% Calculations
    Cbuf_max = min(20e3,20e3*LAI/LAI_Max);   %IliasTsafaras, 11/12/2013
    % Detailed explanation goes here
    Jmax25 = 210*LAI; 
    % PhotoSynthesis
    [Pg_L,~] = Photosynthesis3(Tcan, R_PARcan, CO2air/eta_ppm_mgm3,Jmax25,THETA,LAI);       % micromol CO2.m-2.s-1
    % Potential photosynthesis thus no inhibition by temperature effects or buffer size
    MCairbuf_pot = 30e-3*Pg_L;            % Conversion from micromol co2 to mg CH2O.m-2.s-1
    % Inhibition function by buffer "fullniss" or empty or full
    %     h_Buf_MCBufOrg = SmoothIfElse(Cbuf, 0.05*Cbuf_max, -5e-3);     % Buffer empty
    % Updated to avoid model that Cbuf becomes negetaive, Sept 24 October
    h_Buf_MCAirBuf = SmoothIfElse(Cbuf, Cbuf_max, 5e-4);           % Buffer full
    % By Photo inhibition the photosynthesis decreases
    MCAirBuf    = h_Buf_MCAirBuf*MCairbuf_pot;          % mg{CH20 }

end

