%% Doxygen Documentation
%> @file BramVanthoorCarboBuffer.m
%======================================================================
%> @brief Bram Vanthoor Carbohydrate Buffer Model
%>
%> This block describes the Bram Vanthoor *Carbohydrate Buffer Model*.
%>
%> @param MCAirBuf, MCbuffruit_tot, MCbufleaf, MCbuf_stemroot, MCbufair_g, Cbuf <br><br>
%> <table>
%> <tr><th> Input Name <th> Category <th> Type <th> Unit
%> <tr><td> \b MCAirBuf <td> ConnectionCategory::Photosynthesis <td> ConnectionType::PhotoSynthesis_kg_CH2O_per_m2_per_s <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\;s^{-1}\right]\f$
%> <tr><td> \b MCbuffruit_tot <td> ConnectionCategory::DMFlow <td> ConnectionType::Fruit_Flow     <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\;s^{-1}\right]\f$
%> <tr><td> \b MCbufleaf      <td> ConnectionCategory::DMFlow <td> ConnectionType::Leaf_Flow      <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\;s^{-1}\right]\f$
%> <tr><td> \b MCbuf_stemroot <td> ConnectionCategory::DMFlow <td> ConnectionType::StemRoots_Flow <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\;s^{-1}\right]\f$
%> <tr><td> \b MCbufair_g     <td> ConnectionCategory::GrowthResp <td> ConnectionType::Growth_Resp_kg_CH2O_per_m2_per_s <td>  \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\;s^{-1}\right]\f$
%> <tr><td> \b Cbuf           <td> <td> StateType::AssimilateBuffer <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\right]\f$
%> </table>
%>
%> @retval Cbuf <br><br>
%> <table>
%> <tr><th> Output Name <th> Category <th> Type <th> Unit
%> <tr><td> \b Cbuf           <td> <td> StateType::AssimilateBuffer <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\right]\f$
%> </table>
%>
%> ### Unit Conversions
%> The Framework works with \f$ kg \f$ instead of \f$ mg \f$.
%======================================================================
function [Cbufdot] = BramVanthoorCarboBuffer(MCAirBuf, MCbuffruit_tot, MCbufleaf, MCbuf_stemroot, MCbufair_g)
    %% Calculations
    % Assimilate Buffer Update
    Cbufdot = MCAirBuf - MCbuffruit_tot - MCbufleaf - MCbuf_stemroot - MCbufair_g;     % Carbon buffer CH20 of the plant mg CH20 /m2
end

