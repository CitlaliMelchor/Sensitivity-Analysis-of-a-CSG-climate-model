%% Doxygen Documentation
%> @file BramVanthoorTomatoHarvest.m
%======================================================================
%> @brief Bram Vanthoor Tomato Harvest Model
%>
%> This block describes the Bram Vanthoor *Tomato Harvest Model*.
%>
%> @param MCfruit, Char <br><br>
%> <table>
%> <tr><th> Input Name <th> Category <th> Type <th> Unit
%> <tr><td> \b MCfruit <td> ConnectionCategory::HarvestedTomato <td> ConnectionType::Har_Tomato_kg_CH2O_per_m2_per_s (CHANGE!) <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\;s^{-1}\right]\f$
%> <tr><td> \b Char    <td> <td> StateType::HarvestedTomatoes <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\right]\f$
%> </table>
%>
%> @retval Char <br><br>
%> <table>
%> <tr><th> Output Name <th> Category <th> Type <th> Unit
%> <tr><td> \b Char    <td> <td> StateType::HarvestedTomatoes <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\right]\f$
%> </table>
%>
%> ### Unit Conversions
%> The Framework works with \f$ kg \f$ instead of \f$ mg \f$.
%======================================================================
function [Char] = BramVanthoorTomatoHarvest(MCfruit,Char)
    global Simulation
    
    %% Input Conversion
    MCfruit = MCfruit*(10^6);
    Char    = Char*(10^6);
    
    %% Calculations
	Char = Char + MCfruit(end)*Simulation.SimRes;  % Final crop harvest, mg{DM} m-2
    
    %% Output Conversion
    Char = Char/(10^6);
end

