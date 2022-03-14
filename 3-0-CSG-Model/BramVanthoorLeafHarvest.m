%% Doxygen Documentation
%> @file BramVanthoorLeafHarvest.m
%======================================================================
%> @brief Bram Vanthoor Leaf Harvest Model
%>
%> This block describes the Bram Vanthoor *Leaf Harvest Model*.
%>
%> @param Cleaf <br><br>
%> <table>
%> <tr><th> Input Name <th> Category <th> Type <th> Unit
%> <tr><td> \b Cleaf            <td> <td> StateType::LeafBiomass          <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\right]\f$
%> </table>
%>
%> @retval MCleafhar <br><br>
%> <table>
%> <tr><th> Output Name <th> Category <th> Type <th> Unit
%> <tr><td> \b MCleafhar <td> ConnectionCategory::HarvestedLeaf <td> ConnectionType::Har_Leaf_kg_CH2O_per_m2_per_s     <td> \f$\left[kg\;[CH_{2}O]\;m^{-2}\;[ground]\;s^{-1}\right]\f$
%> </table>
%>
%> ### Unit Conversions
%> The Framework works with \f$ kg \f$ instead of \f$ mg \f$.
%======================================================================
function [MCleafhar] = BramVanthoorLeafHarvest(Cleaf)
    global LAI_Max SLA

    %% Calculations
    % --------------------Fruit Harvest and leaf pruning-------------
    
    % Vegetative harvest, leaf pruning
    Cleaf_Max = LAI_Max/SLA ;
       
    % Continuous appraoch of "if-then-else structure", first value is to avoid strong discontinuty    
    MCleafhar =0.001*SmoothIfElse(Cleaf, Cleaf_Max, -5e-5)*(Cleaf - Cleaf_Max);        %0.001*
    MCleafhar = max(MCleafhar,0);

end