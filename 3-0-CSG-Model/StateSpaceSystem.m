function [StateOut,OutputFlows, LAIOut] = StateSpaceSystem(StatesIn,Input,Time, LAIIn)
%INTEGRATIONSTEP State Space Representation of the Tomato Crop Yield Model
%   This file implements the integration for the tomato crop yield model of
%   Bram Vanthoor. This file is a dressed out version containing only the
%   integration of the tomato crop yield model and not the greenhouse
%   model. 
global NrBox
%% Parse State Vector
    Cbuf         = StatesIn(1);
    Cstemroot    = StatesIn(2);
    Cleaf        = StatesIn(3);
    Cfruit       = StatesIn(4:4+(NrBox-1));
    Nfruit       = StatesIn(4+(NrBox-1)+1:4+2*(NrBox-1)+1);
    Char         = StatesIn(4+2*(NrBox-1)+2);
    Tsum_1       = StatesIn(4+2*(NrBox-1)+3);
    Tcan24_1     = StatesIn(4+2*(NrBox-1)+4);
    Tsum_2       = StatesIn(4+2*(NrBox-1)+5);
    Tcan24_2     = StatesIn(4+2*(NrBox-1)+6);
%% Parse Input Vector
    Tcan         = Input(1);
    R_PARcan     = Input(2);
    CO2air       = Input(3);  
%% Run the Model
    % PhotoSynthesis Model
    [MCairbuf] = BramVanthoorPhotoSynthesis(Tcan, R_PARcan, CO2air, LAIIn, Cbuf);
    % Assimilate Partitioning
    [MCbuffruit_tot,MCbufleaf,MCbuf_stemroot,Tcan24_1,Tsum_1] = BramVanthoorPartioning(LAIIn, Tcan, Cbuf, Tcan24_1, Tsum_1);
    % Growth Respiration
    [MCbufair_g] = BramVanthoorGrowthRespiration(MCbuffruit_tot,MCbufleaf,MCbuf_stemroot,LAIIn,Cbuf);
    % CarboHydrate Buffer 
    [Cbufdot] = BramVanthoorCarboBuffer(MCairbuf, MCbuffruit_tot, MCbufleaf, MCbuf_stemroot, MCbufair_g);
    % Maintenance Respiration
    [MCstemroot_air_m, MCleafair_m, MCfruitair_m] = BramVanthoorMaintenanceRespiration(MCbufleaf,Tcan,Cstemroot,Cleaf,Cfruit);
    % Leaf Harvest
    [MCleafhar] = BramVanthoorLeafHarvest(Cleaf);
    % Plant Buffers
    [Chardot,LAIdot,Cstemrootdot,Cleafdot,Cfruitdot,Nfruitdot] = BramVanthoorPlantBuffer(MCbuffruit_tot,MCbufleaf,MCbuf_stemroot,MCstemroot_air_m, MCleafair_m, MCfruitair_m,MCleafhar,Tcan,Tsum_2,Tcan24_2,Cfruit,Nfruit);
    % Tsum & Tcan24
    Tsumdot    = Tcan/86400;
    Tcan24dot  = 1/86400*(Tcan-Tcan24);
%% Compose State Vector
    StateOut = [Cbuf Cstemroot Cleaf Cfruit Nfruit Char Tsum_1 Tcan24_1 Tsum_2 Tcan24_2];
%% Additional Output Variables
    % Save the value of some of the Mass Flows to be Plotted afterwards. 
    OutputFlows = [
        MCleafhar,MCairbuf,MCbuffruit_tot,MCbufleaf,MCbuf_stemroot,...
        MCstemroot_air_m,MCleafair_m,MCfruit,MCfruitair_m,MCbufair_g];
    LAIOut = LAI;
end

