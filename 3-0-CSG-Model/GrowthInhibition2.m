% Yield reduction by temperature influences
% Calculated for 2 approaches:
%   - effect of instantaneous temperature
%   - effect of daily temperature sum

function [h_Tcan_inst,h_Tcan_24,h_Tcan_inst_der,h_Tcan_24_der,MCaircrp_mean,h_DIF,h_Tcan_photo,h_Tcan_photo_der] = ...
         GrowthInhibition2(Tcan,Tcan24,T24_S1,T24_b1,T24_S2,T24_b2,Tinst_S1,Tinst_b1,Tinst_S2,Tinst_b2) % ,Tmax,Tmin)
      
% effect of instantaneous temperature % based on Boote(2006) although it is an average for the last 24 hours
% the last 24 hours are chosen, because night temperatures also influences
% daily growth!

% Calculate the instantaneous influence of Tcan upon crop growth: h_cropyield_T
[h_Tcan_24, h_Tcan_24_der]= inhibition_continue(Tcan24, T24_S1,T24_b1,T24_S2,T24_b2);

% Determine inhibition on photosynthesis by low 24 hour mean temperature, only inhibition at low temperatures
% [h_Tcan_photo, h_Tcan_photo_der]= inhibition_continue(Tcan24, T24_S1,T24_b1,T24_S2,200);
[h_Tcan_photo, h_Tcan_photo_der]= inhibition_continue(Tcan24, T24_S1,T24_b1,T24_S2,T24_b2);
% [h_Tcan_inst ]= inhibition_continue(Tcan, 2, 12, 28, 45);
[h_Tcan_inst, h_Tcan_inst_der]= inhibition_continue(Tcan, Tinst_S1,Tinst_b1,Tinst_S2,Tinst_b2);

% T1 = Tcan; T2 = Tcan24;
% S1 = Tinst_S1; S2 = Tinst_S2; b1 = Tinst_b1; b2 = Tinst_b2;
% S3 = T24_S1; S4 = T24_S2; b3 = T24_b1;b4 = T24_b2;

% ----------------For the first element-----------------------------------
% h1 = 1./(1 + exp(-S1*(T1-b1))) * 1./(1 + exp(-S2*(T1-b2)));
% 
% h2 = 1./(1 + exp(-S3*(T2-b3))) * 1./(1 + exp(-S4*(T2-b4)));

% h_diff = (1+exp(-S1*(T-b1))).^-2.*S1.*exp(-S1*(T-b1)).* 1./(1 + exp(-S2*(T-b2)))...
%          + 1./(1 + exp(-S1*(T-b1))).*(1+exp(-S2*(T-b2))).^-2.*S2.*exp(-S2*(T-b2));

% -----------------------For both effects, Tcan + Tcan24------------------
% h_diff1 = (1+exp(-S1*(T1-b1))).^-2.*S1.*exp(-S1*(T1-b1)).* 1./(1 + exp(-S2*(T1-b2)))...
%          + 1./(1 + exp(-S1*(T1-b1))).*(1+exp(-S2*(T1-b2))).^-2.*S2.*exp(-S2*(T1-b2));
%      
% h_diff2 = ((1+exp(-S3*(T2-b3))).^-2.*S3.*exp(-S3*(T2-b3)).* 1./(1 + exp(-S4*(T2-b4)))...
%          + 1./(1 + exp(-S3*(T2-b3))).*(1+exp(-S4*(T2-b4))).^-2.*S4.*exp(-S4*(T2-b4)));
     
    
% h_diff_tot = h1.*h_diff2*Tcan24_Tcan_der + h2.*h_diff1;

% Photosynthesis rate based upon mean momentanuos inhibition
%MCaircrp_mean_inst = h_Tcan_inst*h_Tcan_24*MCaircrp_pot;

% Only mean inhibition
% MCaircrp_mean = h_Tcan_24*MCaircrp_pot;
MCaircrp_mean = 0;

% Mean inhibition and DIF impact 
% Tdif = Tmax - Tmin;

% [h_DIF] = inhibition_continue(Tdif, TDIF_S1,TDIF_b1,TDIF_S2,TDIF_b2);
h_DIF = 0;
% MCaircrp= h_Tcan_24*h_DIF*MCaircrp_pot;
%MCaircrp= MCaircrp_pot;

