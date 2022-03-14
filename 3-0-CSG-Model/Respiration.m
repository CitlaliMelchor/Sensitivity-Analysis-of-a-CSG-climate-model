
% function to calculate maintenance and crop respiration
%
% function [MCcrpair_m, MCcrpair_g] = Respiration(Ccrp, Tcan, MCaircrp, rm_crp, Q10m_crp, c_g)
% MCcrpair_m, MCcrpair_g in kg[C02].m-2.s-1


function [MCcrpair_m, MCcrpair_g, MCcrpair_m_der] = Respiration(Ccrp, Tcan, MCaircrp, rm25, Q10m, c_g)

Tref = 25;
% Maintenance respiration of the crop
MCcrpair_m = rm25 * Q10m^(0.1*(Tcan - Tref))*Ccrp;                          % original
%Ilias Tsafaras 16/12/2013
%MCcrpair_m = min(rm25 * Q10m^(0.1*(Tcan - Tref))*Ccrp,MCaircrp);             % mg[C02].m-2.s-1,
%MCcrpair_m = min(rm25 * Q10m^(0.1*(Tcan - Tref))*Ccrp,MCbufleaf);             % mg[C02].m-2.s-1,

% For sensitivty equations
% MCcrpair_m_der =    Ccrp * rm25 * Q10m.^(0.1* -Tref).* 0.1.*Q10m.^(0.1*Tcan)*log(Q10m);
MCcrpair_m_der =0;

% Growth respiration replaced by smoothing function
% if MCaircrp > MCcrpair_m
%     MCcrpair_g = c_g*(MCaircrp - MCcrpair_m);                                           % kg[C02].m-2.s-1
% else
%     MCcrpair_g =0;
% end

% Smoothing function 
Slope = 10e7;
ValueIfElse = SmoothIfElse(MCaircrp, MCcrpair_m, Slope);
MCcrpair_g = c_g*(MCaircrp - MCcrpair_m)*(ValueIfElse);