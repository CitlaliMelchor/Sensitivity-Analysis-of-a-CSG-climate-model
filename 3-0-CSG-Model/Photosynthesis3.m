
% Photosynthesis model based upon Faraquhar (1980 & 1982) reconstructed out of Intkam model 
% by Bram Vanthoor October 2008
% 
% Explanation of names Intkam <-> Oothegem
% Inputs
% Tcan:     oC canopy temperature 
% PARabs:   W/m2 Absorbed PAR transmission
% CO2air:   ppm
% Pg_L:     micromol {CO2}.m-2.s-1  canopy photosynthesis

function [Pg_L,VC2,FRESP,J] = Photosynthesis3(Tcan, PARABS, CO2air,JMAX25,THETA,LAI)

% Tcan = 25;
% PARABS = 200;                 % W/m2 Absorbed PAR transmission
% CO2air = 350;                 % ppm, micromol(CO2)/mol (air)          

% Photosynthesis Parameters

% -----Copied from Intkam Intput file: HortiTomaat.par------%
% Values don't agree with Oothegem --> Find literature reference

% JMAX25 = 386.78;%;%386.78;%250;      % Maximal rate of electron transport [mumol electrons m-2 s-1] at 25 oC
VCMAX25 = 0.5*JMAX25;%200;%175;%125;     % Maximal rate of carboxylation [mumol CO2 m-2 s-1] at 25 oC
RD25 = 1.25;       % Rate of dark respiration of leaf at 25 oC [mumol CO2 m-2 s-1]

% Values agreed with Oothegem
KC25 = 310;         % M.M. constant of Rubisco for binding to CO2 at 25 oC (mumol mol-1) == microbar
KO25 = 155;         % M.M. constant of Rubisco for binding to O2 at 25 oC (mmol mol-1) = milibar
p_O2can = 210;      % O2 partial pressure stomata in milibar
LGHTCON = 4.6;      % Conversion J PAR to mumol photons
R = 8.314;          % Gas constant J.mol-1.K-1
F =0.23;            % F is fraction of photons absorbed by non-photosynthetic tissues
S = 710;            % J mol-1K-1
H = 220000;         % J mol-1
eta_CO2air_CO2can = 0.67;   % Coversion from co2 concentration in air to co2 inside stomato 
EC = 59356;         % J mol-1 Activation energy KC rubisco carboxylation 
EO = 35948;         % J mol-1 Activation energy KO rubisco oxygenation
EJ = 37000;         % J mol-1 Activation energy Jmax, maximum electron transport rate
EVC = 58520;        % J mol-1 Activation energy VCmax, maximum carboxylation rate
MCO2 = 44e-3;       % CO2 Molair mass, kg{CO2}  .mol-1{CO2}

% Conversions and assumptions 
% The intern CO2-concentration is not calculated but assumed to be a fixed 
% fraction of the outdoor concentration--> could be
% implemented....--> becomes a state of the model
CO2can = eta_CO2air_CO2can*CO2air;       % ppm, Internal CO2-concentration!!!!!
% CO2can = CO2can/1.8040;     % ppm == micromol(co2)/micromol(air) Internal CO2-concentration
Tcan = Tcan + 273.15;       % Coversion form °C --> Kelvin
T25 = 25 + 273.15;          % Kelvin Temperature at 25°C

% Michaelis Menten Curves
X = (Tcan - T25)/(Tcan*R*T25);

KC = KC25 * exp(EC*X);
KO = KO25 * exp(EO*X); 
VCMAX = VCMAX25 *exp(EVC*X); 

% CO2 compensation point increases with temperature according to Brooks & Farquhar, 1985
% GAMMA1 = (42.7 + 1.68 * (Tcan - T25) + 0.012 * (Tcan - T25)^2);   % Differs with van Ootegem 
% GAMMA1 = 34.6;%     % 34.6 at 20°C,  42.7 at 25°C
% GAMMA1 = 34.5+0.74*(Tcan -T25+5);

b = 1.7*20*(1-1/LAI);       % Gamma value at T = 0c
% Introduced september 24 2010 to avoid simulations difficulties
% LAI = max(LAI,0.2);
GAMMA1 = 1.7/LAI*(Tcan -273.15) + b;      % Because Jmax depends on LAI, also GAMMA depends on LAI
% GAMMA1 = 1.7*(Tcan-273.15 );

if GAMMA1 < 5
    GAMMA1 = 5;
end

% To calculata Jmax
% TJMAX according to pers. comm. Farquhar to Ad Schapendonk (1985)
D1 = 1. + exp( (S-H/T25) / R );
D2 = 1. + exp( (S-H/Tcan) / R );
TJMAX = exp( EJ * X ) * D1 / D2;

% Temperature dependent potential rate of electron transport mu Eq m-2 s-1 (16.33) + (16.34): F,1982 see also F,1982 fig 16.7
% TJMAX =1;
JMAX = JMAX25 * TJMAX ;     


% Calculate potential rate of electron transport (mumol e- m-2 s-1) 2 electrons per absorbed photon
% alfa = (1. - F)/2;
alfa =0.385;           % F = 0.23 according farquhar--> means alfa of 0.385

EFFRAD = PARABS * LGHTCON*alfa ;       % Conversion of J m-2 s-1 to mumol m-2 s-1 with LGHTCON


J = ( JMAX + EFFRAD -  sqrt( (JMAX+EFFRAD)^2-4.*THETA*EFFRAD*JMAX )  )/ (2. * THETA);   % potential rate of electron transport (mumol e- m-2 s-1)   

%------------------Carboxylationrate-------------------------------------%

% Calculate RuP2-saturated rate of carboxylation (16.59) F, 1982 (is limiting rate of carboxylation)
% Because of canopy photosynthesis this one is not used!

%  VC1 = VCMAX * (CO2can )/ ( CO2can + KC * (1.+p_O2can/KO) );

% Calculate electron transport/photophosphorylation limited rate of RuP2 regeneration
% Division by (4.5*CO2can+10.5*GAMMA)  assumes pseudocyclic electron transport.
%	!ae,18/1/2000: deze levert iets hoger productie

JC = J; 

%  VC2 = JC * (CO2can-GAMMA1) / ( 4.5*CO2can + 10.5*GAMMA1 );  %    Used by van Oothegem

VC2 = JC / 4. * (CO2can- GAMMA1) / ( CO2can + 2. * GAMMA1 );  % Used by Intkam.
if VC2 > 1.4e5
%     p=1; 
end
% % Compute actual carboxylation velocity
% VC = min( VC1, VC2 );

% Rubisco is not limiting
VC=VC2;
%---------End of Carboxylationrate-------------------------------------%

% Compute photorespiration (16.3)+(16.18) F,1982
FRESP = VC * GAMMA1/ CO2can;

% Leaf gross photosynthesis (mumol CO2 m-2 s-1) (16.57) F,1982
FARPHG = VC - FRESP;

Pg_L = FARPHG;















