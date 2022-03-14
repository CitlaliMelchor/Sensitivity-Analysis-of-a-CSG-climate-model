function VPsat = vapourPsat( T )
%vapourPsat: calculate saturate vapour pressure based on Magnus equations
%   T is in oC
VPsat = 10.^(2.7857+9.5*T./(265.5+T)).*switch01(T,1)+10.^(2.7857+7.5*T./(237.3+T)).*switch01(T,-1);


end

