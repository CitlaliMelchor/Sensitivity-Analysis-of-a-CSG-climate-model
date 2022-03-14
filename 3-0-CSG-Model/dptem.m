function dp = dptem( tem, rh )
%dp dew pont temperature.
%   dp = -60.45+7.0322*log(e)+0.3700*log(e)^2         -60<dp<=0
%   dp = -35.957-1.8726*log(e)+1.1689*log(e)^2        0<dp<=70
esat = vapourPsat(tem);
e    = esat.*rh;

dp1 = -60.45+7.0322*log(e)+0.3700*log(e).^2;
dp2 = -35.957-1.8726*log(e)+1.1689*log(e).^2;
sw1 = abs(switch02(dp1,1)-switch01(dp2,-1)); % incase dp1<=0 and dp2>0

dp = (dp1.*switch02(dp1,1)+dp2.*switch01(dp2,-1)).*sw1;

dp = (dp1+dp2)./2.*switch02(sw1,1)+dp.*switch01(sw1,-1); %if sw1=0, dp = (dp1+dp2)/2

end

