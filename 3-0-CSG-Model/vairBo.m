function V_air = vairBo( Pspg )
%the volum of inside air [m3 m-1]

%Geometry of CSG
% AB = Pcsg(1,8)*Pcsg(1,1)/cos(Pcsg(2,8));
BC = Pspg(2,1)-Pspg(1,8)*Pspg(1,1)*tan(Pspg(2,8));
% CD = Pcsg(1,1);
% AD = ((Pcsg(2,1))^2+(Pcsg(1,1)*(1-Pcsg(1,8)))^2)^0.5;
% AC = ((Pcsg(2,1))^2+(Pcsg(1,1)*Pcsg(1,8))^2)^0.5;
% BD = (BC^2+CD^2)^0.5;

Leng_CSG=1;                        %The length of CSG is 1 m.

a6=Pspg.*[0,1,1,1,1,1,1,0;0,1,1,1,1,1,1,0];
a3=Pspg(1)*a6;
a4=Pspg(2)*a6;
a4=a4';
a5=0.5*a3*a4;
Aair=a5(3)+ Pspg(1)*Pspg(1,3)*Pspg(2)*Pspg(2,2) ...
    + Pspg(1)*Pspg(1,4)*Pspg(2)*(Pspg(2,2)+Pspg(2,3))...
    +Pspg(1)*Pspg(1,5)*Pspg(2)*(Pspg(2,2)+Pspg(2,3)...
    +Pspg(2,4))+Pspg(1)*Pspg(1,6)*Pspg(2)*(Pspg(2,2)...
    +Pspg(2,3)+Pspg(2,4)+Pspg(2,5))+Pspg(1)*Pspg(1,7)...
    *Pspg(2)*(Pspg(2,6)+Pspg(2,2)+Pspg(2,3)+Pspg(2,4)...
    +Pspg(2,5))+Pspg(1)*Pspg(1,8)*(Pspg(2)*(Pspg(2,7)...
    +Pspg(2,2)+Pspg(2,3)+Pspg(2,4)+Pspg(2,5)+Pspg(2,6)+BC)/2);

V_air = Aair*Leng_CSG;       %Air volum    [m3 m-1]

end

