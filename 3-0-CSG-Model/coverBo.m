function Shape_cover = coverBo( Pcsg )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


a1     = (Pcsg(1,1)*Pcsg).^2;
a2     = (Pcsg(2,1)*Pcsg).^2;

slop_srof1=atan(Pcsg(2)*Pcsg(2,:)./(Pcsg(1)*Pcsg(1,:)));                % calculate the slope of south pane [rad]
slop_srof = slop_srof1(2:7);                                            % slope of the south roof pane  [rad]
Leng1 = (a1(1,:)+a2(2,:)).^0.5;                                         % calculate the length of each pane [m]
Leng = Leng1(2:7);                                                      % the length of each pane [m]

Shape_cover = [slop_srof',Leng'];

end

