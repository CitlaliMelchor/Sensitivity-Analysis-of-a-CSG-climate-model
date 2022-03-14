function [angle1,length1] = angleBo( c )
%the angle of each panel to south (x); the length of the panel
%   c is the coordinate of the panel start and end points
%   [x1st, y1st, x1end, y1end; 
%    x2st, y2st, x2end, y2end;
%    ...]  with Clockwise direction

length1 = ((c(:,1)-c(:,3)).^2+(c(:,2)-c(:,4)).^2).^0.5;

dx = c(:,1)-c(:,3);
dy = c(:,2)-c(:,4);

angle1 = atan(dy./dx)+pi.*(1-switch01(dy.*dx,-1));

end

