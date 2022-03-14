function view_dif = viewdifBo( angle1,label2 )
%calculate the view facotr of greenhouse panel to diffuse radiation
%  angel1 -- the angle of each panel
%  length1 -- the length of each panel
%  label2 -- the internal or middel layer  [1-inner face, 0-middle, -1-outer face]

a = angle1.*label2;     % outer face positive
a(a<0) = a+pi;          % add pi to outer face

view_dif = a./pi;

% floor ?



end

