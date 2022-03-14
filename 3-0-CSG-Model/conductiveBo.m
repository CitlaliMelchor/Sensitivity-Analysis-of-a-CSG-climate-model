function sens_cond = conductiveBo( tem, thickness, area, heatcond, labelcond, T_soilbound ,n_cover, n_blanketi, n_blankete, u_blanket,I )
%Sensible heat exchange by heat conductive. 
%

labelcond(n_blanketi:n_blankete)= labelcond(n_cover)+(1-u_blanket);


% a = length(tem);
b = zerodivided0(0.5*thickness,heatcond);
% db = b((1:(a-1)),:)+b((2:a),:);
db = movsum(b,2);
HEC = zerodivided0(1,db(2:end));                   %[W / m2 K]
% dT = -tem((1:(a-1)),:)+tem((2:a),:);
dT = diff(tem);
dH = HEC.*dT;

% I = find(labelcond==1,1,'last');           % find the last layer of soil  --- 1

dH_soilb = heatcond(I)/(0.5*thickness(I))*(-tem(I)+T_soilbound);     % the heat exchange between the last soil layer and boundry

%c = labelcond((1:(a-1)),:)-labelcond((2:a),:);    % same label indicate conductive inside
c = diff(labelcond);
% c(c==0)=-100;
% c(c>-100)=0;
% c(c==-100)=1;                                     % conductive c =1 ; none conductive c = 0;
c = ~logical(c);
dH1 = dH.*c;
dH2 = -[0;dH1]+[dH1;0];                            %

dH2(I)=dH2(I)+dH_soilb;


sens_cond = dH2.*area;  %[W per length greenhouse]

end

