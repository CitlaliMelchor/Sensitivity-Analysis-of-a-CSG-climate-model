function FIRtran = FIRtranBo( label1, firtran, label2, n_blanketi, n_blankete, u_blanket , outlayer, n_canopy, n_floor)
%% FIRtranBo calculate the FIR transmission between state. 
% canopy , inside thermal screen, cover material
% label1
% [1-floor, 2-cover, 5-blanket and sky, 0-others]
% firtran, the FIR transmission of material. [canopy, cover, screen]
% label2, the internal or middel layer  [1-inner face, 0-middle, -1-outer face]
% outlayer   [walle 1, roofe 1, blanket 1]
%% 

n = length(label1);
P = repmat(label1,1,n);
P = P+P';
P = P-diag(diag(P));    % diagonal == 0;
firtran = roundn(firtran,-2);
firtran = firtran + 100;       %exchange incase transmission = 0;

P(P==1)=firtran(1);                       % canopy
P(P==2)=firtran(3);                       % inside screen
P(P==3)=firtran(1)*firtran(3);            % canopy and inside screen
P(P==4)=0;                                % between cover 1 and cover 2
P(P==5)=firtran(2)*firtran(3);            % cover and inside screen
P(P==6)=firtran(1)*firtran(2)*firtran(3); % canopy and cover and inside screen
P(P==7)=0;                                % the transmission between cover and blanket
P(P==0)=1;                                % convert 0 to 1
P(P==10)=1;                               % blanket to sky

P(P==firtran(1))=firtran(1)-100;          % convert back
P(P==firtran(3))=firtran(3)-100;
P(P==firtran(1)*firtran(3))=(firtran(1)-100)*(firtran(3)-100);
P(P==firtran(2)*firtran(3))=(firtran(2)-100)*(firtran(3)-100);
P(P==firtran(1)*firtran(2)*firtran(3))=(firtran(1)-100)*(firtran(2)-100)*(firtran(3)-100);

label2(label2<0)=0;
L = label2*label2';          % the FIR exchange layers give the value of 1 (wal_i,rof_i,floor,cover,canopy,blanket,screen,sky), others is 0.

FIRtran1 = P.*L;

b = zeros(n);                          % change the external layer to 1
b(:,n) = outlayer;
b(n,:) = outlayer';

FIRtran = (1-b).*FIRtran1 + b;

a = FIRtran(n_blanketi,n_blanketi);
FIRtran(n_blanketi,:) = FIRtran(n_blanketi,:).*u_blanket;
FIRtran(:,n_blanketi) = FIRtran(:,n_blanketi).*u_blanket;
FIRtran(n_blanketi,n_blanketi) = a.*u_blanket;

FIRtran(n_canopy, n_floor) = 1;         % the tran between canopy and floor is 1
FIRtran(n_floor, n_canopy) = 1;


d = (1-outlayer).*FIRtran(:,n);         % elements to sky with control u_blanket

FIRtran(:,n) = FIRtran(:,n)-d.*u_blanket;
FIRtran(n,:) = FIRtran(n,:)-d'.*u_blanket;
FIRtran(n_blankete,n) = u_blanket;
FIRtran(n,n_blankete) = u_blanket;

FIRtran(n_blanketi,n) = 0;
FIRtran(n,n_blanketi) = 0;

end

