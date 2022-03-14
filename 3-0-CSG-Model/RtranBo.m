function Rtran = RtranBo( label1, radtran, label2 , n_blanketi, n_blankete, u_blanket, outlayer)
%% FIRtranBo calculate the radiation(FIR, NIR, PAR) transmission between state. 
% canopy , inside thermal screen, cover material
% label1
% [1-floor, 2-cover, 5-blanket and sky, 0-others]
% radtran, the Rad transmission of material. [canopy; cover; screen]
% label2, the internal or middel layer  [1-inner face, 0-middle, -1-outer face]
% cover and blanket have only 1 layer.
%% 

n = length(label1);
P = repmat(label1,1,n);
P = P+P';
P = P-diag(diag(P));    % diagonal == 0;
radtran = roundn(radtran,-2);
radtran = radtran + 100;       %exchange incase transmission = 0;

P(P==1)=radtran(1);                       % canopy
P(P==2)=radtran(3);                       % inside screen
P(P==3)=radtran(1)*radtran(3);            % canopy and inside screen
P(P==4)=0;                                % between cover 1 and cover 2
P(P==5)=radtran(2)*radtran(3);            % cover and inside screen
P(P==6)=radtran(1)*radtran(2)*radtran(3); % canopy and cover and inside screen
P(P==7)=0;                                % the transmission between cover and blanket
P(P==0)=1;                                % convert 0 to 1
P(P==10)=1;                               % blanket to sky

P(P==radtran(1))=radtran(1)-100;          % convert back
P(P==radtran(3))=radtran(3)-100;
P(P==radtran(1)*radtran(3))=(radtran(1)-100)*(radtran(3)-100);
P(P==radtran(2)*radtran(3))=(radtran(2)-100)*(radtran(3)-100);
P(P==radtran(1)*radtran(2)*radtran(3))=(radtran(1)-100)*(radtran(2)-100)*(radtran(3)-100);

label2(label2<0)=0;
L = label2*label2';          % the Rad exchange layers give the value of 1 (wal_i,rof_i,floor,cover,canopy,blanket,screen,sky), others is 0.

Rtran1 = P.*L;

b = zeros(n);                          % change the external layer to 1
b(:,n) = outlayer;
b(n,:) = outlayer';

Rtran = (1-b).*Rtran1 + b;

a = Rtran(n_blanketi,n_blanketi);
Rtran(n_blanketi,:) = Rtran(n_blanketi,:).*u_blanket;
Rtran(:,n_blanketi) = Rtran(:,n_blanketi).*u_blanket;
Rtran(n_blanketi,n_blanketi) = a.*u_blanket;

Rtran(n_blankete,n) = u_blanket;
Rtran(n,n_blankete) = u_blanket;

Rtran(n_blanketi,n) = 0;
Rtran(n,n_blanketi) = 0;

end





