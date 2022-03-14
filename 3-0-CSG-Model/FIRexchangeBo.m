function FIR = FIRexchangeBo( vf, tem,n_screen,n_canopy,n_cover,n_blanketi,n_blankete,n_floor,u_blanket,spg )
%% FIRexchangeBo calculate the FIR heat exchange between the elements
%   vf is the area * view factor               [n*n]
%   emis is the list of FIR emission factor    [n]
%   tem is the list of temperature             [n]
%   FIRtran is the transmission of FIR related to the control vector  [n*n]
%   (inside and outside screen)
%% 
emis = spg.FIRabs;
  firtran = [spg.FIRtran(n_canopy,1),spg.FIRtran(n_cover,1),spg.FIRtran(n_screen,1)];
  FIRtran = FIRtranBo( spg.label1, firtran, spg.label2, n_blanketi, n_blankete, u_blanket,spg.outlayer,n_canopy, n_floor );

global Sigma

n = length(tem);          % the number of elements
tem1 = (tem+273.15).^4;   % oC 2 K
tem2 = repmat(tem1,1,n);  % repate the list to matrix
tem3 = tem2' - tem2;      % temperature difference

emis1 = emis * emis';     % emission array to matrix

FIR1 = Sigma.*vf.*emis1.*FIRtran.*tem3;

%%% if there is no screen inside, thermal screen = 0 %%%
FIR1(n_screen,:) = 0;
FIR1(:,n_screen) = 0;


% FIR = diag(FIR1);
FIR = -sum(FIR1)';

end

