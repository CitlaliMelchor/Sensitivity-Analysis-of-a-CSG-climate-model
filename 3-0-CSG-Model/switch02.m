function sw = switch02(x,q)
%switch02: switch the results between 0 and 1.
%while q=-1, if x>=0, ans=1; else x<0 ans=0.
%while q=1,  if x>0, ans=0; else x<=0 ans=1.
k  = 1e16;
a  = 1./(1+exp(q*k*x));
sw = round(a);

end