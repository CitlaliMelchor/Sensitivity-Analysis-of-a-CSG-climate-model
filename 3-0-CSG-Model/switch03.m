function sw = switch03(x,q)
%switch03: switch the results between 0 and 1. smooth.
%while q=-1, if x>0, ans=1; else x<=0 ans=0.
%while q=1,  if x>=0, ans=0; else x<0 ans=1.
k  = 10;
a  = 1./(1+exp(q*k*x));
sw = a;

end

