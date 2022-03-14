function div = zerodivided0(a,b)
%zerodivided0: a divided by b (element by element).
%   if the value of b is 0, the answer is 0; else ans=a/b. 
c = switch01(abs(b),-1);
d = switch02(abs(b),1);

div = a./(b+d).*c;

end

