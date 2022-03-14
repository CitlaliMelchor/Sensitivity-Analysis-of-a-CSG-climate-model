% function to describe temperature filter based upon Avissar 1985
% function [h_c h_c_der] = inhibition_continue(X, S1,b1,S2,b2)

function [h_c, h_c_der] = inhibition_continue(X, S1,b1,S2,b2)

% the first flank of the curve
h_c1 = 1/ (1 + exp(-S1*(X - b1)));

% the second and last flank of the curve
h_c2 = 1/ (1 + exp(-S2*(X - b2)));

% the total inhibition
h_c = h_c1*h_c2;

% ----------------------------The derivitave of h_c --> X
h_c_der = (1+exp(-S1*(X-b1))).^-2.*S1.*exp(-S1*(X-b1)).* 1./(1 + exp(-S2*(X-b2)))...
         + 1./(1 + exp(-S1*(X-b1))).*(1+exp(-S2*(X-b2))).^-2.*S2.*exp(-S2*(X-b2));
     

     
     
 
