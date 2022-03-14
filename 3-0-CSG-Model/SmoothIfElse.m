% function to replace if then else statement with a motthed version
% based upon PhD of Van Oothegem 2007 page 10, by Bram Vanthoor October
% 2007
% 
% function [SmoothValue] = SmoothIfElse(Value, SwitchValue, Slope)
% Value:            Value of the state      e.g. 10 m/s (actual windspeed)
% Switch Value:     Value where to switch   e.g.  4 m/s (switch)  
% Slope:            Slope of the function   e.g.  10 BUT watch out depends of the absolute values of the switch value !
% Smoothvalue:      Output of the function, lies between [0.....1]

function [SmoothValue] = SmoothIfElse(Value, SwitchValue, Slope)

% SmoothValue = 1/(1 + 10^(-Slope* (Value - SwitchValue)));

SmoothValue = 1./(1 + exp(Slope.*(Value - SwitchValue)));   % Changed at 8 december 2008
