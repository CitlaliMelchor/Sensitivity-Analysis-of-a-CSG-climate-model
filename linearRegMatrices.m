%% Multiple Linear Regression by matrix algebra
% MSc Thesis Sensitivity analysis of a Chinese Solar Greenhouse
% Author: Citlali Melchor Ramírez
% December 2021

clear; clc
load GlobalSA_varMethod_v5.mat A B y_A y_B y_C
Y=[y_A;y_B];
K =[A;B];

M = [ones(600,1), K]; % Parameters MonteCarlo
% Coefficients vector
B = (M'*M)\M'*Y;
% y estimate
% epsilon = Y - M*B; 
Y_hat = M*B;

SS_tot = (Y - mean(Y))'*(Y - mean(Y));          % Total sum of squares
SS_reg = (Y_hat - mean(Y))'*(Y_hat - mean(Y));  % Regresion sum of squares
SS_res = (Y_hat - Y)'*(Y_hat - Y);              % Residual sum of squares

R_sq = SS_reg/SS_tot;

lm = fitlm(K,Y,'RobustOpt','off');
%%
SY = std(Y);
for i = 1:24
    SK(i) = std(K(:,i));
    BETA(i) = B(i+1)*(SK(i)/SY);
end

sum(BETA)
plotResiduals(lm)
% Orthogonality of regressors
%https://stats.stackexchange.com/questions/223818/can-standardized-beta-coefficients-in-linear-regression-be-used-to-estimate-t#comment424193_223818
%% Check for correlations and collinearity
% Correlation matrix
R0 = corrcoef(K);

VIF = diag(inv(R0))';

R_sq = lm.Rsquared;

R_sq_i = 1-(1./VIF);

% Get the eigenvalues
 [V,D]=eig(R0); 
 
 %variance proportions
 %lambda = eigenvalue D
 %w~ eigenvector V pairs for M = X'X
 c= sum(V.^2./diag(D));
 for i =1:50
     for j = 1:50 
         p(i,j) = V(i,j)^2 /(D(i,i) *c(j));
     end
 end