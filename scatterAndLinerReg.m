%% Multiple Linear Regression by Matlab functions
% MSc Thesis Sensitivity analysis of a Chinese Solar Greenhouse
% Author: Citlali Melchor Ramírez
% December 2021
%% Import data
clear
clc
load GlobalSA_varMethod_v5.mat
load S_temperature.mat
S_temperature(2,:) = [];
S_temperature(5,:) = [];
S_temperature(15,:) = [];
S_temperature(26,:) = [];
S_temperature(71,:) = [];
S_temperature(85,:) = [];
S_temperature(103,:) = [];
ParamNames=readtable('ParamNames.xlsx');
for i=1:50
S_temperature.Properties.RowNames{i} = char(ParamNames.ParamName(i));
end

K = [A; B];
Y = [y_A; y_B];

% for i=1:50
%     C=B;
%     C(:,i)=A(:,i);
%     K=[K;C];
%     Y=[Y; y_C(:,i)];
% end
%% Distribution of SA output
orangeColor = [234 150 81]./255;
figure
set(gcf, 'units','centimeters','Position',  [0, 0, 18.3, 8.9])
set(axes,'YGrid','on')
hold on
histogram(Y,20,'FaceColor',orangeColor, 'EdgeColor', orangeColor)
hold on
mu = mean(Y);
sigma = std(Y);
y = min(Y)-0.5:0.1:max(Y);
f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
plot(y,700*f,'LineWidth',1.5)
ylabel('Frequency, [-]','FontSize', 12)
xlabel('Air Temperature at noon, [°C]','FontSize', 12);
xlim([min(Y)-0.5 max(Y)])
%% Scatter plots
blueColor =[195 215 232]./255;
blueD =[62 117 166]./255;
figure
set(gcf, 'units','centimeters','Position',  [0, 0, 20, 20])
for i=1:25
    P1 = subplot(5,5,i, 'align');
    scatter(K(:,i),Y,'.','MarkerEdgeColor', blueColor ,'MarkerFaceColor', blueColor)
    ref = refline;
    ref.LineWidth = 1; ref.Color=blueD;
    xlim([min(K(:,i)) max(K(:,i))])
    label = char(S_temperature.Properties.RowNames(i));
    set( get(P1,'XLabel'), 'String', label );
    if i == 1
        ylabel('T_{air}, [°C]')
    end
end
%%
figure
set(gcf, 'units','centimeters','Position',  [0, 0, 20, 25])
for i=1:25
    P2 = subplot(5,5,i, 'align');
    scatter(K(:,25+i),Y,'.','MarkerEdgeColor', blueColor ,'MarkerFaceColor', blueColor)
    ref = refline;
    ref.LineWidth = 1; ref.Color=blueD;
    xlim([min(K(:,25+i)) max(K(:,25+i))])
    label = char(S_temperature.Properties.RowNames(25+i));
    set( get(P2,'XLabel'), 'String', label );
    if i ==1
        ylabel('T_{air}, [°C]')
    end
end
%% linear regression

% lm = stepwiselm(K,Y, 'quadratic')
% anova(lm,'summary')
lm_linear = fitlm(K,Y,'linear');

figure
set(gcf, 'units','centimeters','Position',  [0, 0, 18.3, 8.9])
ph = plot(lm_linear);
set(ph(1), 'MarkerEdgeColor',  blueD);
set(ph(2), 'Color', orangeColor);
grid off
xlabel('$Adjusted\ whole\ model$','interpreter','latex','FontSize',12);
ylabel('$Adjusted\ y$','interpreter','latex','FontSize',12);
title('')
% set(get(gca,'XLabel'),'interpreter','latex','FontSize',12);
% set(get(gca,'YLabel'),'interpreter','latex','FontSize',12);


an_linear= anova(lm_linear);

%% Check for multicollinearity
% No scaled parameters
R0 = corrcoef(K); % correlation matrix
V=diag(inv(R0));


%% Standardize linear model
K_st=K-mean(K);

lm_inter_st = stepwiselm(K_st,Y,'Interactions');
an_inter_st = anova(lm_inter_st); 

R0_st = corrcoef(K_st); % correlation matrix
V_st=diag(inv(R0_st));
% for i=
% betas = lm.coefficients.*(std()/std(Y))
%%
% save MLRModels_v1.mat lm_linear lm_inter an_linear