%% Pretty plots
% MSc Thesis Sensitivity analysis of a Chinese Solar Greenhouse
% Author: Citlali Melchor Ramírez
% December 2021

%% Local Sensitivity Analysis

% Nice green and red pallette for the sensitivities
negativeColor = ([[242 191 150];[234 150  81];[210 109 25];[140 72 16]])./255;
positiveColor = ([[209 218 180];[179 194 131];[143 163 80];[95 108 53]])./255;

% Temperature
load S_temperature.mat
% S_temperature.Properties.RowNames{21} = 'T_{min\_night}';
% S_temperature.Properties.RowNames{32} = 'T_{min\_day}';
% S_temperature.Properties.RowNames{8} = 'sky&Air_{FIRabs}';
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
%%
figure
L1 = plot([NaN],[NaN],'s', 'MarkerFaceColor', negativeColor(1,:),...
                           'MarkerEdgeColor', negativeColor(1,:));hold on
L2 = plot([NaN],[NaN],'s', 'MarkerFaceColor', negativeColor(2,:),...
                           'MarkerEdgeColor', negativeColor(2,:));
L3 = plot([NaN],[NaN],'s', 'MarkerFaceColor', negativeColor(3,:),...
                           'MarkerEdgeColor', negativeColor(3,:));
L4 = plot([NaN],[NaN],'s', 'MarkerFaceColor', negativeColor(4,:),...
                           'MarkerEdgeColor', negativeColor(4,:));
L5 = plot([NaN],[NaN],'s', 'MarkerFaceColor', positiveColor(1,:),...
                           'MarkerEdgeColor', positiveColor(1,:));
L6 = plot([NaN],[NaN],'s', 'MarkerFaceColor', positiveColor(2,:),...
                           'MarkerEdgeColor', positiveColor(2,:));
L7 = plot([NaN],[NaN],'s', 'MarkerFaceColor', positiveColor(3,:),...
                           'MarkerEdgeColor', positiveColor(3,:));
L8 = plot([NaN],[NaN],'s', 'MarkerFaceColor', positiveColor(4,:),...
                           'MarkerEdgeColor', positiveColor(4,:));
set(gcf, 'units','centimeters','Position',  [0, 0, 18.3, 8.9])
b = bar(table2array(S_temperature(1:20,:)),'FaceColor','flat','EdgeColor','none');
for i=1:20
    for ii=1:4
        if table2array(S_temperature(i,ii))<0
            newColor(ii,:) = negativeColor(ii,:);
        else
            newColor(ii,:) = positiveColor(ii,:);
        end
    end
    for k=1:4
        b(k).CData(i,:) = newColor(k,:);
    end
end
%title('Local SA: Air Temperature')
set(groot,'defaultAxesTickLabelInterpreter','default');  
xticks([1:1:290])
xticklabels(S_temperature.Properties.RowNames)
% a = get(gca,'XTickLabel');
% set(gca,'fontsize',10)
xtickangle(45)
ylim([-0.3 0.3])
ylabel('$\bar{S},[-]$','Interpreter','Latex','FontSize', 12)
box on
L = legend([L1,L2,L3,L4,L5,L6,L7,L8],...
    {' 0.8',' 0.9',' 1.1',' 1.2','','','',''},...
    'NumColumns',2, 'Location','eastoutside');
L.ItemTokenSize(1) = 2;
legend boxoff
title(L,['\theta     ',char(160)]);
% exportgraphics(gcf,'SA_Temperature.png','Resolution',600)
%%
%Relative Humidity
load S_rh.mat
S_rh.Properties.RowNames{3} = 'T_{min\_night}';
S_rh.Properties.RowNames{9} = 'T_{min\_day}';
S_rh.Properties.RowNames{16} = 'sky&Air_{FIRabs}';
figure
L1 = plot([NaN],[NaN],'s', 'MarkerFaceColor', negativeColor(1,:),...
                           'MarkerEdgeColor', negativeColor(1,:));hold on
L2 = plot([NaN],[NaN],'s', 'MarkerFaceColor', negativeColor(2,:),...
                           'MarkerEdgeColor', negativeColor(2,:));
L3 = plot([NaN],[NaN],'s', 'MarkerFaceColor', negativeColor(3,:),...
                           'MarkerEdgeColor', negativeColor(3,:));
L4 = plot([NaN],[NaN],'s', 'MarkerFaceColor', negativeColor(4,:),...
                           'MarkerEdgeColor', negativeColor(4,:));
L5 = plot([NaN],[NaN],'s', 'MarkerFaceColor', positiveColor(1,:),...
                           'MarkerEdgeColor', positiveColor(1,:));
L6 = plot([NaN],[NaN],'s', 'MarkerFaceColor', positiveColor(2,:),...
                           'MarkerEdgeColor', positiveColor(2,:));
L7 = plot([NaN],[NaN],'s', 'MarkerFaceColor', positiveColor(3,:),...
                           'MarkerEdgeColor', positiveColor(3,:));
L8 = plot([NaN],[NaN],'s', 'MarkerFaceColor', positiveColor(4,:),...
                           'MarkerEdgeColor', positiveColor(4,:));
set(gcf, 'units','centimeters','Position',  [0, 0, 18.3, 8.9])
b = bar(table2array(S_rh(1:20,:)),'FaceColor','flat','EdgeColor','none');
for i=1:20
    for ii=1:4
        if table2array(S_rh(i,ii))<0
            newColor(ii,:) = negativeColor(ii,:);
        else
            newColor(ii,:) = positiveColor(ii,:);
        end
    end
    for k=1:4
        b(k).CData(i,:) = newColor(k,:);
    end
end
%title('Local SA: Relative Humidity')
xticks([1:1:290])
xticklabels(S_rh.Properties.RowNames)
xtickangle(45)
ylim([-0.3 0.3])
ylabel('$\bar{S},[-]$','Interpreter','Latex','FontSize', 12)
box on
L = legend([L1,L2,L3,L4,L5,L6,L7,L8],...
    {' 0.8',' 0.9',' 1.1',' 1.2','','','',''},...
    'NumColumns',2, 'Location','eastoutside');
L.ItemTokenSize(1) = 2;
legend boxoff
title(L,['\theta     ',char(160)]);
% exportgraphics(gcf,'SA_RH.png','Resolution',600)

%Soil temperature
load S_soilT.mat
S_soilT.Properties.RowNames{29} = 'T_{min\_night}';
S_soilT.Properties.RowNames{38} = 'T_{min\_day}';
S_soilT.Properties.RowNames{7} = 'sky&Air_{FIRabs}';
figure
L1 = plot([NaN],[NaN],'s', 'MarkerFaceColor', negativeColor(1,:),...
                           'MarkerEdgeColor', negativeColor(1,:));hold on
L2 = plot([NaN],[NaN],'s', 'MarkerFaceColor', negativeColor(2,:),...
                           'MarkerEdgeColor', negativeColor(2,:));
L3 = plot([NaN],[NaN],'s', 'MarkerFaceColor', negativeColor(3,:),...
                           'MarkerEdgeColor', negativeColor(3,:));
L4 = plot([NaN],[NaN],'s', 'MarkerFaceColor', negativeColor(4,:),...
                           'MarkerEdgeColor', negativeColor(4,:));
L5 = plot([NaN],[NaN],'s', 'MarkerFaceColor', positiveColor(1,:),...
                           'MarkerEdgeColor', positiveColor(1,:));
L6 = plot([NaN],[NaN],'s', 'MarkerFaceColor', positiveColor(2,:),...
                           'MarkerEdgeColor', positiveColor(2,:));
L7 = plot([NaN],[NaN],'s', 'MarkerFaceColor', positiveColor(3,:),...
                           'MarkerEdgeColor', positiveColor(3,:));
L8 = plot([NaN],[NaN],'s', 'MarkerFaceColor', positiveColor(4,:),...
                           'MarkerEdgeColor', positiveColor(4,:));
set(gcf, 'units','centimeters','Position',  [0, 0, 18.3, 8.9])
b = bar(table2array(S_soilT(1:20,:)),'FaceColor','flat','EdgeColor','none');
for i=1:20 
    for ii=1:4
        if table2array(S_soilT(i,ii))<0
            newColor(ii,:) = negativeColor(ii,:);
        else
            newColor(ii,:) = positiveColor(ii,:);
        end
    end
    for k=1:4
        b(k).CData(i,:) = newColor(k,:);
    end
end
%title('Local SA: Soil Temperature')
xticks([1:1:290])
xticklabels(S_soilT.Properties.RowNames)
xtickangle(45)
ylim([-0.3 0.3])
ylabel('$\bar{S},[-]$','Interpreter','Latex','FontSize', 12)
box on
L = legend([L1,L2,L3,L4,L5,L6,L7,L8],...
    {' 0.8',' 0.9',' 1.1',' 1.2','','','',''},...
    'NumColumns',2, 'Location','eastoutside');
L.ItemTokenSize(1) = 2;
legend boxoff
title(L,['\theta     ',char(160)]);
% exportgraphics(gcf,'SA_soilT.png','Resolution',600)


%% Global Sensitivity Analysis
load GlobalSA_varMethod_v5.mat%GlobalSAv4.mat
valuesMatrix = [A; B];
y2 = [y_A; y_B];
% Distribution of SA output
orangeColor = [234 150 81]./255;
figure
set(gcf, 'units','centimeters','Position',  [0, 0, 18.3, 8.9])
set(axes,'YGrid','on')
hold on
histogram(y2,20,'FaceColor',orangeColor, 'EdgeColor', orangeColor)
grid off
box on
hold on
mu = mean(y2);
sigma = std(y2);
y = min(y2)-0.5:0.1:max(y2);
f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
plot(y,100*f,'LineWidth',1.5)
ylabel('$Frequency,\ [-]$','Interpreter','Latex','FontSize', 12)
xlabel('$Air\ Temperature\ at\ noon,\ [^{\circ}C]$','Interpreter','Latex','FontSize', 12)
xlim([min(y2)-0.5 max(y2)])
% exportgraphics(gcf,'GlobalSA_frequency.png','Resolution',600)
fprintf('Output mean: %2f\n', mu)
fprintf('Output standard deviation: %2f\n', sigma)


% Correlations
blueColor =[195 215 232]./255;
blueD =[62 117 166]./255;
figure
set(gcf, 'units','centimeters','Position',  [0, 0, 18.3, 18.3])
subplot(5,4,1, 'align')
scatter(valuesMatrix(:,1),y2,'.','MarkerEdgeColor', blueColor,'MarkerFaceColor', blueColor)
ref = refline;
ref.LineWidth = 1; ref.Color=blueD;
ylabel('$T_{Air},\ [^{\circ}C]$','Interpreter','Latex','FontSize', 12)
xlabel('$rsmin$','Interpreter','Latex','FontSize', 12)
xlim([min(valuesMatrix(:,1)) max(valuesMatrix(:,1))])
ylim([min(y2) max(y2)])
box on

subplot(5,4,2, 'align')
scatter(valuesMatrix(:,2),y2,'.','MarkerEdgeColor', blueColor,'MarkerFaceColor', blueColor)
ref = refline;
ref.LineWidth = 1; ref.Color=blueD;
yticklabels({''})
xlabel('$Cd$','Interpreter','Latex','FontSize', 12)
xlim([min(valuesMatrix(:,2)) max(valuesMatrix(:,2))])
ylim([min(y2) max(y2)])
box on

subplot(5,4,3, 'align')
scatter(valuesMatrix(:,3),y2,'.','MarkerEdgeColor', blueColor,'MarkerFaceColor', blueColor)
ref = refline;
ref.LineWidth = 1; ref.Color=blueD;
yticklabels({''})
xlabel('$r_{net}$','Interpreter','Latex','FontSize', 12)
xlim([min(valuesMatrix(:,3)) max(valuesMatrix(:,3))])
ylim([min(y2) max(y2)])
box on

subplot(5,4,4, 'align')
scatter(valuesMatrix(:,4),y2,'.','MarkerEdgeColor', blueColor,'MarkerFaceColor', blueColor)
ref = refline;
ref.LineWidth = 1; ref.Color=blueD;
yticklabels({''})
xlabel('$Aside_{vent}$','Interpreter','Latex','FontSize', 12)
xlim([min(valuesMatrix(:,4)) max(valuesMatrix(:,4))])
ylim([min(y2) max(y2)])
box on

subplot(5,4,5, 'align')
scatter(valuesMatrix(:,5),y2,'.','MarkerEdgeColor', blueColor,'MarkerFaceColor', blueColor)
ref = refline;
ref.LineWidth = 1; ref.Color=blueD;
yticklabels({''})
xlabel('$\tau_{NIR,ThScr}$','Interpreter','Latex','FontSize', 12)
xlim([min(valuesMatrix(:,5)) max(valuesMatrix(:,5))])
ylim([min(y2) max(y2)])
box on

subplot(5,4,6, 'align')
scatter(valuesMatrix(:,6),y2,'.','MarkerEdgeColor', blueColor,'MarkerFaceColor', blueColor)
ref = refline;
ref.LineWidth = 1; ref.Color=blueD;
yticklabels({''})
xlabel('$a_{FIR, Sky\&Air}$','Interpreter','Latex','FontSize', 12)
xlim([min(valuesMatrix(:,6)) max(valuesMatrix(:,6))])
ylim([min(y2) max(y2)])
box on

subplot(5,4,7, 'align')
scatter(valuesMatrix(:,7),y2,'.','MarkerEdgeColor', blueColor,'MarkerFaceColor', blueColor)
ref = refline;
ref.LineWidth = 1; ref.Color=blueD;
yticklabels({''})
xlabel('$cout1$','Interpreter','Latex','FontSize', 12)
xlim([min(valuesMatrix(:,7)) max(valuesMatrix(:,7))])
ylim([min(y2) max(y2)])
box on

subplot(5,4,8, 'align')
scatter(valuesMatrix(:,8),y2,'.','MarkerEdgeColor', blueColor,'MarkerFaceColor', blueColor)
ref = refline;
ref.LineWidth = 1; ref.Color=blueD;
yticklabels({''})
xlabel('$cin$','Interpreter','Latex','FontSize', 12)
xlim([min(valuesMatrix(:,8)) max(valuesMatrix(:,8))])
ylim([min(y2) max(y2)])
box on

subplot(5,4,9, 'align')
scatter(valuesMatrix(:,9),y2,'.','MarkerEdgeColor', blueColor,'MarkerFaceColor', blueColor)
ref = refline;
ref.LineWidth = 1; ref.Color=blueD;
xlim([min(valuesMatrix(:,9)) max(valuesMatrix(:,9))])
yticklabels({''})
xlabel('$htop$','Interpreter','Latex','FontSize', 12)
ylim([min(y2) max(y2)])
box on

subplot(5,4,10, 'align')
scatter(valuesMatrix(:,10),y2,'.','MarkerEdgeColor', blueColor,'MarkerFaceColor', blueColor)
ref = refline;
ref.LineWidth = 1; ref.Color=blueD;
yticklabels({''})
xlabel('$\theta_{w}$','Interpreter','Latex','FontSize', 12)
xlim([min(valuesMatrix(:,10)) max(valuesMatrix(:,10))])
ylim([min(y2) max(y2)])
box on

subplot(5,4,11, 'align')
scatter(valuesMatrix(:,11),y2,'.','MarkerEdgeColor', blueColor,'MarkerFaceColor', blueColor)
ref = refline;
ref.LineWidth = 1; ref.Color=blueD;
yticklabels({''})
xlabel('$\tau_{FIR, Cov}$','Interpreter','Latex','FontSize', 12)
xlim([min(valuesMatrix(:,11)) max(valuesMatrix(:,11))])
box on

subplot(5,4,12, 'align')
scatter(valuesMatrix(:,12),y2,'.','MarkerEdgeColor', blueColor,'MarkerFaceColor', blueColor)
ref = refline;
ref.LineWidth = 1; ref.Color=blueD;
yticklabels({''})
xlabel('$n$','Interpreter','Latex','FontSize', 12)
xlim([min(valuesMatrix(:,12)) max(valuesMatrix(:,12))])
box on

subplot(5,4,13, 'align')
scatter(valuesMatrix(:,13),y2,'.','MarkerEdgeColor', blueColor,'MarkerFaceColor', blueColor)
ref = refline;
ref.LineWidth = 1; ref.Color=blueD;
yticklabels({''})
xlabel('$\tau_{FIR, ThScr}$','Interpreter','Latex','FontSize', 12)
xlim([min(valuesMatrix(:,13)) max(valuesMatrix(:,13))])
ylim([min(y2) max(y2)])
box on

subplot(5,4,14, 'align')
scatter(valuesMatrix(:,14),y2,'.','MarkerEdgeColor', blueColor,'MarkerFaceColor', blueColor)
ref = refline;
ref.LineWidth = 1; ref.Color=blueD;
xlim([min(valuesMatrix(:,14)) max(valuesMatrix(:,14))])
yticklabels({''})
xlabel('$a_{FIR,Cov}$','Interpreter','Latex','FontSize', 12)
ylim([min(y2) max(y2)])
box on

subplot(5,4,15, 'align')
scatter(valuesMatrix(:,15),y2,'.','MarkerEdgeColor', blueColor,'MarkerFaceColor', blueColor)
ref = refline;
ref.LineWidth = 1; ref.Color=blueD;
yticklabels({''})
xlabel('$Cst$','Interpreter','Latex','FontSize', 12)
xlim([min(valuesMatrix(:,15)) max(valuesMatrix(:,15))])
ylim([min(y2) max(y2)])
box on

subplot(5,4,16, 'align')
scatter(valuesMatrix(:,18),y2,'.','MarkerEdgeColor', blueColor,'MarkerFaceColor', blueColor)
ref = refline;
ref.LineWidth = 1; ref.Color=blueD;
yticklabels({''})
xlabel('$T_{min\_night}$','Interpreter','Latex','FontSize', 12)
xlim([min(valuesMatrix(:,18)) max(valuesMatrix(:,18))])
ylim([min(y2) max(y2)])
box on

subplot(5,4,17, 'align')
scatter(valuesMatrix(:,16),y2,'.','MarkerEdgeColor', blueColor,'MarkerFaceColor', blueColor)
ref = refline;
ref.LineWidth = 1; ref.Color=blueD;
yticklabels({''})
xlabel('$\tau^{dif}_{PAR,ThScr}$','Interpreter','Latex','FontSize', 12)
xlim([min(valuesMatrix(:,16)) max(valuesMatrix(:,16))])
ylim([min(y2) max(y2)])
box on

subplot(5,4,18, 'align')
scatter(valuesMatrix(:,17),y2,'.','MarkerEdgeColor', blueColor,'MarkerFaceColor', blueColor)
ref = refline;
ref.LineWidth = 1; ref.Color=blueD;
xlim([min(valuesMatrix(:,17)) max(valuesMatrix(:,17))])
yticklabels({''})
xlabel('$\tau^{dif}_{PAR,Cov}$','Interpreter','Latex','FontSize', 12)
ylim([min(y2) max(y2)])
box on

subplot(5,4,19, 'align')
scatter(valuesMatrix(:,19),y2,'.','MarkerEdgeColor', blueColor,'MarkerFaceColor', blueColor)
ref = refline;
ref.LineWidth = 1; ref.Color=blueD;
yticklabels({''})
set(gca, 'XTick', [1910 2112 2300], 'XTickLabel', {'1910', '2112' , '2300'})
xlabel('$C_{abs}$','Interpreter','Latex','FontSize', 12)
xlim([min(valuesMatrix(:,19)) max(valuesMatrix(:,19))])
ylim([min(y2) max(y2)])
box on

subplot(5,4,20, 'align')
scatter(valuesMatrix(:,20),y2,'.','MarkerEdgeColor', blueColor,'MarkerFaceColor', blueColor)
ref = refline;
ref.LineWidth = 1; ref.Color=blueD;
yticklabels({''})
xlabel('$a_{NIR,Wall,i}$','Interpreter','Latex','FontSize', 12)
xlim([min(valuesMatrix(:,20)) max(valuesMatrix(:,20))])
ylim([min(y2) max(y2)])
box on

% exportgraphics(gcf,'GlobalSA_scatter.png','Resolution',600)


%% Parameter estimation plots
% % load PE_datav2.mat
% startDate = datenum( '10-Oct-2019 00:00:00');
% endDate = datenum( '11-Oct-2019 23:59:00');
% TimeData = linspace(startDate,endDate,193);
% orangeColor = [210 109 25]./255;
% blueColor =[62 117 166]./255;
% 
% % Temperature 
% figure
% set(gcf, 'units','centimeters','Position',  [0, 0, 18.3, 8.9])
% plot(TimeData,xe(1:193,1),'-','Color', blueColor,'LineWidth',1); hold on
% plot(TimeData,x(1:193,1),'--','Color',orangeColor ,'LineWidth',1)
% plot(TimeData,xreal(:,1),':k','LineWidth',1)
% %title('Air Temperature')
% xlabel('$Time, [H]','FontSize', 12)
% ylabel('$Temperature, [°C]','FontSize', 12)
% datetick('x','HH','keepticks')
% legend ({'Callibrated','Uncallibrated','Real data'},...
%         'Location','eastoutside')
% exportgraphics(gcf,'PE_temperature.png','Resolution',600)
% 
% % Vapour pressure 
% vpsat_real      = vapourPsat(xreal(:,1));
% vp_real         = vpsat_real.*xreal(:,2);
% vpsat_original  = vapourPsat(x(:,1));
% vp_original     = vpsat_original.*x(:,2);
% vpsat_e         = vapourPsat(xe(:,1));
% vp_e            = vpsat_e.*xe(:,2);
% figure
% set(gcf, 'units','centimeters','Position',  [0, 0, 18.3, 8.9])
% plot(TimeData, vp_original(1:193),'-','Color',blueColor ,'LineWidth',1); hold on
% plot(TimeData, vp_e(1:193),'--', 'Color', orangeColor,'LineWidth',1)
% plot(TimeData,vp_real,':k','LineWidth',1)
% %title('Vapour Pressure')
% xlabel('$Time, [H]','FontSize', 12)
% ylabel('$Vapour Pressure, [kPa]','FontSize', 12)
% datetick('x','HH','keepticks')
% legend ({'Callibrated','Uncallibrated','Real data'},...
%         'Location','eastoutside')
% % exportgraphics(gcf,'PE_rh.png','Resolution',600)
% 
% % Soil temperature
% figure
% set(gcf, 'units','centimeters','Position',  [0, 0, 18.3, 8.9])
% plot(TimeData,xe(1:193,6),'-','Color', blueColor,'LineWidth',1); hold on
% plot(TimeData,x(1:193,6),'--','Color', orangeColor, 'LineWidth',1)
% plot(TimeData,xreal(:,3),':k','LineWidth',1.5)
% %title('Soil Temperature')
% xlabel('$Time, [H]','FontSize', 12)
% ylabel('$Soil Temperature, [°C]','FontSize', 12)
% datetick('x','HH','keepticks')
% legend ({'Callibrated','Uncallibrated','Real data'},...
%         'Location','eastoutside')
% % exportgraphics(gcf,'PE_soilT.png','Resolution',600)

%% Plot the sensitivity change over a day for each parameter
load S20.mat
load S_temperature.mat
load xnom.mat
startDate = datenum( '7-Oct-2019 00:00:00');
endDate = datenum( '7-Oct-2019 23:59:00');
TimeData = linspace(startDate,endDate,97);

figure
set(gcf, 'units','centimeters','Position',  [0, 0, 18.3, 8.9])

names = S_temperature.Properties.RowNames;
ParamNames=readtable('ParamNames.xlsx');
greyColor =[203 221 235]./255;
% for i=1:50
% S_temperature.Properties.RowNames{i} = char(ParamNames.ParamName(i));
% end

yyaxis left
for i=1:361 %number of parameters (half)
%     subplot(5,10,i)
    temp =S.(names{i}); % temporal variable with sensitivities from parameter i
    for j=1:4 % relative changes
        plot(TimeData, temp(:,1,j),'-','Color',greyColor, 'LineWidth',0.8)%,...
%             [1:97], ones(97,1).*table2array(S_temperature(i,j)),'-k') %sensitivities for temperature
        hold on
%         xlabel(char(ParamNames.ParamName(i)),'Interpreter','Latex','FontSize', 12)
    end
end
ylim([-0.5 0.5])
ylabel('$\bar{S},\ [-]$','Interpreter','Latex','FontSize', 12)
datetick('x','HH','keepticks')

yyaxis right
plot(TimeData,x_nom(:,1),'-','Color', orangeColor,'LineWidth',1)
ylim([15 32])
ylabel('$T_{Air},\ [^{\circ}C]$','Interpreter','Latex','FontSize', 12)
datetick('x','HH','keepticks')
xlabel('$Time,\ [h]$','Interpreter','Latex','FontSize', 12)
% ylim([-0.5 0.5])
%% Cumulative areas Si STi
startDate = datenum( '7-Oct-2019 00:00:00');
endDate = datenum( '7-Oct-2019 23:59:00');
TimeData = linspace(startDate,endDate,108);

load GlobalSA_cumulative.mat
% First order sensitivities
% Si = [1 5 3 4 6;...
%     3 2 7 3 2;...
%     1 5 3 4 6];

figure
set(gcf, 'units','centimeters','Position',  [0, 0, 18.3, 8.9])
a = area(TimeData, Si_t);
datetick('x','HH','keepticks')
xlabel('$Time,\ [h]$','Interpreter','Latex','FontSize', 12)
ylabel('$S_i,\ [-]$','Interpreter','Latex','FontSize', 12)
ylim([0 0.25])
box on

a(1).FaceColor = [90 77 123]./255;
a(2).FaceColor = [111 95 151]./255;
a(3).FaceColor = [156 145 185]./255;
a(4).FaceColor = [188 180 208]./255;
a(5).FaceColor = [211 206 224]./255;

legend(names,'Location', 'eastoutside')

%Total effect sensitivities
% STi = [1 5 4 4 6;...
%     3 1 2 3 2;...
%     1 7 3 8 6];
figure
set(gcf, 'units','centimeters','Position',  [0, 0, 18.3, 8.9])
a = area(TimeData, STi_t);
datetick('x','HH','keepticks')
xlabel('$Time,\ [h]$','Interpreter','Latex','FontSize', 12)
ylabel('$S_{T_i},\ [-]$','Interpreter','Latex','FontSize', 12)
ylim([0 1.2])
box on

a(1).FaceColor = [90 77 123]./255;
a(2).FaceColor = [111 95 151]./255;
a(3).FaceColor = [156 145 185]./255;
a(4).FaceColor = [188 180 208]./255;
a(5).FaceColor = [211 206 224]./255;

legend(names,'Location', 'eastoutside')
%% scatter Beta vs Sbar parameters
ParamTable = readtable('ParamsScatter.xlsx');
values = table2array(ParamTable(:,2:3));
distance = sqrt((values(:,1).^2)+(values(:,2).^2));
ParamTable.distance=distance;
ParamTable=sortrows(ParamTable,'distance','descend');
values = table2array(ParamTable(:,2:3));

blueColor =[195 215 232]./255;
blueD =[62 117 166]./255;

figure
set(gcf, 'units','centimeters','Position',  [0, 0, 18.3, 8.9])
scatter(values(1:15,1),values(1:15,2),20)
xline(0,'Color', [0.8,0.8,0.8])
yline(0,'Color', [0.8,0.8,0.8])
xlim([-0.5 0.7])
% grid on
% ylim([-0.6 0.6])
axis equal
ylabel('$\hat{S},\ [-]$','Interpreter','Latex','FontSize', 12)
xlabel('$\hat{\beta},\ [-]$','Interpreter','Latex','FontSize', 12)
box on
text( values(1:15,1)+0.01,values(1:15,2)+0.02, string(1:1:15),'Color',blueD, 'fontweight', 'bold' )

%% Boxplot 














%% create a table like in saltelli
tstep=10;
S_Tab = table;
S_Tab.Time=datetime(TimeData(1:tstep:end),'ConvertFrom','datenum','Format','HH')';

for i=1:3
    S_Tab.('S_{'+string(names(i))+'}')=Si_t(1:tstep:108,i);
    S_Tab.('ST_{'+string(names(i))+'}')=STi_t(1:tstep:108,i);
    S_Tab.('ST_{'+string(names(i))+'}-''S_{'+string(names(i))+'}')=STi_t(1:tstep:108,i)-Si_t(1:tstep:108,i);
    
% S_Tab.S_tauNIRThScr= Si_t(1:10:108,1);
% S_Tab.ST_tauNIRThScr= STi_t(1:10:108,1);
% S_Tab.S_Tminday= Si_t(1:10:108,2);
% S_Tab.ST_Tminday= STi_t(1:10:108,2);
% S_Tab.S_Cd= Si_t(1:10:108,3);
% S_Tab.ST_Cd= STi_t(1:10:108,3);
% S_Tab.('hi')=[1:11]'
end


%% Parameter space plot
% A comparison of the parameter space proved in the LSA and GSA 
% (two parameters)
load('GlobalSA_varMethod_latin5_A.mat')
load 'S_temperature.mat'
load 'xnom.mat'
load 'S20.mat'
%colors
orangeColor = [234 150 81]./255;
blueColor =[195 215 232]./255;
blueD =[62 117 166]./255;
greyColor =[203 221 235]./255;
% Time data
startDate = datenum( '7-Oct-2019 00:00:00');
endDate = datenum( '7-Oct-2019 23:59:00');
TimeData = linspace(startDate,endDate,97);

Cd_nom = 0.75;
r_net_nom = 0.6;

Cd_LSA =[0.8 0.9 1.1 1.2].*Cd_nom;
r_net_LSA = [0.8 0.9 1.1 1.2].*r_net_nom;

Cd_GSA = A(:,2);
r_net_GSA = A(:,3);

% figure

% Local SA parameter space
figure
set(gcf, 'units','centimeters','Position',  [0, 0, 9.15, 9.15])

scatter(ones(1,4).*r_net_nom, Cd_LSA, 'MarkerEdgeColor', blueD)
hold on;
scatter(r_net_LSA,ones(1,4).*Cd_nom, 'MarkerEdgeColor', blueD)
plot(r_net_nom,Cd_nom,'s','MarkerSize',10,...
    'MarkerEdgeColor',orangeColor,...
    'MarkerFaceColor',orangeColor)
xlabel('$r_{net},\ [-]$', 'Interpreter', 'Latex','FontSize', 12)
ylabel('$Cd,\ [-]$', 'Interpreter', 'Latex','FontSize', 12)
xlim([0.3 0.72])
ylim([0 1])
box on

% Global SA parameter space
figure
set(gcf, 'units','centimeters','Position',  [0, 0, 9.15, 9.15])
scatter(Cd_GSA,r_net_GSA,'MarkerEdgeColor', greyColor)
hold on
plot(r_net_nom,Cd_nom,'s','MarkerSize',10,...
    'MarkerEdgeColor',orangeColor,...
    'MarkerFaceColor',orangeColor)
xlabel('$r_{net},\ [-]$', 'Interpreter', 'Latex','FontSize', 12)
ylabel('$Cd,\ [-]$', 'Interpreter', 'Latex','FontSize', 12)
xlim([0.3 0.72])
ylim([0 1])
box on

%Model outputs
delta =[0.8 0.9 1.1 1.2];
Y_Cd = (delta.*squeeze(S.Cd(:,1,:))).*(x_nom(:,1).*(delta-1))+x_nom(:,1);
% Y_Cd=(delta.*table2array(S_temperature({'Cd'},:))).*(x_nom(:,1).*(delta-1))+x_nom(:,1);
% Y_Cd=x_nom(:,1).*ones(1,4)+x_nom(:,1).*Y_Cd;
Y_rnet=(delta.*squeeze(S.r_net(:,1,:))).*(x_nom(:,1).*(delta-1))+x_nom(:,1);
% Y_rnet=x_nom(:,1).*ones(1,4)+x_nom(:,1).*Y_rnet;


figure
set(gcf, 'units','centimeters','Position',  [0, 0, 9.15, 9.15])
plot(([TimeData].*ones(4,1))', Y_Cd,'Color',blueColor); hold on
plot(([TimeData].*ones(4,1))', Y_rnet,'LineStyle','-','Color',blueColor,'LineWidth',1.5)
plot(TimeData, x_nom(:,1),'-','Color',orangeColor,'LineWidth',1.5)
datetick('x','HH','keepticks')
xlabel('$Time,\ [h]$','Interpreter','Latex','FontSize', 12)
ylabel('$Temperature, \ [^\circ C]$','Interpreter','Latex','FontSize', 12)
ylim([10 37])
box on

vq1 = interp1(1:108,y_A(:,:,1)',linspace(1,108,97));
figure
set(gcf, 'units','centimeters','Position',  [0, 0, 9.15, 9.15])
plot(TimeData, vq1,'Color', blueColor,'LineWidth',1.5); hold on
plot(TimeData, x_nom(:,1),'-','Color',orangeColor,'LineWidth',1.5)
datetick('x','HH','keepticks')
xlabel('$Time,\ [h]$','Interpreter','Latex','FontSize', 12)
ylabel('$Temperature, \ [^\circ C]$','Interpreter','Latex','FontSize', 12)
ylim([10 37])
box on


