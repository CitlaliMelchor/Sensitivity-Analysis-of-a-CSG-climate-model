%% Useful graphs
clear all
clc
clf
% Beijing historical weather data [www.weatherbase.com]

Months = {'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A','S', 'O', 'N', 'D','J'};
T_avg = [-3, 0, 6,13,20,24,26,25,20,13,5,-1,-3]; % Average temperature,[°C]
Rad_avg=[9.9, 13.2, 17.1, 21, 22.6, 21.7, 19, 16.4, 16.3, 13.3, 10.1, 8.6, 9.9]; % Average solar radiation,[MJ m^-2]
RH_avg=[33,32,32,33,40,47,62,63,51,47,47,42,33]; %Average relative humidity,[%]

% Amsterdam historical weather data [www.weatherbase.com]

T_avgNL = [2.6, 3, 5.8, 8.2, 12.6,  15, 17.3, 17.1, 14.2, 10.3, 6.2, 3.8, 2.6]; % Average temperature,[°C]
Rad_avgNL=[2.8, 5.3, 9.4, 14.8, 19.8, 20, 19.5, 17, 10.8, 6, 3.1, 1.9, 2.8]; % Average solar radiation,[MJ m^-2]
RH_avgNL=[]; %Average relative humidity,[%]


% Plot
figure(1)
set(gcf,'Position',[100 100 800 500])
bluish = [0, 0.4470, 0.7410];
orange = [0.8500, 0.3250, 0.0980];
plot(T_avg,Rad_avg, '-s','Color',orange, 'LineWidth',1,'MarkerFaceColor', orange)
text(T_avg,Rad_avg,Months,'VerticalAlignment','bottom','HorizontalAlignment','right', 'Color', orange)
hold on
plot(T_avgNL,Rad_avgNL, '--o','Color',bluish, 'LineWidth',1,'MarkerFaceColor', bluish)
text(T_avgNL,Rad_avgNL,Months,'VerticalAlignment','bottom','HorizontalAlignment','right', 'Color', bluish)
% grid on
xlabel('$Temperature,\ [^\circ C]$' ,'Interpreter','Latex','FontSize', 12)
ylabel('$Radiation,\ [MJ m^{-2}day^{-1}]$','Interpreter','Latex','FontSize', 12)
% plot requirements
area([-5,8],[27 27],"FaceAlpha", 0.06,'LineStyle','none','FaceColor',[0.9290, 0.6940, 0.1250])
area([8,22],[27 27],"FaceAlpha", 0.06,'LineStyle',':','FaceColor',[0.3010, 0.7450, 0.9330])
area([22,30],[27 27],"FaceAlpha", 0.06,'LineStyle','none','FaceColor',[0, 0.75, 0.75])
text(1.5,24.5,{'Heating','required'}, 'FontSize', 12,'HorizontalAlignment','Center', 'Color',[0.9290, 0.6940, 0.1250])
text(15,24.5,{'Natural','ventilation'}, 'FontSize', 12, 'HorizontalAlignment','Center','Color',[0.3010, 0.7450, 0.9330])
text(26,24.5,{'Cooling','required'}, 'FontSize', 12,'HorizontalAlignment','Center', 'Color',[0, 0.75, 0.75])
ylim([0,27])
legend('Beijing, China','De Bilt, NL','Location','SE', 'FontSize', 10)
set(gcf, 'units','centimeters','Position',  [0, 0, 18.3, 8.9])
exportgraphics(gcf,'BeijingWeatherR.png','Resolution',600)

    