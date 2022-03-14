%% Bootstrap method for Sobol' indices
% MSc Thesis Sensitivity analysis of a Chinese Solar Greenhouse
% Author: Citlali Melchor Ramírez
% December 2021

  clear; clc;
 %%
load GlobalSA_varMethod_latin5v2_A.mat
load GlobalSA_varMethod_latin5v2_B.mat
load GlobalSA_varMethod_latin5v2_C.mat
load GlobalSA_varMethod_latin5v2_D.mat
yAfull = y_A; 
yBfull = y_B; 
yCfull = y_C; 
yDfull = y_D; 
load GlobalSA_varMethod_latin5_A.mat
load GlobalSA_varMethod_latin5_B.mat
load GlobalSA_varMethod_latin5_C.mat
load GlobalSA_varMethod_latin5_D.mat
yAfull = [yAfull ; y_A]; 
yBfull = [yBfull; y_B]; 
yCfull = [yCfull; y_C]; 
yDfull = [yDfull; y_D];
y_A = yAfull;
y_B = yBfull;
y_C = yCfull;
y_D = yDfull;

%% Find and delete null values
[rowA,colA]=find(isnan(y_A(:,:,1,:)));
[rowB,colB]=find(isnan(y_B(:,:,1,:)));
[rowC,colC]=find(isnan(y_C(:,:,1,:)));
[rowD,colD]=find(isnan(y_D(:,:,1,:)));

AreNan=unique([rowA;rowB;rowC;rowD]');

y_A(AreNan,:,:)=[];
y_B(AreNan,:,:)=[];
y_C(AreNan,:,:,:)=[];
y_D(AreNan,:,:,:)=[];

%%
% The bootstrap consists in taking a sample from the data and performing
% the statistic calculations n times, with different samples. 
N = size(y_A,1);
k = 5;
N_prop=N-5; 
S_prop = zeros(1,k);
ST_prop = zeros(1,k); 
Sij_prop = zeros(5,5);

for boot = 1:500000
    Sij_props = zeros(5,5);
    S_est =[];
    ST_est =[];
    Sij_est=[];
    VY_est =[];
    rows = randperm(N,N_prop);% Select 999 random rows     
    for i = 1:108 %time steps               
        % Create the prop matrices 
        y_A_prop = y_A(rows,i,1);%y_A(rows,i);%
        y_B_prop = y_B(rows,i,1);%y_B(rows,i);%
        y_C_prop = squeeze(y_C(rows,i,1,:));%y_C(rows,i,:);%
        y_D_prop = squeeze(y_D(rows,i,1,:));%y_D(rows,i,:);%

        %Calculate first-order indices
    %     f_0sq =((1/N_prop)*sum(y_A_prop))^2; % mean^2
        EY = 1/N_prop*sum(y_A_prop); %unconditional mean Saltelli2004
    %     VY_prop = (1/N_prop)*sum(y_A_prop.^2) - f_0sq; % Variance(Y)
        VY_prop = 1/(N_prop-1)*sum(y_A_prop.^2)-EY^2;  %unconditional variance Saltelli2004
        
        for j = 1:k 
            U = 1/(N_prop-1)*sum(y_A_prop.*y_C_prop(:,j)); 
            VX = U - EY^2;
            S_prop(j) = VX/VY_prop; %Saltelli 2004
    %         S_prop(j) = ((1/N_prop)*sum(y_A_prop.*y_C_prop(:,j))- f_0sq)/VY_prop;
        end
        % 2. Compute the total effect indices 
        for j = 1:k 
            U = 1/(N_prop-1)*sum(y_A_prop.*y_D_prop(:,j)); 
            VX = U - EY^2;
            ST_prop(j) = 1-(VX/VY_prop);
    %         ST_prop(m) = 1-(((1/N_prop)*sum(y_B_prop.*y_C_prop(:,m))- f_0sq)/VY_prop );
        end

        % 3. Compute the second order sensitivities


        for m = 1:k 
            for n = m:k
                Um = 1/(N_prop-1)*sum(y_A_prop.*y_C_prop(:,m)); 
                Vm = Um - EY^2;
                Un = 1/(N_prop-1)*sum(y_B_prop.*y_D_prop(:,n)); 
                Vn = Un - EY^2;
                Umn = 1/(N_prop-1)*sum(y_D_prop(:,m).*y_C_prop(:,n));
                Vmn = Umn - EY^2;
                Sij_props(m,n)= (Vmn-Vm-Vn)/(VY_prop);
            end
        end
        upper  = logical(triu(ones(5)));
        Sij_prop=Sij_props(upper);

        S_est = [S_est, S_prop'];
        ST_est = [ST_est, ST_prop'];
        Sij_est = [Sij_est, Sij_prop];
        VY_est = [VY_est; VY_prop];

    end
    if ~isnan(S_est)
        S_boot(:,:,:,boot) = S_est;
    end
    
    if ~isnan(ST_est)
        ST_boot(:,:,:,boot)  = ST_est;
    end
    
    if ~isnan(Sij_est)
        Sij_boot(:,:,:,boot)  = Sij_est;
    end


end

%%
Si  = squeeze(mean(S_boot,4));
STi = squeeze(mean(ST_boot,4));
Sij = squeeze(mean(Sij_boot,4));

%% Bar plots
t=57;
orangeColor = [234 150 81]./255;
blueColor =[195 215 232]./255;
blueD =[62 117 166]./255;

ParamNames=readtable('ParamBoundsFive.xlsx');
ParamNames.Si=Si(:,t);
ParamNames.STi=STi(:,t);
ParamNames.Simean=mean(Si,2);
ParamNames = sortrows(ParamNames,{'Si'},{'descend'});

figure
set(gcf, 'units','centimeters','Position',  [0, 0, 18.3, 8.9])
ax = nexttile;
b=bar(ax,[ParamNames.Si(1:5),ParamNames.STi(1:5)],'stacked');
xticks([1:1:5])
xticklabels(ParamNames.ParamName(1:5))
xtickangle(45)
b(1).FaceColor = 'flat';
b(1).CData = blueD;
b(2).FaceColor = 'flat';
b(2).CData = blueColor;
ylabel('S_i, S_{T_i}')
% ylim([0 0.55])
%%
names = ParamNames.ParamName; 
Si_t=Si';
STi_t=STi';
Si_t = [Si_t(:,3), Si_t(:,5),Si_t(:,1),Si_t(:,4),Si_t(:,2)]; %ordered Si 
STi_t = [STi_t(:,3), STi_t(:,5),STi_t(:,1),STi_t(:,4),STi_t(:,2)]; %ordered Si 

%% Cumulative areas Si STi
startDate = datenum( '7-Oct-2019 00:00:00');
endDate = datenum( '7-Oct-2019 23:59:00');
TimeData = linspace(startDate,endDate,108);

% load GlobalSA_cumulative.mat
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
 ylim([0 5.1])
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
%%
%  save BootstrapGSA.mat A B S_est Si ST_est STi VY_est

%% Second-order plot
clear SOrder SOnames
names2=readtable('ParamBoundsFive.xlsx');
names2 = names2.ParamName; 
n=1;
for i=1:5
    for j= i:5
        SOnames(n,:) =  append(names2(i),',',names2(j));
        n=n+1;
    end
end
SOrder=array2table(Sij);
SOrder.parameters = SOnames; 
SOrder = sortrows(SOrder,{'Sij55'},{'descend'});

figure
set(gcf, 'units','centimeters','Position',  [0, 0, 20.3, 8.9])
a=area(TimeData,table2array(SOrder(1:5,1:1:108))');
datetick('x','HH','keepticks')
xlabel('$Time,\ [h]$','Interpreter','Latex','FontSize', 12)
ylabel('$S_{ij},\ [-]$','Interpreter','Latex','FontSize', 12)
a(1).FaceColor = [90 77 123]./255;
a(2).FaceColor = [111 95 151]./255;
a(3).FaceColor = [156 145 185]./255;
a(4).FaceColor = [188 180 208]./255;
a(5).FaceColor = [211 206 224]./255;
legend(SOrder.parameters(1:5),'Location', 'eastoutside','Orientation','Vertical')
% xlim([TimeData(1) TimeData(107)])
ylim([-5.1 0.8])

a2 = axes();
a2.Position = [0.4 0.25 0.27 0.32]; % xlocation, ylocation, xsize, ysize
a2=area(a2,TimeData(40:1:75),table2array(SOrder(1:5,40:1:75))'); axis tight
datetick('x','HH','keepticks')
a2(1).FaceColor = [90 77 123]./255;
a2(2).FaceColor = [111 95 151]./255;
a2(3).FaceColor = [156 145 185]./255;
a2(4).FaceColor = [188 180 208]./255;
a2(5).FaceColor = [211 206 224]./255;
ylim([0 0.6])

%% Bar plots
figure 
% subplot (2,2,[1,3])
n = 5;
t=57;
orangeColor = [234 150 81]./255;
blueColor =[195 215 232]./255;
blueD =[62 117 166]./255;

SecOrder = [(SOrder.Sij57(3)+ SOrder.Sij57(4)),...
    (SOrder.Sij57(2)+ SOrder.Sij57(3))+SOrder.Sij57(5),...
    0,...
    SOrder.Sij57(5),...
   (SOrder.Sij57(1)+ SOrder.Sij57(2)+ SOrder.Sij57(4))];

x1 = [ParamNames.Si(1:5),SecOrder']';
x2 = [ParamNames.STi(1:5)]';

b1 = bar(x1','stacked');
hold on;
b1(1).BarWidth = .4;
b1(1).XData = (1:n) + .2; % move stacked bars right
b1(2).XData = (1:n) + .2;
b2 = bar(x2);
b2.BarWidth = 0.4;
b2.XData = (1:n) - .2; % move single bars left
xticks([1:1:5])
xticklabels(ParamNames.ParamName(1:5))
xtickangle(45)

b1(1).FaceColor = 'flat';
b1(1).CData = blueD;
b1(2).FaceColor = 'flat';
b1(2).CData = blueColor;

% subplot(2,2,4)
% SOsum=sum(SOrder.Sij57(1:5));
% pie([SOsum,1-SOsum])
% colormap([blueD;      %// red
%           .5 .5 .5])  %// grey
% subplot(2,2,2)
% Sisum=sum(ParamNames.Si(1:5));
% pie([Sisum, 1-Sisum])
% % colormap([blueColor;      %// red
% %           .5 .5 .5])  %// grey
% hPieComponentHandles = pie([Sisum, 1-Sisum]);%pie(ax, ones(1,2));
% set(hPieComponentHandles(k*2-1), 'FaceColor', [blueColor;      %// red
%           .5 .5 .5]); 
figure
set(gcf, 'units','centimeters','Position',  [0, 0, 18.3, 8.9])
ax = nexttile;
b=bar(ax,[ParamNames.Si(1:5),ParamNames.STi(1:5),SecOrder'],'stacked');
xticks([1:1:5])
xticklabels(ParamNames.ParamName(1:5))
xtickangle(45)
b(1).FaceColor = 'flat';
b(1).CData = blueD;
b(2).FaceColor = 'flat';
b(2).CData = blueColor;
ylabel('S_i, S_{T_i}')
% ylim([0 0.55])
%%
save GlobalSA_cumulative_final.mat Si_t STi_t names SOrder 
save('GlobalSA_cumulative_Sboot.mat', 'ST_boot', '-v7.3')
save('GlobalSA_cumulative_STboot.mat', 'ST_boot', '-v7.3')
save('GlobalSA_cumulative_Sijboot.mat', 'ST_boot', '-v7.3')
