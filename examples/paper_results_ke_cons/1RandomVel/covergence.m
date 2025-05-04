clear
format longG

string="10";
filename = "results/"+string+"subitr.csv";
startRow = 2;
formatSpec = '%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
subitr=abs([dataArray{1:end-1}]);
time_10=subitr(:,1);
dt_10=subitr(:,2);
dKEdt_10=subitr(:,3);
KE_10=subitr(:,4);
VF_10=subitr(:,5);
rho_10=subitr(:,6);
rhoU_10=subitr(:,7);
rhoV_10=subitr(:,8);

string="20";
filename = "results/"+string+"subitr.csv";
startRow = 2;
formatSpec = '%24f%24f%24f%24f%24f%24f%24f%24f%24f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
subitr=abs([dataArray{1:end-1}]);
time_20=subitr(:,1);
dt_20=subitr(:,2);
dKEdt_20=subitr(:,3);
KE_20=subitr(:,4);
VF_20=subitr(:,5);
rho_20=subitr(:,6);
rhoU_20=subitr(:,7);
rhoV_20=subitr(:,8);



filename = "results/2subitr.csv";
startRow = 2;
formatSpec = '%24f%24f%24f%24f%24f%24f%24f%24f%24f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
subitr=abs([dataArray{1:end-1}]);
time_2=subitr(:,1);
dt_2=subitr(:,2);
dKEdt_2=subitr(:,3);
KE_2=subitr(:,4);
VF_2=subitr(:,5);
rho_2=subitr(:,6);
rhoU_2=subitr(:,7);
rhoV_2=subitr(:,8);

filename = "results/conservation_1000_old";
startRow = 2;
formatSpec = '%24f%24f%24f%24f%24f%24f%24f%24f%24f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
subitr=abs([dataArray{1:end-1}]);
time_de20=subitr(:,2);
KE_de20=subitr(:,7);

filename = "results/conservation_1000";
startRow = 2;
formatSpec = '%24f%24f%24f%24f%24f%24f%24f%24f%24f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
subitr=abs([dataArray{1:end-1}]);
time_de2017=subitr(:,2);
KE_de2017=subitr(:,7);


%%
blue=[71,135,224]/255;
red =[218,62,32]/255;
green= [82,147,47]/255;
yellow= [241,178,56]/255;
pink=[149,20,83]/255;

LW2=1.0;
KE0=KE_20(1);
VF0=3.8757845850374766;
rho0=3887.4119387925894;
rhoU0=-297.66125637821875;
rhoV0=124.13953954228096;
umax=6.6286086170103822;
dx=2*pi/16;

mvwindow=1;
LW2=0.8;
%%
close all
mvwindow=100;
numbertoplot=2; 
LW2=2.0;
figure
semilogy(time_2(1:numbertoplot:end)*umax/dx,movmean(dKEdt_2(1:numbertoplot:end).*dt_2(1:numbertoplot:end)/KE0,mvwindow),'LineWidth',LW2,'LineStyle','-.','Color',blue)
hold on
semilogy(time_20(1:numbertoplot:end)*umax/dx,movmean(dKEdt_20(1:numbertoplot:end).*dt_20(1:numbertoplot:end)/KE0,mvwindow),'LineWidth',LW2,'LineStyle','-','Color',blue)
hold on
semilogy(time_10(1:numbertoplot:end)*umax/dx,movmean(dKEdt_10(1:numbertoplot:end).*dt_10(1:numbertoplot:end)/KE0,mvwindow),'LineWidth',LW2,'LineStyle',':','Color',blue)
xlim([0,80])
ylim([0.7*10^(-16),10^(-3)])
set(gca,'Fontsize',20)
set(gca,'fontname','Times New Roman')
set(gca,'LineWidth',2.0)

%% Conservation properties
figure 
mvwindow=100;
numbertoplot=2; 
start=20;
LW2=2.0;
semilogy(time_20(1:numbertoplot:end)*umax/dx,movmean(dKEdt_20(1:numbertoplot:end).*dt_20(1:numbertoplot:end)/KE0,mvwindow),'LineWidth',LW2,'LineStyle','-','Color',blue)
hold on
semilogy(time_20(start:numbertoplot:end)*umax/dx,movmean(VF_20(start:numbertoplot:end).*dt_20(start:numbertoplot:end)/VF0,mvwindow),'LineWidth',LW2,'LineStyle','-','Color',red)
hold on
semilogy(time_20(1:numbertoplot:end)*umax/dx,movmean(rho_20(1:numbertoplot:end).*dt_20(1:numbertoplot:end)/rho0,mvwindow),'LineWidth',LW2,'LineStyle','-','Color',green)
hold on
semilogy(time_20(1:numbertoplot:end)*umax/dx,movmean(rhoU_20(1:numbertoplot:end).*dt_20(1:numbertoplot:end)/abs(rhoU0),mvwindow),'LineWidth',LW2,'LineStyle','-','Color',yellow)
hold on
semilogy(time_20(1:numbertoplot:end)*umax/dx,movmean(rhoV_20(1:numbertoplot:end).*dt_20(1:numbertoplot:end)/abs(rhoV0),mvwindow),'LineWidth',LW2,'LineStyle','-','Color',pink)
xlim([0,80])
ylim([0.7*10^(-16),1.1*10^(-12)])
set(gca,'Fontsize',20)
set(gca,'fontname','Times New Roman')
set(gca,'LineWidth',2.0)
%%
mvwindow=1;
LW2=2.5;
numbertoplot=300;
start = 150;
figure
plot(time_2(1:150:end)*umax/dx,movmean(KE_2(1:150:end)/KE0,mvwindow),'LineWidth',LW2,'LineStyle','-','Color',blue,'Marker','diamond','MarkerSize',10,'MarkerFaceColor',blue)
hold on
plot(time_10(start:numbertoplot:end)*umax/dx,movmean(KE_10(start:numbertoplot:end)/KE0,mvwindow),'LineWidth',LW2,'LineStyle','-','Color',blue,'Marker','square','MarkerSize',10,'MarkerFaceColor',blue)
hold on
plot(time_20(1:numbertoplot:end)*umax/dx,movmean(KE_20(1:numbertoplot:end)/KE0,mvwindow),'LineWidth',LW2,'LineStyle','-','Color',blue,'Marker','o','MarkerSize',10,'MarkerFaceColor',blue)
hold off 
yticks([0.997 0.998 0.999 1.0 1.001])
xlim([0,80])
ylim([0.997,1.001])
set(gca,'Fontsize',20)
set(gca,'fontname','Times New Roman')
set(gca,'LineWidth',2.0)
%%
mvwindow =1; 
figure
plot(time_de20,movmean(KE_de20,mvwindow),'LineWidth',LW2,'LineStyle','-','Color','k')
hold on
plot(time_de2017,movmean(KE_de2017,mvwindow),'LineWidth',LW2,'LineStyle','--','Color','k')
xlim([0,80])
ylim([0.0,1.01])
set(gca,'Fontsize',20)
set(gca,'fontname','Times New Roman')
set(gca,'LineWidth',2.0)


