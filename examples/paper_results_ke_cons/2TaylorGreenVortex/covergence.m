clear
format longG


string="10";
filename = "results/"+string+"subitr.csv";
startRow = 2;
formatSpec = '%24f%24f%24f%24f%24f%24f%24f%24f%24f%f%[^\n\r]';
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


filename = "results/default_20sub.csv";
startRow = 3;
formatSpec = '%24f%24f%24f%24f%24f%24f%24f%24f%24f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
subitr=abs([dataArray{1:end-1}]);
time_de20=subitr(:,1);
dt_de20=subitr(:,2);
dKEdt_de20=subitr(:,3);
KE_de20=subitr(:,4);
VF_de20=subitr(:,5);
rho_de20=subitr(:,6);
rhoU_de20=subitr(:,7);
rhoV_de20=subitr(:,8);


blue=[71,135,224]/255;
red =[218,62,32]/255;
green= [82,147,47]/255;
yellow= [241,178,56]/255;
pink=[149,20,83]/255;

KE0=KE_20(1);
VF0=78.956835208714452 ;
rho0=79125.928586947819;
rhoU0=1.0;
rhoV0=1.0;
umax = 1;
dx=2*pi/64;

%%
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
ylim([10^(-17),10^(-10)])
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
ylim([5*10^(-18),5*10^(-14)])
set(gca,'Fontsize',20)
set(gca,'fontname','Times New Roman')
set(gca,'LineWidth',2.0)
%%
mvwindow=1;
LW2=2.5;
numbertoplot=1200;
start = 400;
figure
plot(time_2(800:numbertoplot:end)*umax/dx,movmean(KE_2(800:numbertoplot:end)/KE0,mvwindow),'LineWidth',LW2,'LineStyle','-','Color',blue,'Marker','diamond','MarkerSize',10,'MarkerFaceColor',blue)
hold on
plot(time_10(start:numbertoplot:end)*umax/dx,movmean(KE_10(start:numbertoplot:end)/KE0,mvwindow),'LineWidth',LW2,'LineStyle','-','Color',blue,'Marker','square','MarkerSize',10,'MarkerFaceColor',blue)
hold on
plot(time_20(1:numbertoplot:end)*umax/dx,movmean(KE_20(1:numbertoplot:end)/KE0,mvwindow),'LineWidth',LW2,'LineStyle','-','Color',blue,'Marker','o','MarkerSize',10,'MarkerFaceColor',blue)
hold off 
xlim([0,80])
ylim([0.999,1.001])
set(gca,'Fontsize',20)
set(gca,'fontname','Times New Roman')
set(gca,'LineWidth',2.0)


%%
figure
plot(time_de20*umax/dx,movmean(KE_de20/KE0,mvwindow),'LineWidth',LW2,'LineStyle','-','Color','k')
xlim([0,80])
ylim([0.0,1.01])
set(gca,'Fontsize',20)
set(gca,'fontname','Times New Roman')
set(gca,'LineWidth',2.0)
