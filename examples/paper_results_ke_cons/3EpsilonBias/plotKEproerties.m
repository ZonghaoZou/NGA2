%% getting from default case
itr = "2sub";
type = "KE";
string = "1em1";
filename = "results/"+type+"_"+itr+"_"+string+".csv";
startRow = 3;
formatSpec = '%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
simulation = [dataArray{1:end-1}];
clearvars filename startRow formatSpec fileID dataArray ans;
% close all
start_index=2;
time_default=simulation(start_index:end,1);
KEdiff_default=abs(simulation(start_index:end,3));
KE_default=abs(simulation(start_index:end,4));

string = "1em3";
filename = "results/"+type+"_"+itr+"_"+string+".csv";
startRow = 3;
formatSpec = '%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
simulation = [dataArray{1:end-1}];
clearvars filename startRow formatSpec fileID dataArray ans;
start_index=2;
time_1em3=simulation(start_index:end,1);
KEdiff_1em3=abs(simulation(start_index:end,3));
KE_1em3=abs(simulation(start_index:end,4));


string = "1em2";
filename = "results/"+type+"_"+itr+"_"+string+".csv";
startRow = 3;
formatSpec = '%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
simulation = [dataArray{1:end-1}];
clearvars filename startRow formatSpec fileID dataArray ans;
% close all
start_index=2;
time_1em2=simulation(start_index:end,1);
KEdiff_1em2=abs(simulation(start_index:end,3));
KE_1em2=abs(simulation(start_index:end,4));

filename = "results/conservation_1000_old";
startRow = 2;
formatSpec = '%24f%24f%24f%24f%24f%24f%24f%24f%24f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
subitr=abs([dataArray{1:end-1}]);
time_dnga=subitr(:,2);
KE_dnga=subitr(:,7);
%% Plotting
figure('Color','w')

LW1=1.2;
LW2=2.5;
blue=[71,135,224]/255;
red =[218,62,32]/255;
green= [82,147,47]/255;
yellow= [241,178,56]/255;
pink=[149,20,83]/255;
u=6;
dx=2*pi/16;

plot(time_default*u/dx,KE_default/(KE_default(1)),'LineWidth',LW2,'LineStyle','--','color',blue)
hold on
plot(time_1em3*u/dx,KE_1em3/(KE_default(1)),'LineWidth',LW2,'LineStyle','-.','color',blue)
hold on
plot(time_1em2*u/dx,KE_1em2/(KE_default(1)),'LineWidth',LW2,'LineStyle',':','color',blue)
hold on
plot(time_dnga,KE_dnga,'LineWidth',LW2,'LineStyle','-','color','k')
xlim([0,80])
ylim([0,1.1])
set(gca,'Fontsize',20)
set(gca,'fontname','Times New Roman')
set(gca,'LineWidth',2)
%%

figure 
ratio =2.5;
plot(time_default*u/dx,KE_default/(KE_default(1)),'LineWidth',LW2*ratio,'LineStyle','--','color',blue)
hold on
plot(time_1em3*u/dx,KE_1em3/(KE_default(1)),'LineWidth',LW2*ratio,'LineStyle','-.','color',blue)
hold on
plot(time_1em2*u/dx,KE_1em2/(KE_default(1)),'LineWidth',LW2*ratio,'LineStyle',':','color',blue)
xlim([70,80])
ylim([0.9,1.01])
yticks([0.9 0.925 0.95 0.975 1.0])
set(gca,'Fontsize',20*ratio)
set(gca,'fontname','Times New Roman')
set(gca,'LineWidth',2*ratio)

