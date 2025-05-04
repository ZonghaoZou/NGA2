%% getting from KE conserving case
clear
type_file="unity";
meshsize="8";
filename="results/simulation_"+type_file+"_"+meshsize;
% filename="monitor/simulation";
startRow = 3;
formatSpec = '%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
simulation = [dataArray{1:end-1}];
clearvars filename startRow formatSpec fileID dataArray ans;
start_index=1;
time_8=simulation(start_index:end-1,2);
Asim_8=abs(simulation(start_index+1:end,15));
% 
meshsize="16";
filename="results/simulation_"+type_file+"_"+meshsize;
startRow = 3;
formatSpec = '%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
simulation = [dataArray{1:end-1}];
clearvars filename startRow formatSpec fileID dataArray ans;
time_16=simulation(start_index:end-1,2);
Asim_16=abs(simulation(start_index+1:end,15));

meshsize="32";
filename="results/simulation_"+type_file+"_"+meshsize;
startRow = 3;
formatSpec = '%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
simulation = [dataArray{1:end-1}];
clearvars filename startRow formatSpec fileID dataArray ans;
time_32=simulation(start_index:end-1,2);
Asim_32=abs(simulation(start_index+1:end,15));
% % 
meshsize="64";
filename="results/simulation_"+type_file+"_"+meshsize;
startRow = 3;
formatSpec = '%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
simulation = [dataArray{1:end-1}];
clearvars filename startRow formatSpec fileID dataArray ans;
start_index=2;
time_64=simulation(start_index:end-1,2);
Asim_64=abs(simulation(start_index+1:end,15));


%% StandingWave Amplitude
close all
blue=[71,135,224]/255;
red =[218,62,32]/255;
green= [82,147,47]/255;
yellow= [241,178,56]/255;
LW1=1.0;
LW2=2.50;

rho1=1.0;
rho2=1.0;
beta=rho1*rho2/((rho1+rho2)^2);
nu=0.064720863;
sigma=2;
omega0=sqrt(sigma/(rho1+rho2));
A0=0.01*2*pi;
r=roots([1.0 -4*beta*sqrt(nu) 2*(1-6*beta)*nu 4*(1-3*beta)*nu^(3/2) (1-4*beta)*nu^2+omega0^2]);
z1=r(1);
z2=r(2);
z3=r(3);
z4=r(4);
Z1 = (z2 - z1)*(z3 - z1)*(z4 - z1);
Z2 = (z1 - z2)*(z3 - z2)*(z4 - z2);
Z3 = (z1 - z3)*(z2 - z3)*(z4 - z3);
Z4 = (z1 - z4)*(z2 - z4)*(z3 - z4);


t = 0:0.3:20;
Aex=((4*(1-4*beta)*nu^2)/(8*(1-4*beta)*nu^2 + omega0^2)*A0*erfc(sqrt(nu*t)) + (z1*A0*omega0^2)/(Z1*(z1^2-nu))*exp((z1^2-nu)*t).*erfc(real(z1)*sqrt(t)) ...
  + (z2*A0*omega0^2)/(Z2*(z2^2-nu))*exp((z2^2-nu)*t).*erfc(real(z2)*sqrt(t)) ...
  + (z3*A0*omega0^2)/(Z3*(z3^2-nu))*exp((z3^2-nu)*t).*erfc(real(z3)*sqrt(t)) ...
  + (z4*A0*omega0^2)/(Z4*(z4^2-nu))*exp((z4^2-nu)*t).*erfc(real(z4)*sqrt(t)) )/(2*pi);

figure('Color','w')
plot(time_8,(Asim_8-pi)/(2*pi),'LineWidth',LW2,'LineStyle','--','Color','k')
hold on
plot(time_16,(Asim_16-pi)/(2*pi),'LineWidth',LW2,'LineStyle',':','Color','k')
hold on
plot(time_32,(Asim_32-pi)/(2*pi),'LineWidth',LW2,'LineStyle','-.','Color','k')
hold on
plot(time_64,(Asim_64-pi)/(2*pi),'LineWidth',LW2,'LineStyle','-','Color','k')
hold on
plot(t,Aex,'LineWidth',LW1,'Marker','.','Markersize',14.0,'Color','k')
xlim([0,20])
ylim([-0.01,0.01])
set(gca,'Fontsize',20)
set(gca,'fontname','Times New Roman')
set(gca,'LineWidth',2.0)

%% StandingWave Error
perturb=0.01;
rho1=1.0;
rho2=1.0;
beta=rho1*rho2/((rho1+rho2)^2);
nu=0.064720863;
sigma=2;
omega0=sqrt(sigma/(rho1+rho2));
A0=perturb*2*pi;
r=roots([1.0 -4*beta*sqrt(nu) 2*(1-6*beta)*nu 4*(1-3*beta)*nu^(3/2) (1-4*beta)*nu^2+omega0^2]);
z1=r(1);
z2=r(2);
z3=r(3);
z4=r(4);
Z1 = (z2 - z1)*(z3 - z1)*(z4 - z1);
Z2 = (z1 - z2)*(z3 - z2)*(z4 - z2);
Z3 = (z1 - z3)*(z2 - z3)*(z4 - z3);
Z4 = (z1 - z4)*(z2 - z4)*(z3 - z4);

t = time_8;
Aex=((4*(1-4*beta)*nu^2)/(8*(1-4*beta)*nu^2 + omega0^2)*A0*erfc(sqrt(nu*t)) + (z1*A0*omega0^2)/(Z1*(z1^2-nu))*exp((z1^2-nu)*t).*erfc(real(z1)*sqrt(t)) ...
  + (z2*A0*omega0^2)/(Z2*(z2^2-nu))*exp((z2^2-nu)*t).*erfc(real(z2)*sqrt(t)) ...
  + (z3*A0*omega0^2)/(Z3*(z3^2-nu))*exp((z3^2-nu)*t).*erfc(real(z3)*sqrt(t)) ...
  + (z4*A0*omega0^2)/(Z4*(z4^2-nu))*exp((z4^2-nu)*t).*erfc(real(z4)*sqrt(t)) )/(2*pi);
error_8= ((Asim_8-pi)/(2*pi)-Aex)/perturb;

t = time_16;
Aex=((4*(1-4*beta)*nu^2)/(8*(1-4*beta)*nu^2 + omega0^2)*A0*erfc(sqrt(nu*t)) + (z1*A0*omega0^2)/(Z1*(z1^2-nu))*exp((z1^2-nu)*t).*erfc(real(z1)*sqrt(t)) ...
  + (z2*A0*omega0^2)/(Z2*(z2^2-nu))*exp((z2^2-nu)*t).*erfc(real(z2)*sqrt(t)) ...
  + (z3*A0*omega0^2)/(Z3*(z3^2-nu))*exp((z3^2-nu)*t).*erfc(real(z3)*sqrt(t)) ...
  + (z4*A0*omega0^2)/(Z4*(z4^2-nu))*exp((z4^2-nu)*t).*erfc(real(z4)*sqrt(t)) )/(2*pi);
error_16= ((Asim_16-pi)/(2*pi)-Aex)/perturb;

t = time_32;
Aex=((4*(1-4*beta)*nu^2)/(8*(1-4*beta)*nu^2 + omega0^2)*A0*erfc(sqrt(nu*t)) + (z1*A0*omega0^2)/(Z1*(z1^2-nu))*exp((z1^2-nu)*t).*erfc(real(z1)*sqrt(t)) ...
  + (z2*A0*omega0^2)/(Z2*(z2^2-nu))*exp((z2^2-nu)*t).*erfc(real(z2)*sqrt(t)) ...
  + (z3*A0*omega0^2)/(Z3*(z3^2-nu))*exp((z3^2-nu)*t).*erfc(real(z3)*sqrt(t)) ...
  + (z4*A0*omega0^2)/(Z4*(z4^2-nu))*exp((z4^2-nu)*t).*erfc(real(z4)*sqrt(t)) )/(2*pi);
error_32= ((Asim_32-pi)/(2*pi)-Aex)/perturb;

t = time_64;
Aex=((4*(1-4*beta)*nu^2)/(8*(1-4*beta)*nu^2 + omega0^2)*A0*erfc(sqrt(nu*t)) + (z1*A0*omega0^2)/(Z1*(z1^2-nu))*exp((z1^2-nu)*t).*erfc(real(z1)*sqrt(t)) ...
  + (z2*A0*omega0^2)/(Z2*(z2^2-nu))*exp((z2^2-nu)*t).*erfc(real(z2)*sqrt(t)) ...
  + (z3*A0*omega0^2)/(Z3*(z3^2-nu))*exp((z3^2-nu)*t).*erfc(real(z3)*sqrt(t)) ...
  + (z4*A0*omega0^2)/(Z4*(z4^2-nu))*exp((z4^2-nu)*t).*erfc(real(z4)*sqrt(t)) )/(2*pi);
error_64= ((Asim_64-pi)/(2*pi)-Aex)/perturb;



skipcell=1;

figure
plot(time_8(1:skipcell:end),error_8(1:skipcell:end),'LineWidth',LW2,'LineStyle','--','Color','k')
hold on
plot(time_16(1:skipcell:end),error_16(1:skipcell:end),'LineWidth',LW2,'LineStyle',':','Color','k')
hold on
plot(time_32(1:skipcell:end),error_32(1:skipcell:end),'LineWidth',LW2,'LineStyle','-.','Color','k')
hold on
plot(time_64(1:skipcell:end),error_64(1:skipcell:end),'LineWidth',LW2,'LineStyle','-','Color','k')

xlim([0,20])
ylim([-0.2,0.2])
set(gca,'Fontsize',20)
set(gca,'fontname','Times New Roman')
set(gca,'LineWidth',2)
%% RMS convergence
% close all
rms_8=rms(error_8);
rms_16=rms(error_16);
rms_32=rms(error_32);
rms_64=rms(error_64);


fprintf('%.15f\n', rms_8)
fprintf('%.15f\n', rms_16)
fprintf('%.15f\n', rms_32)
fprintf('%.15f\n', rms_64)
