clear
format longG
foldername1="1DefaultNGA2_SG";
foldername2="2SpatialTemporalKEcons_SG/twice_VOFtransport";
foldername3="3SLmomentum_CG";
folderpick=foldername3;
string=folderpick+"/Run_La1.2e+2_nx";
ca_8_e2= readmatrix(string+"8/8.csv");
ca_16_e2= readmatrix(string+"16/16.csv");
ca_32_e2= readmatrix(string+"32/32.csv");
ca_64_e2= readmatrix(string+"64/64.csv");

string=folderpick+"/Run_La1.2e+4_nx";
ca_8_e4= readmatrix(string+"8/8.csv");
ca_16_e4= readmatrix(string+"16/16.csv");
ca_32_e4= readmatrix(string+"32/32.csv");
ca_64_e4= readmatrix(string+"64/64.csv");

string=folderpick+"/Run_La1.2e+6_nx";
ca_8_e6= readmatrix(string+"8/8.csv");
ca_16_e6= readmatrix(string+"16/16.csv");
ca_32_e6= readmatrix(string+"32/32.csv");
ca_64_e6= readmatrix(string+"64/64.csv");

blue=[71,135,224]/255;
red =[218,62,32]/255;
green= [82,147,47]/255;
yellow= [241,178,56]/255;
LW1=1.0;

LW2=2.5;

figure('Color','w')
cmax_8=max(ca_8_e2, [], 2);
cmax_16=max(ca_16_e2, [], 2);
cmax_32=max(ca_32_e2, [], 2);
cmax_64=max(ca_64_e2, [], 2);
if size(cmax_8)>2
    ca_rms_U_e2=[cmax_8(1) cmax_16(1) cmax_32(1) cmax_64(1)];
    ca_max_U_e2=[cmax_8(3) cmax_16(3) cmax_32(3) cmax_64(3)];
else
    ca_rms_U_e2=[cmax_8(1) cmax_16(1) cmax_32(1) cmax_64(1)];
    ca_max_U_e2=[cmax_8(2) cmax_16(2) cmax_32(2) cmax_64(2)];
end

cmax_8=max(ca_8_e4, [], 2);
cmax_16=max(ca_16_e4, [], 2);
cmax_32=max(ca_32_e4, [], 2);
cmax_64=max(ca_64_e4, [], 2);
if size(cmax_8)>2
    ca_rms_U_e4=[cmax_8(1) cmax_16(1) cmax_32(1) cmax_64(1)];
    ca_max_U_e4=[cmax_8(3) cmax_16(3) cmax_32(3) cmax_64(3)];
else
    ca_rms_U_e4=[cmax_8(1) cmax_16(1) cmax_32(1) cmax_64(1)];
    ca_max_U_e4=[cmax_8(2) cmax_16(2) cmax_32(2) cmax_64(2)];
end


cmax_8=max(ca_8_e6, [], 2);
cmax_16=max(ca_16_e6, [], 2);
cmax_32=max(ca_32_e6, [], 2);
cmax_64=max(ca_64_e6, [], 2);

if size(cmax_8)>2
    ca_rms_U_e6=[cmax_8(1) cmax_16(1) cmax_32(1) cmax_64(1)];
    ca_max_U_e6=[cmax_8(3) cmax_16(3) cmax_32(3) cmax_64(3)];
else
    ca_rms_U_e6=[cmax_8(1) cmax_16(1) cmax_32(1) cmax_64(1)];
    ca_max_U_e6=[cmax_8(2) cmax_16(2) cmax_32(2) cmax_64(2)];
end

cdx_D = sqrt(2)*[1/8 1/16 1/32 1/64]/0.4;

MS=10;
sy=3*10^(-3);
sy2=2*10^(-5);
cdx_DD= [0.5 0.25 0.125 0.125/2 0.125/4];
firstorder = sy*[1 1/2 1/4  1/8 1/16];
secondorder = sy2*[1 1/4 1/16 1/64 1/256];
loglog(cdx_D,ca_rms_U_e2,'LineWidth',LW2,'LineStyle','--','Color',blue,'Marker',"square",'MarkerSize',MS,'MarkerFaceColor',blue)
hold on
loglog(cdx_D,ca_max_U_e2,'LineWidth',LW2,'LineStyle','-','Color',blue,'Marker','o','MarkerSize',MS,'MarkerFaceColor',blue)
hold on
loglog(cdx_D,ca_rms_U_e4,'LineWidth',LW2,'LineStyle','--','Color',red,'Marker',"square",'MarkerSize',MS,'MarkerFaceColor',red)
hold on
loglog(cdx_D,ca_max_U_e4,'LineWidth',LW2,'LineStyle','-','Color',red,'Marker','o','MarkerSize',MS,'MarkerFaceColor',red)
hold on
loglog(cdx_D,ca_rms_U_e6,'LineWidth',LW2,'LineStyle','--','Color',green,'Marker',"square",'MarkerSize',MS,'MarkerFaceColor',green)
hold on
loglog(cdx_D,ca_max_U_e6,'LineWidth',LW2,'LineStyle','-','Color',green,'Marker','o','MarkerSize',MS,'MarkerFaceColor',green)
hold on
loglog(cdx_DD,firstorder,'LineWidth',3,'LineStyle',':','Color','k')
hold on
loglog(cdx_DD,secondorder,'LineWidth',3,'LineStyle','--','Color','k')

xlim([0.05,0.5])
xticks([0.05,0.1,0.5])
% xticklabels([10^(-2)])
ylim([10^(-7),10^(-2)])
% hold on
% loglog(cdx_D,ca_rms_Uhat_e2,'LineWidth',LW2,'LineStyle','--')
% hold on
% loglog(cdx_D,ca_max_Uhat_e2,'LineWidth',LW2,'LineStyle','--')
% 
% xlim([0.01,1])
% ylim([10^(-9),10^(-2)])
% ylim([-0.01,0.01])
% xlabel("$C\Delta = \sqrt{2}\Delta/D$", 'Interpreter', 'latex')
% ylabel("max(Ca)", 'Interpreter', 'latex')
% legend({"Ca RMS (U) La=1.2e2","Ca Max(U) La=1.2e2","Ca RMS (U) La=1.2e4","Ca Max(U) La=1.2e4","Ca RMS (U) La=1.2e6","Ca Max(U) La=1.2e6"}, 'Location', 'southeast', 'Interpreter', 'latex');
set(gca,'Fontsize',20)
set(gca,'fontname','Times New Roman')
set(gca,'LineWidth',2.0)

disp('Ca rms e2')
disp(ca_rms_U_e2)
disp('Ca max e2')
disp(ca_max_U_e2)
disp('Ca rms e4')
disp(ca_rms_U_e4)
disp('Ca max e4')
disp(ca_max_U_e4)
disp('Ca rms e6')
disp(ca_rms_U_e6)
disp('Ca max e6')
disp(ca_max_U_e6)