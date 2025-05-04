%% getting from KE conserving case
clear
% close all

filename = "results/simulation_KE_32";
startRow = 3;
formatSpec = '%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
simulation = [dataArray{1:end-1}];
clearvars filename startRow formatSpec fileID dataArray ans;
start_index=1;
time_32=simulation(start_index:end-1,2);
xmax_32=abs(simulation(start_index+1:end,15));

filename = "results/simulation_KE_64";
startRow = 3;
formatSpec = '%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
simulation = [dataArray{1:end-1}];
clearvars filename startRow formatSpec fileID dataArray ans;
start_index=1;
time_64=simulation(start_index:end-1,2);
xmax_64=abs(simulation(start_index+1:end,15));
% 
filename = "results/simulation_KE_96";
startRow = 3;
formatSpec = '%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
simulation = [dataArray{1:end-1}];
clearvars filename startRow formatSpec fileID dataArray ans;
start_index=1;
time_96=simulation(start_index:end-1,2);
xmax_96=abs(simulation(start_index+1:end,15));
% % 
filename = "results/simulation_KE_128";
startRow = 3;
formatSpec = '%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
simulation = [dataArray{1:end-1}];
clearvars filename startRow formatSpec fileID dataArray ans;
start_index=1;
time_128=simulation(start_index:end-1,2);
xmax_128=abs(simulation(start_index+1:end,15));
% 
% 
filename = "results/simulation_KE_160";
startRow = 3;
formatSpec = '%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
simulation = [dataArray{1:end-1}];
clearvars filename startRow formatSpec fileID dataArray ans;
start_index=1;
time_160=simulation(start_index:end-1,2);
xmax_160=abs(simulation(start_index+1:end,15));



%% StandingWave Amplitude

LW1=1.0;
LW2=2.0;
t=0:0.1:20; 
b =0.0136958;
omega= 0.995237;
Aex= 2.00+0.02*exp(-b*t).*cos(omega*t);



% % % perod_32 =diff([0 (6.118+6.156)/2.0 12.084 18.164])/(2*pi);
% perod_48 =diff([0 6.08 (12.274+12.312)/2.0 (18.468+18.506)/2.0])/(2*pi);
% perod_64 =diff([0 (6.042+6.08)/2.0 (12.274+12.236)/2.0 (18.43+18.468)/2.0])/(2*pi);
% perod_96 =diff([0 6.286 12.572 18.858])/(2*pi);
% perod_128 =diff([0 (6.258+6.3)/2.0 (12.544+12.586)/2.0 (18.844+18.858)/2.0])/(2*pi);
% perod_160 =diff([0 (6.31+6.32)/2.0 (12.6+12.67)/2.0 (18.92+18.98)/2.0])/(2*pi);
% 
% % fprintf('%.15f\n', abs(mean(perod_32)-omega))
% fprintf('%.15f\n', abs(mean(perod_48)-omega))
% fprintf('%.15f\n', abs(mean(perod_64)-omega))
% fprintf('%.15f\n', abs(mean(perod_96)-omega))
% fprintf('%.15f\n', abs(mean(perod_128)-omega))
% fprintf('%.15f\n', abs(mean(perod_160)-omega))

figure('Color','w')
plot(time_32,xmax_32,'LineWidth',LW2,'LineStyle','-')
hold on
% plot(time_48,xmax_48,'LineWidth',LW2,'LineStyle','-')
% hold on
plot(time_64,xmax_64,'LineWidth',LW2,'LineStyle','-')
hold on
plot(time_96,xmax_96,'LineWidth',LW2,'LineStyle','-')
hold on
plot(time_128,xmax_128,'LineWidth',LW2,'LineStyle','-')
hold on
plot(time_160,xmax_160,'LineWidth',LW2,'LineStyle','-')
hold on
plot(t,Aex,'LineWidth',LW1,'Marker','*')
xlim([0,20])
ylim([1.965,2.02])
% ylim([-0.01,0.01])
legend({"64","96","128","160","Prosperetti"}, 'Location', 'southeast', 'Interpreter', 'latex');
set(gca,'Fontsize',15)
set(gca,'fontname','Times New Roman')
set(gca,'LineWidth',1.7)

%%
% % error_32=abs(mean(perod_32)-omega);
% % error_48=abs(mean(perod_48)-omega);
% error_64=abs(mean(perod_64)-omega);
% error_96=abs(mean(perod_96)-omega);
% error_128=abs(mean(perod_128)-omega);
% error_160=abs(mean(perod_160)-omega);
% % dx=[20/32 20/64 20/96 20/128 20/160];
% % err_list=[error_32 error_64 error_96 error_128 error_160];
% dx=[ 20/64 20/96 20/128 20/160];
% err_list=[ error_64 error_96 error_128 error_160];
% 
% x_1st=[20/64 20/(64*1.5)];
% y_1st=[0.02  0.02/1.5];
% x_2nd=[20/64 20/(64*1.5)];
% y_2nd=[0.01  0.01/(1.5^2)];
% 
% 
% figure
% loglog(dx, err_list,'LineWidth',LW2)
% hold on
% loglog(x_1st, y_1st,'LineWidth',LW2)
% hold on 
% loglog(x_2nd, y_2nd,'LineWidth',LW2)
% xlim([0.1 0.4])
% ylim([0.003,0.03])
% set(gca,'Fontsize',15)
% set(gca,'fontname','Times New Roman')
% set(gca,'LineWidth',1.7)


%% Oscillation Error
% arith_Uhat_64=diff([0 (6.15+6.18)/2.0 12.33 18.52])/(2*pi);
% arith_Uhat_96=diff([0 6.272 12.558 (18.802+18.827)/2.0])/(2*pi);
% arith_Uhat_128=diff([0 6.28 (12.53+12.58)/2.0 18.84])/(2*pi);
% arith_Uhat_160=diff([0 6.307 12.614 (18.907+18.935)/2.0])/(2*pi);
% fprintf('%.15f\n', abs(mean(arith_Uhat_64)-omega))
% fprintf('%.15f\n', abs(mean(arith_Uhat_96)-omega))
% fprintf('%.15f\n', abs(mean(arith_Uhat_128)-omega))
% fprintf('%.15f\n', abs(mean(arith_Uhat_160)-omega))

% harm_Uhat_64=diff([0 (6.12+6.15)/2.0 (12.27+12.3)/2.0 (18.42+18.45)/2.0])/(2*pi);
% harm_Uhat_96=diff([0 (6.216+6.23)/2.0 (12.432+12.46)/2.0 (18.676+18.69)/2.0])/(2*pi);
% harm_Uhat_128=diff([0 6.23 (12.44+12.49)/2.0 (18.68+18.72)/2.0])/(2*pi);
% harm_Uhat_160=diff([0 (6.237+6.293)/2.0 (12.509+12.572)/2.0 (18.774+18.837)/2.0])/(2*pi);
% fprintf('%.15f\n', abs(mean(harm_Uhat_64)-omega))
% fprintf('%.15f\n', abs(mean(harm_Uhat_96)-omega))
% fprintf('%.15f\n', abs(mean(harm_Uhat_128)-omega))
% fprintf('%.15f\n', abs(mean(harm_Uhat_160)-omega))

% %% harmonic vs arithmetic viscousity
% mul = 0.01;
% mug = 0.0001;
% alpha=0:0.01:1;
% mu1 = (mul *mug)./(mug*alpha + mul*(1 -alpha));
% mu2 = mul*alpha + mug*(1 - alpha);
% 
% figure
% plot(alpha,mu1/mul,'LineWidth',LW2)
% hold on
% plot(alpha,mu2/mul,'LineWidth',LW2)
% xlabel("$\alpha$", 'Interpreter', 'latex')
% ylabel("$\mu/\mu_l$", 'Interpreter', 'latex')
% legend({"Harmonic","Arithmetic"}, 'Location', 'northwest', 'Interpreter', 'latex');
% set(gca,'Fontsize',15)
% set(gca,'fontname','Times New Roman')
% set(gca,'LineWidth',1.7)

%% A better way to get the frequency


model_fun = @(params, x) params(1) + params(2) * exp(-params(3) * x) .* cos(params(4)*x);
initial_guess=[2.0,0.02,b,omega];
opts = optimset('Display', 'iter'); % Show iterations
lb = [-Inf, -Inf, 0, 0]; % Lower bounds
ub = [Inf, Inf, Inf, Inf]; % Upper bounds
fitted_params_32= lsqcurvefit(model_fun, initial_guess, time_32, xmax_32, lb, ub, opts);



model_fun = @(params, x) params(1) + params(2) * exp(-params(3) * x) .* cos(params(4)*x);
initial_guess=[2.0,0.02,b,omega];
opts = optimset('Display', 'iter'); % Show iterations
lb = [-Inf, -Inf, 0, 0]; % Lower bounds
ub = [Inf, Inf, Inf, Inf]; % Upper bounds
fitted_params_64 = lsqcurvefit(model_fun, initial_guess, time_64, xmax_64, lb, ub, opts);


model_fun = @(params, x) params(1) + params(2) * exp(-params(3) * x) .* cos(params(4)*x);
initial_guess=[2.0,0.02,b,omega];
opts = optimset('Display', 'iter'); % Show iterations
lb = [-Inf, -Inf, 0, 0]; % Lower bounds
ub = [Inf, Inf, Inf, Inf]; % Upper bounds
fitted_params_96 = lsqcurvefit(model_fun, initial_guess, time_96, xmax_96, lb, ub, opts);


model_fun = @(params, x) params(1) + params(2) * exp(-params(3) * x) .* cos(params(4)*x);
initial_guess=[2.0,0.02,b,omega];
opts = optimset('Display', 'iter'); % Show iterations
lb = [-Inf, -Inf, 0, 0]; % Lower bounds
ub = [Inf, Inf, Inf, Inf]; % Upper bounds
fitted_params_128 = lsqcurvefit(model_fun, initial_guess, time_128, xmax_128, lb, ub, opts);


model_fun = @(params, x) params(1) + params(2) * exp(-params(3) * x) .* cos(params(4)*x);
initial_guess=[2.0,0.02,b,omega];
opts = optimset('Display', 'iter'); % Show iterations
lb = [-Inf, -Inf, 0, 0]; % Lower bounds
ub = [Inf, Inf, Inf, Inf]; % Upper bounds
fitted_params_160 = lsqcurvefit(model_fun, initial_guess, time_160, xmax_160, lb, ub, opts);


fprintf('frequency error %.15f\n', abs(fitted_params_32(4)/omega-1))
fprintf('frequency error %.15f\n', abs(fitted_params_64(4)/omega-1))
fprintf('frequency error %.15f\n', abs(fitted_params_96(4)/omega-1))
fprintf('frequency error %.15f\n', abs(fitted_params_128(4)/omega-1))
fprintf('frequency error %.15f\n', abs(fitted_params_160(4)/omega-1))


% fprintf('frequency error %.15f\n', abs(fitted_params_64(4)-omega))
% fprintf('frequency error %.15f\n', abs(fitted_params_96(4)-omega))
% fprintf('frequency error %.15f\n', abs(fitted_params_128(4)-omega))
% fprintf('frequency error %.15f\n', abs(fitted_params_160(4)-omega))

% 
% fprintf('damping error %.15f\n', abs(fitted_params_64(3)-b)/b)
% fprintf('damping error %.15f\n', abs(fitted_params_96(3)-b)/b)
% fprintf('damping error %.15f\n', abs(fitted_params_128(3)-b)/b)
% fprintf('damping error %.15f\n', abs(fitted_params_160(3)-b)/b)


%Plot just in case I am curious how well the fitting is
% figure;
% plot(time_64, xmax_64, 'bo', 'MarkerSize', 5); % Data
% hold on;
% plot(time_64, model_fun(fitted_params, time_64), 'r-', 'LineWidth', 2); % Fitted curve
% legend('Data', 'Fitted Curve');
% title('Curve Fitting: a + b * exp(-c * x) * cos(d * x)');
% grid on;
