%% getting from KE conserving case
clear

filename = "results/simulation_DE_32";
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

filename = "results/simulation_DE_64";
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

% % 
filename = "results/simulation_DE_128";
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


%% StandingWave Amplitude

LW1=1.0;
LW2=2.0;
t=0:0.1:20; 
b =0.0136958;
omega= 0.995237;
Aex= 2.00+0.02*exp(-b*t).*cos(omega*t);


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
fitted_params_128 = lsqcurvefit(model_fun, initial_guess, time_128, xmax_128, lb, ub, opts);


fprintf('frequency error %.15f\n', abs(fitted_params_32(4)/omega-1))
fprintf('frequency error %.15f\n', abs(fitted_params_64(4)/omega-1))
fprintf('frequency error %.15f\n', abs(fitted_params_128(4)/omega-1))