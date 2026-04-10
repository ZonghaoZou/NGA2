% read_and_plot_kh.m
% Script to read and plot the Kelvin-Helmholtz growth rate data

% 1. Read the data from the CSV file
% readmatrix automatically detects and skips the text header row
close all; clc; clear;

blue   = [71, 135, 224] / 255;
red    = [218, 62, 32] / 255;
green  = [82, 147, 47] / 255;
purple = [119, 43, 160] / 255;
orange = [230, 137, 23] / 255;
teal   = [29, 177, 186] / 255;
pink   = [219, 76, 152] / 255;

% %% CASE A
% filename = 'sol/sol_casea.csv';
% data = readmatrix(filename);
% 
% filename = 'results/caseA/1/sweep_results.csv';
% data1 = readmatrix(filename);
% 
% filename = 'results/caseA/2/sweep_results.csv';
% data2 = readmatrix(filename);
% 
% filename = 'results/caseA/3/sweep_results.csv';
% data3 = readmatrix(filename);
% 
% % 2. Extract the columns
% alpha_nd = data(:, 1);
% growth_rate = data(:, 2);
% 
% alpha_nd1 = data1(:, 1);
% growth_rate1 = data1(:, 2);
% 
% alpha_nd2 = data2(:, 1);
% growth_rate2 = data2(:, 2);
% 
% alpha_nd3 = data3(:, 1);
% growth_rate3 = data3(:, 2);
% 
% % 3. Create the plot
% figure('Name', 'KH Dispersion Curve', 'Color', 'w');
% plot(alpha_nd, growth_rate, 'k-', 'LineWidth', 2.0);
% hold on
% scatter(alpha_nd1, growth_rate1, 80, 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', blue,'MarkerFaceColor',blue);
% hold on
% scatter(alpha_nd2, growth_rate2, 80, 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', red,'MarkerFaceColor',red);
% hold on
% scatter(alpha_nd3, growth_rate3, 80, 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', green,'MarkerFaceColor',green);
% % 4. Format the plot aesthetics
% grid on;
% box on;
% 
% % Set axis limits (starts at 0, gives 10% headroom on the peak)
% xlim([0, 2.0]);
% ylim([0, 0.25]);
% 
% % Add labels using the TeX interpreter for greek letters
% xlabel('$\alpha \delta_g$', 'FontSize', 14, 'Interpreter', 'latex');
% ylabel('$\alpha c_i \delta_g / U_{g,\infty}$', 'FontSize', 14, 'Interpreter', 'latex');
% title('Theoretical Kelvin-Helmholtz Growth Rate Case A', 'FontSize', 14);
% legend("Orr--Sommerfeld","Default NGA2", "KEcon", "SLmom", 'FontSize', 12, 'Location','best');
% % Increase tick label font size for readability
% set(gca, 'FontSize', 12);
% 
% 
% %% CASE B
% filename = 'sol/sol_caseb.csv';
% data = readmatrix(filename);
% 
% filename = 'results/caseB/1/sweep_results.csv';
% data1 = readmatrix(filename);
% 
% filename = 'results/caseB/2/sweep_results.csv';
% data2 = readmatrix(filename);
% 
% filename = 'results/caseB/3/sweep_results.csv';
% data3 = readmatrix(filename);
% 
% % 2. Extract the columns
% alpha_nd = data(:, 1);
% growth_rate = data(:, 2);
% 
% alpha_nd1 = data1(:, 1);
% growth_rate1 = data1(:, 2);
% 
% alpha_nd2 = data2(:, 1);
% growth_rate2 = data2(:, 2);
% 
% alpha_nd3 = data3(:, 1);
% growth_rate3 = data3(:, 2);
% 
% % 3. Create the plot
% figure('Name', 'KH Dispersion Curve', 'Color', 'w');
% plot(alpha_nd, growth_rate, 'k-', 'LineWidth', 2.0);
% hold on
% scatter(alpha_nd1, growth_rate1, 80, 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', blue,'MarkerFaceColor',blue);
% hold on
% scatter(alpha_nd2, growth_rate2, 80, 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', red,'MarkerFaceColor',red);
% hold on
% scatter(alpha_nd3, growth_rate3, 80, 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', green,'MarkerFaceColor',green);
% % 4. Format the plot aesthetics
% grid on;
% box on;
% 
% % Set axis limits (starts at 0, gives 10% headroom on the peak)
% xlim([0, 2.0]);
% ylim([0, 0.25]);
% 
% % Add labels using the TeX interpreter for greek letters
% xlabel('$\alpha \delta_g$', 'FontSize', 14, 'Interpreter', 'latex');
% ylabel('$\alpha c_i \delta_g / U_{g,\infty}$', 'FontSize', 14, 'Interpreter', 'latex');
% title('Theoretical Kelvin-Helmholtz Growth Rate Case B', 'FontSize', 14);
% legend("Orr--Sommerfeld","Default NGA2", "KEcon", "SLmom", 'FontSize', 12, 'Location','best');
% % Increase tick label font size for readability
% set(gca, 'FontSize', 12);


%% CASE C
filename = 'sol/sol_casec.csv';
data = readmatrix(filename);

filename = 'results/caseC/1/sweep_results.csv';
data1 = readmatrix(filename);

filename = 'results/caseC/2/sweep_results.csv';
data2 = readmatrix(filename);

filename = 'results/caseC/3/sweep_results.csv';
data3 = readmatrix(filename);

% 2. Extract the columns
alpha_nd = data(:, 1);
growth_rate = data(:, 2);

alpha_nd1 = data1(:, 1);
growth_rate1 = data1(:, 2);

alpha_nd2 = data2(:, 1);
growth_rate2 = data2(:, 2);

alpha_nd3 = data3(:, 1);
growth_rate3 = data3(:, 2);

% 3. Create the plot
figure('Name', 'KH Dispersion Curve', 'Color', 'w');
plot(alpha_nd, growth_rate, 'k-', 'LineWidth', 2.0);
hold on
scatter(alpha_nd1, growth_rate1, 80, 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', blue,'MarkerFaceColor',blue);
hold on
scatter(alpha_nd2, growth_rate2, 80, 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', red,'MarkerFaceColor',red);
hold on
scatter(alpha_nd3, growth_rate3, 80, 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', green,'MarkerFaceColor',green);
% 4. Format the plot aesthetics
grid on;
box on;

% Set axis limits (starts at 0, gives 10% headroom on the peak)
xlim([0, 4.0]);
ylim([0, 0.25]);

% Add labels using the TeX interpreter for greek letters
xlabel('$\alpha \delta_g$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\alpha c_i \delta_g / U_{g,\infty}$', 'FontSize', 14, 'Interpreter', 'latex');
title('Theoretical Kelvin-Helmholtz Growth Rate Case C', 'FontSize', 14);
legend("Orr--Sommerfeld","Default NGA2", "KEcon", "SLmom", 'FontSize', 12, 'Location','best');
% Increase tick label font size for readability
set(gca, 'FontSize', 12);


%% CASE D
filename = 'sol/sol_cased.csv';
data = readmatrix(filename);

filename = 'results/caseD/1/sweep_results.csv';
data1 = readmatrix(filename);

filename = 'results/caseD/2/sweep_results.csv';
data2 = readmatrix(filename);

filename = 'results/caseD/3/sweep_results.csv';
data3 = readmatrix(filename);

% 2. Extract the columns
alpha_nd = data(:, 1);
growth_rate = data(:, 2);

alpha_nd1 = data1(:, 1);
growth_rate1 = data1(:, 2);

alpha_nd2 = data2(:, 1);
growth_rate2 = data2(:, 2);

alpha_nd3 = data3(:, 1);
growth_rate3 = data3(:, 2);

% 3. Create the plot
figure('Name', 'KH Dispersion Curve', 'Color', 'w');
plot(alpha_nd, growth_rate, 'k-', 'LineWidth', 2.0);
hold on
scatter(alpha_nd1, growth_rate1, 80, 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', blue,'MarkerFaceColor',blue);
hold on
scatter(alpha_nd2, growth_rate2, 80, 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', red,'MarkerFaceColor',red);
hold on
scatter(alpha_nd3, growth_rate3, 80, 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', green,'MarkerFaceColor',green);
% 4. Format the plot aesthetics
grid on;
box on;

% Set axis limits (starts at 0, gives 10% headroom on the peak)
xlim([0, 4.0]);
ylim([0, 0.25]);

% Add labels using the TeX interpreter for greek letters
xlabel('$\alpha \delta_g$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\alpha c_i \delta_g / U_{g,\infty}$', 'FontSize', 14, 'Interpreter', 'latex');
title('Theoretical Kelvin-Helmholtz Growth Rate Case D', 'FontSize', 14);
legend("Orr--Sommerfeld","Default NGA2", "KEcon", "SLmom", 'FontSize', 12, 'Location','best');
% Increase tick label font size for readability
set(gca, 'FontSize', 12);