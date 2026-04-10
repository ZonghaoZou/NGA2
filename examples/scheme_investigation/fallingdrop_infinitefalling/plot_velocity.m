% plot_velocity.m
%
% Script to plot falling droplet velocity for 3 schemes vs theoretical solution
%
clear; clc; close all;

% --- Parameters ---
g       = 1.0;         % Gravity (physical)
rho_l   = 1000.0;       % Liquid density
rho_g   = 1.0;          % Gas density
T_max   = 15.0;          % Max simulation time to plot




% Theoretical acceleration for inviscid sphere (with added mass)
% Equation: (rho_l + Cm * rho_g) * a = (rho_l - rho_g) * g
% Cm = 0.5 for sphere
a_added_mass = g * (rho_l - rho_g) / (rho_l + 1.0* rho_g);

% Cases definition
cases  = {'1', '2', '3'};
names  = {'Scheme 1', 'Scheme 2', 'Scheme 3'};
colors = {[0 0.4470 0.7410], ...     % blue  — Scheme 1
          [0.8500 0.3250 0.0980], ... % red   — Scheme 2
          [0.4660 0.6740 0.1880]};    % green — Scheme 3
% colors = {'r', 'g', 'b'};     % Colors for lines
styles = {'-', '-', '-'};     % Line styles

% Prepare Figure
figure('Name', 'Falling Droplet Velocity', 'Color', 'w', 'Position', [100, 100, 800, 600]);
hold on; box on; grid on;

% --- Loop over schemes ---
for i = 1:length(cases)
    filepath = fullfile(cases{i}, 'monitor', 'simulation');
    
    if exist(filepath, 'file')
        % Import data (assumes NGA monitor format with headers)
        try
            raw = importdata(filepath);
            if isstruct(raw) && isfield(raw, 'data')
                data = raw.data;
            else
                data = raw;
            end
            
            % Check if data is valid
            if ~isempty(data)
                % Column 2 is Time
                time = data(:, 2);
                
                % Last Column is V_drop (Assumed from implementation plan)
                v_drop = data(:, end);
                
                % Plot simulation data
                plot(time, v_drop, 'Color', colors{i}, 'LineStyle', styles{i}, ...
                     'LineWidth', 2, 'DisplayName', names{i});
            else
                fprintf('Warning: Data file empty for %s\n', names{i});
            end
        catch ME
            fprintf('Error reading file %s: %s\n', filepath, ME.message);
        end
    else
        fprintf('Warning: File not found: %s\n', filepath);
    end
end

% --- Plot Theory ---
t_theory = linspace(0, T_max, 100);
v_theory = -a_added_mass * t_theory; % Negative because falling down
plot(t_theory, v_theory, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Theory (Added Mass)');

% --- Formatting ---
xlabel('Time [s]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Vertical Velocity [m/s]', 'FontSize', 12, 'FontWeight', 'bold');
title('Falling Droplet Velocity Comparison', 'FontSize', 14);
legend('Location', 'SouthWest', 'FontSize', 10);
set(gca, 'FontSize', 12);
ylim([min(v_theory)*1.1, 0]); % Adjust Y limit to show falling
xlim([0, T_max]);

hold off;
fprintf('Plotting complete.\n');
