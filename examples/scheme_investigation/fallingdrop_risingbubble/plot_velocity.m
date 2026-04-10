% plot_velocity.m
%
% Script to plot falling droplet velocity for 3 schemes vs theoretical solution
%
clear; clc; close all;

T_max=10;
% Cases definition
testcase=1;
if testcase==1
    foldername='result';
else
    
    foldername='result_p25_case3';
end 
blue   = [71, 135, 224] / 255;
red    = [218, 62, 32] / 255;
green  = [82, 147, 47] / 255;
purple = [119, 43, 160] / 255;
orange = [230, 137, 23] / 255;
teal   = [29, 177, 186] / 255;
pink   = [219, 76, 152] / 255;
cases  = {strcat(foldername,'/1'), strcat(foldername,'/2'), strcat(foldername,'/3'), strcat(foldername,'/4'), strcat(foldername,'/5'), strcat(foldername,'/6'), strcat(foldername,'/7')};
names  = {'Default', 'KE con', 'SL (full band)', 'SL (fewer face)', 'SL (fewer face) 2 subitr', 'SL (fewer face) 2 subitr widedomain', 'Droplet'};
colors = {blue, ...     % blue  — Scheme 1
          red, ... % red   — Scheme 2
          green, ...
          purple, ...
          orange, ...
          teal, ...
          pink};    % green — Scheme 3
% colors = {'r', 'g', 'b'};     % Colors for lines
styles = {'-', '-', '-','-','-','-','-'};     % Line styles

% Prepare Figure
figure('Name', 'Falling Droplet Velocity', 'Color', 'w', 'Position', [100, 100, 800, 600]);
hold on; box on; grid on;

% --- Loop over schemes ---
for i = 1:length(cases)
    filepath = fullfile(cases{i}, 'bubble');
    
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
                vrise= data(:, 4);
                vinflow = data(:, 5);
                
                % Plot simulation data
                plot(time, vinflow, 'Color', colors{i}, 'LineStyle', styles{i}, ...
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
% 
if testcase==1
    t_theory = linspace(0, T_max, 100);
    v_theory1 = ones(100,1)*0.72;    
    plot(t_theory, v_theory1, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Experiment');
    hold on
    v_theory1 = ones(100,1)*0.7185; 
    plot(t_theory, v_theory1, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Massey (3D sphere)');
    hold on
    v_theory2 = ones(100,1)*0.7064; 
   
    plot(t_theory, v_theory2, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Schiller and Naumann (3D sphere)');
    hold on
    v_theory2 = ones(100,1)*0.947589; 
    
    plot(t_theory, v_theory2, 'g--', 'LineWidth', 1.5, 'DisplayName', 'White (2D cylinder)');
    hold on
    v_theory2 = ones(100,1)*1.0417; 
    
    plot(t_theory, v_theory2, 'm--', 'LineWidth', 1.5, 'DisplayName', 'Sucker and Brauer (2D cylinder)');

else
    t_theory = linspace(0, T_max, 100);
    % v_theory1 = ones(100,1)*0.72; 
    v_theory1 = ones(100,1)*4.03; 
    plot(t_theory, v_theory1, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Experiment');
    % hold on
    % v_theory1 = ones(100,1)*0.7185; 
    % plot(t_theory, v_theory1, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Massey (3D sphere)');
    hold on
    % v_theory2 = ones(100,1)*0.7064; 
    v_theory2 = ones(100,1)*3.864; 
    plot(t_theory, v_theory2, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Schiller and Naumann (3D sphere)');
    hold on
    % v_theory2 = ones(100,1)*0.947589; 
    v_theory2 = ones(100,1)*3.1588; 
    plot(t_theory, v_theory2, 'g--', 'LineWidth', 1.5, 'DisplayName', 'White (2D cylinder)');
    hold on
    % v_theory2 = ones(100,1)*1.0417; 
    v_theory2 = ones(100,1)*3.148; 
    plot(t_theory, v_theory2, 'm--', 'LineWidth', 1.5, 'DisplayName', 'Sucker and Brauer (2D cylinder)');
end

% --- Formatting ---
xlabel('Time [s]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Abs Vertical Velocity [m/s]', 'FontSize', 12, 'FontWeight', 'bold');
title('Falling Droplet Velocity Comparison', 'FontSize', 14);
legend('Location', 'SouthWest', 'FontSize', 10);
set(gca, 'FontSize', 12);
if testcase==1
    ylim([0, 1.1]); % Adjust Y limit to show falling
else
    ylim([0, 5.0]); % Adjust Y limit to show falling
end
xlim([0, 0.6]);

hold off;
fprintf('Plotting complete.\n');
