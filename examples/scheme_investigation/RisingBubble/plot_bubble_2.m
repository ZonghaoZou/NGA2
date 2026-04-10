% plot_bubble.m — Rising Bubble benchmark analysis (Hysing et al. 2009)
% Plots rise velocity, center of mass, and circularity vs. time
% Overlays Featflow reference data (Group 1, level 7)
% Style matches spurious current analysis script
%
% Expected data layout:
%   RisingBubble/{scheme}/{nx}/benchmark.csv
%   Reference:  RisingBubble/reference/c1g1l7.txt

clear; close all; clc;

% =========================================================================
%  Style settings (matching spurious current studies)
% =========================================================================
colors = {[0 0.4470 0.7410], ...     % blue  — Scheme 1
          [0.8500 0.3250 0.0980], ... % red   — Scheme 2
          [0.4660 0.6740 0.1880]};    % green — Scheme 3

line_styles = {'-', '--', '-.'};
marker_types = {'s', 'o', 'd'};  % square, circle, diamond
line_width = 2.5;
marker_size = 10;
font_size = 20;
font_name = 'Times New Roman';

casetype="result/case2/";
scheme_names = {'Default NGA2 (SG)', 'KE Cons. (SG)', 'SL Mom. (CG)'};
scheme_dirs  = {strcat(casetype,'1'), strcat(casetype,'2'), strcat(casetype,'3')};
nx_values    = [40, 80, 160, 320];

% =========================================================================
%  Load reference data (Featflow: Group 1, Level 7)
% =========================================================================
ref_file = fullfile('reference', 'c2g1l8.txt');
if exist(ref_file, 'file')
    ref_data = load(ref_file);
    ref_t    = ref_data(:,1);
    ref_area = ref_data(:,2);
    ref_circ = ref_data(:,3);
    ref_CoM  = ref_data(:,4);
    ref_Vr   = ref_data(:,5);
    has_ref  = true;
    fprintf('Loaded reference data: %d points\n', length(ref_t));
else
    has_ref = false;
    fprintf('WARNING: Reference file not found: %s\n', ref_file);
end

% =========================================================================
%  Load simulation data
% =========================================================================
data = struct();
for s = 1:length(scheme_dirs)
    for n = 1:length(nx_values)
        nx = nx_values(n);
        csv_file = fullfile(scheme_dirs{s}, sprintf('benchmark_%d.csv', nx));
        if exist(csv_file, 'file')
            raw = load(csv_file);
            data(s,n).t    = raw(:,1);
            data(s,n).area = raw(:,2);
            data(s,n).circ = raw(:,3);
            data(s,n).CoM  = raw(:,4);
            data(s,n).Vr   = raw(:,5);
            data(s,n).valid = true;
            fprintf('Loaded: %s (nx=%d, %d points)\n', scheme_dirs{s}, nx, length(raw(:,1)));
        else
            data(s,n).valid = false;
            fprintf('Missing: %s\n', csv_file);
        end
    end
end

% =========================================================================
%  Figure 1: Rise Velocity vs. Time
% =========================================================================
figure('Position', [100, 100, 800, 600]);
hold on; box on; grid on;

if has_ref
    plot(ref_t, ref_Vr, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Hysing (TP2D, L7)');
end

% Plot finest resolution for each scheme
for s = 1:length(scheme_dirs)
    % Find finest available resolution
    for n = length(nx_values):-1:1
        if isfield(data, 'valid') || true
            if data(s,n).valid
                plot(data(s,n).t, data(s,n).Vr, ...
                    'Color', colors{s}, 'LineStyle', line_styles{s}, ...
                    'LineWidth', line_width, ...
                    'DisplayName', sprintf('%s (nx=%d)', scheme_names{s}, nx_values(n)));
                break;
            end
        end
    end
end

xlabel('Time', 'FontSize', font_size, 'FontName', font_name);
ylabel('Rise velocity', 'FontSize', font_size, 'FontName', font_name);
set(gca, 'FontSize', font_size, 'FontName', font_name);
legend('Location', 'northeast', 'FontSize', 14);
title('Rising Bubble — Rise Velocity', 'FontSize', font_size, 'FontName', font_name);
xlim([0 3]);

% =========================================================================
%  Figure 2: Center of Mass vs. Time
% =========================================================================
figure('Position', [100, 100, 800, 600]);
hold on; box on; grid on;

if has_ref
    plot(ref_t, ref_CoM, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Hysing (TP2D, L7)');
end

for s = 1:length(scheme_dirs)
    for n = length(nx_values):-1:1
        if data(s,n).valid
            plot(data(s,n).t, data(s,n).CoM, ...
                'Color', colors{s}, 'LineStyle', line_styles{s}, ...
                'LineWidth', line_width, ...
                'DisplayName', sprintf('%s (nx=%d)', scheme_names{s}, nx_values(n)));
            break;
        end
    end
end

xlabel('Time', 'FontSize', font_size, 'FontName', font_name);
ylabel('Center of mass (y)', 'FontSize', font_size, 'FontName', font_name);
set(gca, 'FontSize', font_size, 'FontName', font_name);
legend('Location', 'northwest', 'FontSize', 14);
title('Rising Bubble — Center of Mass', 'FontSize', font_size, 'FontName', font_name);
xlim([0 3]);

% =========================================================================
%  Figure 3: Circularity vs. Time
% =========================================================================
figure('Position', [100, 100, 800, 600]);
hold on; box on; grid on;

if has_ref
    plot(ref_t, ref_circ, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Hysing (TP2D, L7)');
end

for s = 1:length(scheme_dirs)
    for n = length(nx_values):-1:1
        if data(s,n).valid
            plot(data(s,n).t, data(s,n).circ, ...
                'Color', colors{s}, 'LineStyle', line_styles{s}, ...
                'LineWidth', line_width, ...
                'DisplayName', sprintf('%s (nx=%d)', scheme_names{s}, nx_values(n)));
            break;
        end
    end
end

xlabel('Time', 'FontSize', font_size, 'FontName', font_name);
ylabel('Circularity', 'FontSize', font_size, 'FontName', font_name);
set(gca, 'FontSize', font_size, 'FontName', font_name);
legend('Location', 'southwest', 'FontSize', 14);
title('Rising Bubble — Circularity', 'FontSize', font_size, 'FontName', font_name);
xlim([0 3]);
ylim([0.5 1.01]);

% =========================================================================
%  Figure 4: Convergence — Plot all resolutions for one scheme
% =========================================================================
for s = 1:length(scheme_dirs)
    figure('Position', [100, 100, 800, 600]);
    hold on; box on; grid on;
    
    if has_ref
        plot(ref_t, ref_Vr, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Hysing (TP2D, L7)');
    end
    
    for n = 1:length(nx_values)
        if data(s,n).valid
            plot(data(s,n).t, data(s,n).Vr, ...
                'Color', colors{s}, 'LineStyle', line_styles{min(n,3)}, ...
                'LineWidth', line_width - 0.5*(length(nx_values)-n), ...
                'DisplayName', sprintf('nx=%d', nx_values(n)));
        end
    end
    
    xlabel('Time', 'FontSize', font_size, 'FontName', font_name);
    ylabel('Rise velocity', 'FontSize', font_size, 'FontName', font_name);
    set(gca, 'FontSize', font_size, 'FontName', font_name);
    legend('Location', 'northeast', 'FontSize', 14);
    title(sprintf('Rise Vel. Convergence — %s', scheme_names{s}), 'FontSize', font_size, 'FontName', font_name);
    xlim([0 3]);
end

% =========================================================================
%  Print summary at t=3
% =========================================================================
fprintf('\n=== Benchmark values at t=3.0 ===\n');
fprintf('%-30s  %10s  %10s  %10s  %10s\n', 'Config', 'Area', 'Circ', 'CoM_y', 'V_rise');
fprintf('%-30s  %10s  %10s  %10s  %10s\n', '------', '----', '----', '-----', '------');

if has_ref
    idx = find(abs(ref_t - 3.0) < 0.01, 1, 'last');
    if ~isempty(idx)
        fprintf('%-30s  %10.4f  %10.4f  %10.4f  %10.4f\n', ...
            'Hysing TP2D L7', ref_area(idx), ref_circ(idx), ref_CoM(idx), ref_Vr(idx));
    end
end

for s = 1:length(scheme_dirs)
    for n = 1:length(nx_values)
        if data(s,n).valid
            idx = find(abs(data(s,n).t - 3.0) < 0.05, 1, 'last');
            if ~isempty(idx)
                fprintf('%-30s  %10.4f  %10.4f  %10.4f  %10.4f\n', ...
                    sprintf('%s nx=%d', scheme_names{s}, nx_values(n)), ...
                    data(s,n).area(idx), data(s,n).circ(idx), ...
                    data(s,n).CoM(idx), data(s,n).Vr(idx));
            end
        end
    end
end
