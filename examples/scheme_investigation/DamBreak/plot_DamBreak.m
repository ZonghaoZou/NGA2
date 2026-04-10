%% plot_DamBreak.m — Compare dam break front position across 3 NGA2 schemes
%  vs. Martin & Moyce (1952) experimental data
%
%  Usage: Run from the DamBreak/ directory
%         >> plot_DamBreak

clear; close all; clc;

%% ===== User Configuration ===================================================
base_dir = '.';  % Run from scheme_investigation/DamBreak/
scheme_dirs = {'1DefaultNGA2_SG', '2SpatialTemporalKEcons_SG', '3SLmomentum_CG'};
scheme_labels = {'Scheme 1: Default NGA2 (SG)', ...
                 'Scheme 2: KE-Conserving (SG)', ...
                 'Scheme 3: SL-Momentum (CG)'};
scheme_colors = {[0.0 0.45 0.74], [0.85 0.33 0.10], [0.47 0.67 0.19]};
resolution = '320';
resolutions_compare = {'160', '320'};  % resolutions for Plot 3

%% ===== Martin & Moyce (1952) Experimental Data ==============================
% Table 2: Rectangular section, n^2 = 2, a = 2 1/4 in = 0.05715 m
% Z = z/a (front position normalized by column width, starts at 1)
% T = n*t*sqrt(g/a) where n = sqrt(n^2) = sqrt(2)
% Mean column from Table 2 for a = 2 1/4 in (digitized from paper)
Z_exp=[1.11,1.22,1.44,1.67,1.89,2.11,2.33,2.56,2.78,3.00,3.22,3.44,3.67,3.89,4.11,5.00,5.44,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0];
T_exp=[0.41,0.84,1.19,1.43,1.63,1.83,1.98,2.20,2.32,2.51,2.65,2.83,2.97,3.11,3.33,4.02,4.44,5.09,5.69,6.30,6.83,7.44,8.08,8.67,9.31];

%% ===== Column layout (18 columns) ===========================================
col_time    = 2;
col_frontx  = 15;
col_columnh = 16;
col_Znd     = 17;
col_Tnd     = 18;
ncols_expected = 18;

%% ===== Read simulation data for primary resolution ===========================
sim_data = read_monitor_data(base_dir, scheme_dirs, resolution, ncols_expected, ...
                             col_time, col_frontx, col_columnh, col_Znd, col_Tnd);

%% ===== Plot 1: Non-dimensional front position Z vs T ========================
figure('Position', [100 100 800 600]);

plot(T_exp, Z_exp, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', [0.5 0.5 0.5], ...
     'DisplayName', 'Martin & Moyce (1952)');
hold on;

for s = 1:length(scheme_dirs)
    if ~isempty(sim_data{s})
        plot(sim_data{s}.T_nd, sim_data{s}.Z_nd, '-', ...
             'Color', scheme_colors{s}, 'LineWidth', 2.0, ...
             'DisplayName', scheme_labels{s});
    end
end

xlabel('T = n \cdot t \cdot \surd(g/a),  n^2 = 2', 'FontSize', 14);
ylabel('Z = z / a', 'FontSize', 14);
title('Dam Break: Surge Front Position (Martin & Moyce 1952, Table 2)', 'FontSize', 16);
legend('Location', 'northwest', 'FontSize', 11);
grid on;
set(gca, 'FontSize', 12);
xlim([0 10]);
ylim([0 14]);

saveas(gcf, 'front_position_ZvsT.png');
fprintf('Saved: front_position_ZvsT.png\n');

%% ===== Plot 2: Time-shifted Z vs T ==========================================
Z_match = 1.44;
T_exp_match = interp1(Z_exp, T_exp, Z_match, 'linear');

T_shift = 0;
for s = 1:length(scheme_dirs)
    if ~isempty(sim_data{s})
        [Z_unique, idx] = unique(sim_data{s}.Z_nd);
        T_unique = sim_data{s}.T_nd(idx);
        if Z_match >= min(Z_unique) && Z_match <= max(Z_unique)
            T_sim_match = interp1(Z_unique, T_unique, Z_match, 'linear');
            T_shift = T_exp_match - T_sim_match;
            fprintf('Time shift: T_shift = %.4f (matched at Z=%.2f)\n', T_shift, Z_match);
            break;
        end
    end
end

figure('Position', [200 100 800 600]);

plot(T_exp, Z_exp, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', [0.5 0.5 0.5], ...
     'DisplayName', 'Martin & Moyce (1952)');
hold on;

for s = 1:length(scheme_dirs)
    if ~isempty(sim_data{s})
        plot(sim_data{s}.T_nd + T_shift, sim_data{s}.Z_nd, '-', ...
             'Color', scheme_colors{s}, 'LineWidth', 2.0, ...
             'DisplayName', scheme_labels{s});
    end
end

xlabel('T = n \cdot t \cdot \surd(g/a),  n^2 = 2', 'FontSize', 14);
ylabel('Z = z / a', 'FontSize', 14);
title(sprintf('Dam Break: Front Position (shifted \\DeltaT = %.2f at Z=%.2f)', T_shift, Z_match), 'FontSize', 16);
legend('Location', 'northwest', 'FontSize', 11);
grid on;
set(gca, 'FontSize', 12);
xlim([0 10]);
ylim([0 14]);

saveas(gcf, 'front_position_ZvsT_shifted.png');
fprintf('Saved: front_position_ZvsT_shifted.png\n');

%% ===== Plot 3: Resolution comparison (160 vs 320) ===========================
% Load data for both resolutions and plot all 3 schemes x 2 resolutions
res_linestyles = {'-', '--'};  % solid for first res, dashed for second

figure('Position', [300 100 900 650]);

% Plot experimental data
plot(T_exp, Z_exp, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', [0.5 0.5 0.5], ...
     'DisplayName', 'Martin & Moyce (1952)');
hold on;

for r = 1:length(resolutions_compare)
    res = resolutions_compare{r};
    res_data = read_monitor_data(base_dir, scheme_dirs, res, ncols_expected, ...
                                 col_time, col_frontx, col_columnh, col_Znd, col_Tnd);
    
    for s = 1:length(scheme_dirs)
        if ~isempty(res_data{s})
            plot(res_data{s}.T_nd, res_data{s}.Z_nd, res_linestyles{r}, ...
                 'Color', scheme_colors{s}, 'LineWidth', 2.0, ...
                 'DisplayName', sprintf('%s (nx=%s)', scheme_labels{s}, res));
        end
    end
end

xlabel('T = n \cdot t \cdot \surd(g/a),  n^2 = 2', 'FontSize', 14);
ylabel('Z = z / a', 'FontSize', 14);
title('Dam Break: Resolution Comparison (solid=160, dashed=320)', 'FontSize', 16);
legend('Location', 'northwest', 'FontSize', 10);
grid on;
set(gca, 'FontSize', 12);
xlim([0 10]);
ylim([0 14]);

saveas(gcf, 'front_position_resolution_comparison.png');
fprintf('Saved: front_position_resolution_comparison.png\n');

fprintf('\nDone. All plots saved.\n');


%% ===== Helper function: read monitor data ====================================
function sim_data = read_monitor_data(base_dir, scheme_dirs, res, ncols_expected, ...
                                      col_time, col_frontx, col_columnh, col_Znd, col_Tnd)
    sim_data = cell(length(scheme_dirs), 1);
    for s = 1:length(scheme_dirs)
        monfile = fullfile(base_dir, scheme_dirs{s}, 'monitor', strcat('simulation_', res));
        if ~exist(monfile, 'file')
            fprintf('  Warning: %s not found\n', monfile);
            sim_data{s} = [];
            continue;
        end
        
        fid = fopen(monfile, 'r');
        raw = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
        fclose(fid);
        lines = raw{1};
        
        data = [];
        for li = 3:length(lines)
            line = strtrim(lines{li});
            if isempty(line); continue; end
            line = regexprep(line, '(\d)(-\d{2,3})\b', '$1E$2');
            line = regexprep(line, '(\d)(\+\d{2,3})\b', '$1E$2');
            vals = str2num(line);  %#ok<ST2NM>
            if length(vals) == ncols_expected
                data = [data; vals];  %#ok<AGROW>
            end
        end
        
        if isempty(data)
            sim_data{s} = [];
            continue;
        end
        
        fprintf('  Read %d timesteps from %s (res=%s)\n', size(data,1), scheme_dirs{s}, res);
        s_data = struct();
        s_data.time     = data(:, col_time);
        s_data.front_x  = data(:, col_frontx);
        s_data.column_h = data(:, col_columnh);
        s_data.Z_nd     = data(:, col_Znd);
        s_data.T_nd     = data(:, col_Tnd);
        sim_data{s} = s_data;
    end
end
