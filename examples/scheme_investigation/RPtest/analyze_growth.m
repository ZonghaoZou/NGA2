%% Rayleigh-Plateau Instability Analysis: Multi-Scheme, Multi-Resolution
%  Reads simulation monitor files for 3 schemes across 4 mesh resolutions
%  and compares growth rates against theory.
%
%  Data layout: result/{schemeID}/simulation_{N}
%  Follows the same pattern as OscillatingDroplet/plot_oscillation.m
clear; clc; close all;

%% --- 1. USER CONFIGURATION ---
% Base schemes
base_schemes = { ...
    1, 'Default NGA2 (SG)'; ...
    2, 'KE Conservative (SG)'; ...
    3, 'SL Momentum (CG)' };

% Resolutions to compare
resolutions = [32, 64, 128, 256];

% Column indices in monitor file
col_time = 2;   % Time
col_rmin = 15;  % Min Radius
col_rmax = 16;  % Max Radius

% Theoretical growth rate
sigma_theory = 0.33065;

% Linear fit region (time window for exponential growth)
t_start = 5.0;
t_end   = 12.0;

%% --- 2. DATA LOADING ---
% Dynamically build list of available data files
scheme_dirs   = {};
scheme_labels = {};
scheme_ids    = [];  % track which base scheme each entry belongs to
res_ids       = [];  % track which resolution each entry belongs to

for i = 1:size(base_schemes, 1)
    id   = base_schemes{i, 1};
    name = base_schemes{i, 2};
    
    for r = 1:length(resolutions)
        res = resolutions(r);
        fpath = fullfile('result', num2str(id), sprintf('simulation_%d', res));
        
        if isfile(fpath)
            scheme_dirs{end+1}   = fpath; %#ok<SAGROW>
            scheme_labels{end+1} = sprintf('%s N=%d', name, res); %#ok<SAGROW>
            scheme_ids(end+1)    = i; %#ok<SAGROW>
            res_ids(end+1)       = r; %#ok<SAGROW>
        end
    end
end

n_data = numel(scheme_dirs);

% Storage
time_all  = cell(n_data, 1);
rmin_all  = cell(n_data, 1);
rmax_all  = cell(n_data, 1);

for s = 1:n_data
    fpath = scheme_dirs{s};
    raw = readmatrix(fpath, 'FileType', 'text', 'NumHeaderLines', 2);
    
    time_all{s} = raw(2:end-1, col_time);
    rmin_all{s} = raw(3:end, col_rmin);
    rmax_all{s} = raw(3:end, col_rmax);
    
    fprintf('Loaded %d timesteps for: %s\n', numel(time_all{s}), scheme_labels{s});
end

%% --- 3. GROWTH RATE FITTING ---
fprintf('\n============================================================\n');
fprintf('%-30s | %-10s | %-10s | %-8s\n', 'Scheme', 'Sigma_sim', 'Sigma_th', 'Err (%)');
fprintf('------------------------------------------------------------\n');

sigma_fit_all = nan(n_data, 1);

for s = 1:n_data
    if isempty(rmax_all{s}); continue; end
    
    % Calculate log-amplitude: ln( (Rmax - Rmin)/2 )
    amplitude = (rmax_all{s} - rmin_all{s}) / 2.0;
    amplitude(amplitude <= 0) = NaN;
    log_amp = log(amplitude);
    
    % Linear fit in the exponential growth window
    mask = (time_all{s} >= t_start) & (time_all{s} <= t_end) & ~isnan(log_amp);
    t_fit = time_all{s}(mask);
    y_fit = log_amp(mask);
    
    if isempty(t_fit); continue; end
    
    p = polyfit(t_fit, y_fit, 1);
    sigma_fit_all(s) = p(1);
    
    err_pct = abs(sigma_fit_all(s) - sigma_theory) / sigma_theory * 100;
    fprintf('%-30s | %.5f    | %.5f    | %.2f%%\n', ...
        scheme_labels{s}, sigma_fit_all(s), sigma_theory, err_pct);
end
fprintf('============================================================\n');

%% --- 4. FIGURE 1: All schemes and resolutions (log-amplitude vs time) ---
LW1 = 1.5;
LW2 = 1.7;
FS  = 15;

res_colors = lines(length(resolutions));  % Consistent per-resolution colors
scheme_linestyles = {'-', '--', '-.'};     % Different line style per scheme

figure('Color','w', 'Position', [100 100 900 600])
hold on; grid on; box on;

for s = 1:n_data
    if isempty(rmax_all{s}); continue; end
    
    amplitude = (rmax_all{s} - rmin_all{s}) / 2.0;
    amplitude(amplitude <= 0) = NaN;
    log_amp = log(amplitude);
    
    lstyle = scheme_linestyles{scheme_ids(s)};
    lcolor = res_colors(res_ids(s), :);
    
    plot(time_all{s}, log_amp, lstyle, 'Color', lcolor, ...
         'LineWidth', LW2, 'DisplayName', scheme_labels{s})
end

% Add theory reference line
x_ref = [t_start, t_end];
y_mid = mean(ylim);
y_ref = sigma_theory .* (x_ref - mean(x_ref)) + y_mid;
plot(x_ref, y_ref, 'g:', 'LineWidth', 2.5, ...
     'DisplayName', sprintf('Theory ($\\sigma=%.4f$)', sigma_theory))

% Highlight fit region
try
    xregion(t_start, t_end, 'FaceColor', [0.9 0.9 0.9], 'FaceAlpha', 0.5, ...
            'DisplayName', 'Fit Window');
catch
    yl = ylim;
    patch([t_start t_end t_end t_start], [yl(1) yl(1) yl(2) yl(2)], ...
          'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

xlabel('$t$', 'Interpreter', 'latex')
ylabel('$\ln(\mathrm{Amplitude})$', 'Interpreter', 'latex')
title('Rayleigh-Plateau Growth Rate: All Schemes \& Resolutions', 'Interpreter', 'latex')
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 10)
set(gca, 'FontSize', FS, 'FontName', 'Times New Roman', 'LineWidth', 1.7)
xlim([0, 20])

%% --- 5. FIGURES 2-4: Per-scheme grid convergence ---
for i = 1:size(base_schemes, 1)
    current_scheme_name = base_schemes{i, 2};
    
    figure('Color','w', 'Position', [100 100 700 500], ...
           'Name', ['Grid Convergence: ' current_scheme_name])
    hold on; grid on; box on;
    
    for s = 1:n_data
        if scheme_ids(s) ~= i; continue; end
        if isempty(rmax_all{s}); continue; end
        
        amplitude = (rmax_all{s} - rmin_all{s}) / 2.0;
        amplitude(amplitude <= 0) = NaN;
        log_amp = log(amplitude);
        
        lcolor = res_colors(res_ids(s), :);
        
        % Plot data
        plot(time_all{s}, log_amp, '-', 'Color', lcolor, ...
             'LineWidth', LW2, 'DisplayName', scheme_labels{s})
        
        % Overlay fitted line
        mask = (time_all{s} >= t_start) & (time_all{s} <= t_end) & ~isnan(log_amp);
        if any(mask)
            p = polyfit(time_all{s}(mask), log_amp(mask), 1);
            plot(time_all{s}(mask), polyval(p, time_all{s}(mask)), ...
                 '--', 'Color', lcolor, 'LineWidth', 2.0, 'HandleVisibility', 'off')
        end
    end
    
    % Theory reference
    x_ref = [t_start, t_end];
    y_mid = mean(ylim);
    y_ref = sigma_theory .* (x_ref - mean(x_ref)) + y_mid;
    plot(x_ref, y_ref, 'g:', 'LineWidth', 2.5, ...
         'DisplayName', sprintf('Theory ($\\sigma=%.4f$)', sigma_theory))
    
    % Fit region shading
    try
        xregion(t_start, t_end, 'FaceColor', [0.9 0.9 0.9], 'FaceAlpha', 0.5, ...
                'DisplayName', 'Fit Window');
    catch
        yl = ylim;
        patch([t_start t_end t_end t_start], [yl(1) yl(1) yl(2) yl(2)], ...
              'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    
    xlabel('$t$', 'Interpreter', 'latex')
    ylabel('$\ln(\mathrm{Amplitude})$', 'Interpreter', 'latex')
    title([current_scheme_name ': Grid Convergence'], 'Interpreter', 'latex')
    legend('Location', 'best', 'Interpreter', 'latex')
    set(gca, 'FontSize', FS, 'FontName', 'Times New Roman', 'LineWidth', 1.7)
    xlim([0, 20])
end

%% --- 6. FIGURE 5: Convergence plot (error vs Dx) ---
figure('Color','w', 'Position', [100 100 600 450])
hold on; grid on; box on;

scheme_markers = {'o-', 's-', 'd-'};  % circle, square, diamond
scheme_colors  = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0 0 0];

% Domain length
L = 10.0;

for i = 1:size(base_schemes, 1)
    current_scheme_name = base_schemes{i, 2};
    
    Dx_vals  = [];
    err_vals = [];
    
    for s = 1:n_data
        if scheme_ids(s) ~= i; continue; end
        if isnan(sigma_fit_all(s)); continue; end
        
        N = resolutions(res_ids(s));
        Dx_vals(end+1)  = L / N; %#ok<SAGROW>
        err_vals(end+1) = abs(sigma_fit_all(s) - sigma_theory) / sigma_theory * 100; %#ok<SAGROW>
    end
    
    if ~isempty(Dx_vals)
        loglog(Dx_vals, err_vals, scheme_markers{i}, ...
               'Color', scheme_colors(i,:), 'MarkerFaceColor', scheme_colors(i,:), ...
               'LineWidth', LW2, 'MarkerSize', 8, 'DisplayName', current_scheme_name)
    end
end

% Add reference convergence slopes anchored to coarsest data point
Dx_ref = [L/max(resolutions), L/min(resolutions)];
Dx_coarse = L/min(resolutions);
% Find error at coarsest resolution from first available scheme
err_anchor = nan;
for s = 1:n_data
    if res_ids(s) == 1 && ~isnan(sigma_fit_all(s))
        err_anchor = abs(sigma_fit_all(s) - sigma_theory) / sigma_theory * 100;
        break;
    end
end
if isnan(err_anchor); err_anchor = 50; end
% 1st order reference: error ~ Dx
err_ref1 = err_anchor * (Dx_ref / Dx_coarse);
loglog(Dx_ref, err_ref1, 'k--', 'LineWidth', 1.0, 'DisplayName', '$O(\Delta x)$')
% 2nd order reference: error ~ Dx^2
err_ref2 = err_anchor * (Dx_ref / Dx_coarse).^2;
loglog(Dx_ref, err_ref2, 'k:', 'LineWidth', 1.0, 'DisplayName', '$O(\Delta x^2)$')

xlabel('$\Delta x$', 'Interpreter', 'latex')
ylabel('Growth rate error (\%)', 'Interpreter', 'latex')
title('Grid Convergence of Growth Rate', 'Interpreter', 'latex')
legend('Location', 'best', 'Interpreter', 'latex')
set(gca, 'FontSize', FS, 'FontName', 'Times New Roman', 'LineWidth', 1.7)