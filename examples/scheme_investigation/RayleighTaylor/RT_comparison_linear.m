%% Rayleigh–Taylor Instability: Scheme Comparison vs Viscous Theory
%  Reads monitor/simulation and output.csv from each scheme.
%  Compares amplitude, KE, and PE growth against exact viscous growth rates.
%  Generates separate plots for each gas viscosity case.
%
%  Usage: Run from the RayleighTaylor/ directory.
clear; clc; close all;

%% -------- Physical parameters --------
rho_l = 1.225;       % heavy (liquid, VF=1)
rho_g = 0.1694;      % light (gas,    VF=0)
g_val = 9.81;        % gravity magnitude
A0    = 0.012;       % perturbation amplitude
Lx    = 1.0;         % domain width

At = (rho_l - rho_g) / (rho_l + rho_g);   % Atwood number
k  = 2*pi / Lx;                           % wavenumber
n_inviscid = sqrt(k * g_val * At);

fprintf('Atwood number At   = %.4f\n', At);
fprintf('Wavenumber k       = %.4f (2*pi)\n', k);
fprintf('Inviscid  n        = %.4f s^-1\n\n', n_inviscid);

%% -------- Case & Scheme Definitions --------
% Case 1: High gas viscosity
cases(1).dir      = 'arith_stats_256_same';
cases(1).mu_g     = 0.00313;
cases(1).n_theory = 6.34882;

% Case 2: Low gas viscosity
cases(2).dir      = 'arith_stats_256';
cases(2).mu_g     = 0.000432834;
cases(2).n_theory = 6.73591;

scheme_labels = {
    'Default NGA2 (SG)', 
    'KE Conserving (SG)', 
    'SL Momentum (CG)', 
    'SL Momentum (CG) sharper', 
    'SL Momentum (CG) 2 sub', 
    'SL Momentum (CG) sharper 2 sub'
    % , 'KE Conserving (SG) uhat visc' 
};
n_schemes = numel(scheme_labels);

% Define 6 distinct colors for the schemes
blue   = [71, 135, 224] / 255;
red    = [218, 62, 32] / 255;
green  = [82, 147, 47] / 255;
purple = [119, 43, 160] / 255;
orange = [230, 137, 23] / 255;
teal   = [29, 177, 186] / 255;
pink   = [219, 76, 152] / 255;
colors = [blue; red; green; purple; orange; teal; pink];

LW1 = 1.5;   % theory line width
LW2 = 1.7;   % simulation line width
FS  = 15;    % font size

%% -------- Linear regime bounds --------
fprintf('--- Linear Regime Limits ---\n');
for c = 1:length(cases)
    cases(c).t_linear_end = acosh(1/(k*A0)) / cases(c).n_theory;
    fprintf('Case %d (mu_g=%.4g): valid until t ~ %.3f s\n', c, cases(c).mu_g, cases(c).t_linear_end);
end
fprintf('\n');

%% -------- Read files --------
time_sim = cell(length(cases), n_schemes);
amp_sim  = cell(length(cases), n_schemes);
time_csv = cell(length(cases), n_schemes);
KE_tot   = cell(length(cases), n_schemes);
PE_drop  = cell(length(cases), n_schemes);

for c = 1:length(cases)
    for s = 1:n_schemes
        base_dir = sprintf('result_%s/%d', cases(c).dir, s);
        
        % 1. Read Amplitude
        mfile = fullfile(base_dir, 'simulation');
        if isfile(mfile)
            raw_sim = readmatrix(mfile, 'FileType', 'text', 'NumHeaderLines', 2);
            time_sim{c, s} = raw_sim(2:end, 2);    
            amp_sim{c, s}  = raw_sim(2:end, 17);   
        else
            warning('Monitor file not found: %s', mfile);
        end
        
        % 2. Read Energies
        csvfile = fullfile(base_dir, 'output.csv');
        if isfile(csvfile)
            raw_csv = readmatrix(csvfile, 'FileType', 'text', 'NumHeaderLines', 1);
            time_csv{c, s} = raw_csv(:, 1);
            KE_tot{c, s}   = raw_csv(:, 3);
            
            PE_l = raw_csv(:, 15);
            PE_g = raw_csv(:, 16);
            PE_drop{c, s} = PE_l(1) + PE_g(1) - PE_l - PE_g;
        else
            warning('output.csv not found: %s', csvfile);
        end
    end
end

t_th = linspace(0, 0.5, 2000);
t_ref = 0.3; % Anchor point for theory lines

%% -------- Generate Plots per Case --------
for c = 1:length(cases)
    
    title_suffix = sprintf('($\\mu_g=%g$)', cases(c).mu_g);
    
    % --- Figure A: Log Amplitude Growth ---
    figure('Color','w', 'Position', [100+c*50, 100+c*50, 800, 600])
    hold on
    for s = 1:n_schemes
        if isempty(amp_sim{c,s}); continue; end
        plot(time_sim{c,s}, log(amp_sim{c,s}), '-', 'Color', colors(s,:), ...
             'LineWidth', LW2, 'DisplayName', scheme_labels{s})
    end
    
    if ~isempty(amp_sim{c,1})
        [~, idx] = min(abs(time_sim{c,1} - t_ref));
        logA_ref = log(amp_sim{c,1}(idx));
        theory_line = logA_ref + cases(c).n_theory * (t_th - t_ref);
        plot(t_th, theory_line, 'k--', 'LineWidth', LW1, ...
             'DisplayName', sprintf('Theory $n=%.3f$', cases(c).n_theory))
    end
    xline(cases(c).t_linear_end, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2, 'HandleVisibility', 'off')
    
    xlabel('$t$ (s)', 'Interpreter', 'latex')
    ylabel('$\ln(A)$', 'Interpreter', 'latex')
    title(sprintf('RT Instability: Log Amplitude Growth %s', title_suffix), 'Interpreter', 'latex')
    legend('Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 10, 'NumColumns', 2)
    set(gca, 'FontSize', FS, 'FontName', 'Times New Roman', 'LineWidth', 1.7); grid on; xlim([0,0.55]);

    % --- Figure B: Log Total Kinetic Energy ---
    figure('Color','w', 'Position', [120+c*50, 80+c*50, 800, 600])
    hold on
    for s = 1:n_schemes
        if isempty(KE_tot{c,s}); continue; end
        plot(time_csv{c,s}, log(KE_tot{c,s}), '-', 'Color', colors(s,:), ...
             'LineWidth', LW2, 'DisplayName', scheme_labels{s})
    end
    
    if ~isempty(KE_tot{c,1})
        [~, idx] = min(abs(time_csv{c,1} - t_ref));
        logKE_ref = log(KE_tot{c,1}(idx));
        theory_line_KE = logKE_ref + 2*cases(c).n_theory * (t_th - t_ref);
        plot(t_th, theory_line_KE, 'k--', 'LineWidth', LW1, ...
             'DisplayName', sprintf('Theory $2n=%.3f$', 2*cases(c).n_theory))
    end
    xline(cases(c).t_linear_end, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2, 'HandleVisibility', 'off')
    
    xlabel('$t$ (s)', 'Interpreter', 'latex')
    ylabel('$\ln(\mathrm{KE})$', 'Interpreter', 'latex')
    title(sprintf('RT Instability: Log Total Kinetic Energy %s', title_suffix), 'Interpreter', 'latex')
    legend('Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 10, 'NumColumns', 2)
    set(gca, 'FontSize', FS, 'FontName', 'Times New Roman', 'LineWidth', 1.7); grid on; xlim([0, 0.55]);

    % --- Figure C: Log Potential Energy Drop ---
    figure('Color','w', 'Position', [140+c*50, 60+c*50, 800, 600])
    hold on
    for s = 1:n_schemes
        if isempty(PE_drop{c,s}); continue; end
        plot(time_csv{c,s}, log(PE_drop{c,s}), '-', 'Color', colors(s,:), ...
             'LineWidth', LW2, 'DisplayName', scheme_labels{s})
    end
    
    if ~isempty(PE_drop{c,1})
        [~, idx] = min(abs(time_csv{c,1} - t_ref));
        logPE_ref = log(PE_drop{c,1}(idx));
        theory_line_PE = logPE_ref + 2*cases(c).n_theory * (t_th - t_ref);
        plot(t_th, theory_line_PE, 'k--', 'LineWidth', LW1, ...
             'DisplayName', sprintf('Theory $2n=%.3f$', 2*cases(c).n_theory))
    end
    xline(cases(c).t_linear_end, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2, 'HandleVisibility', 'off')
    
    xlabel('$t$ (s)', 'Interpreter', 'latex')
    ylabel('$\ln(\Delta\mathrm{PE})$', 'Interpreter', 'latex')
    title(sprintf('RT Instability: Log Potential Energy Drop %s', title_suffix), 'Interpreter', 'latex')
    legend('Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 10, 'NumColumns', 2)
    set(gca, 'FontSize', FS, 'FontName', 'Times New Roman', 'LineWidth', 1.7); grid on; xlim([0, 0.55]);
end

%% -------- Growth Rate Fitting Summary --------
t_fit_start = 0.20;   
t_fit_end   = 0.45;

fprintf('\n===== AMPLITUDE Growth Rate Fitting (%.2f - %.2f s) =====\n', t_fit_start, t_fit_end);
fprintf('%-38s | %-8s | %-8s | %-8s\n', 'Scheme & Case', 'n_fit', 'n_theory', 'Err (%)');
fprintf('----------------------------------------------------------------------\n');
for c = 1:length(cases)
    for s = 1:n_schemes
        if isempty(amp_sim{c,s}); continue; end
        mask = (time_sim{c,s} >= t_fit_start) & (time_sim{c,s} <= t_fit_end);
        if sum(mask) > 5
            p = polyfit(time_sim{c,s}(mask), log(amp_sim{c,s}(mask)), 1);
            err_pct = abs(p(1) - cases(c).n_theory) / cases(c).n_theory * 100;
            label_str = sprintf('%s (mu_g=%.4g)', scheme_labels{s}, cases(c).mu_g);
            fprintf('%-38s | %.4f   | %.4f   | %.2f%%\n', label_str, p(1), cases(c).n_theory, err_pct);
        end
    end
end

fprintf('\n===== KINETIC ENERGY Growth Rate Fitting (%.2f - %.2f s) =====\n', t_fit_start, t_fit_end);
fprintf('%-38s | %-8s | %-8s | %-8s\n', 'Scheme & Case', 'n_fit', 'n_theory', 'Err (%)');
fprintf('----------------------------------------------------------------------\n');
for c = 1:length(cases)
    n_E_theory = 2 * cases(c).n_theory;
    for s = 1:n_schemes
        if isempty(KE_tot{c,s}); continue; end
        mask = (time_csv{c,s} >= t_fit_start) & (time_csv{c,s} <= t_fit_end);
        if sum(mask) > 5
            p = polyfit(time_csv{c,s}(mask), log(KE_tot{c,s}(mask)), 1);
            err_pct = abs(p(1) - n_E_theory) / n_E_theory * 100;
            label_str = sprintf('%s (mu_g=%.4g)', scheme_labels{s}, cases(c).mu_g);
            fprintf('%-38s | %.4f   | %.4f   | %.2f%%\n', label_str, p(1), n_E_theory, err_pct);
        end
    end
end

fprintf('\n===== POTENTIAL ENERGY Growth Rate Fitting (%.2f - %.2f s) =====\n', t_fit_start, t_fit_end);
fprintf('%-38s | %-8s | %-8s | %-8s\n', 'Scheme & Case', 'n_fit', 'n_theory', 'Err (%)');
fprintf('----------------------------------------------------------------------\n');
for c = 1:length(cases)
    n_E_theory = 2 * cases(c).n_theory;
    for s = 1:n_schemes
        if isempty(PE_drop{c,s}); continue; end
        mask = (time_csv{c,s} >= t_fit_start) & (time_csv{c,s} <= t_fit_end);
        if sum(mask) > 5
            p = polyfit(time_csv{c,s}(mask), log(PE_drop{c,s}(mask)), 1);
            err_pct = abs(p(1) - n_E_theory) / n_E_theory * 100;
            label_str = sprintf('%s (mu_g=%.4g)', scheme_labels{s}, cases(c).mu_g);
            fprintf('%-38s | %.4f   | %.4f   | %.2f%%\n', label_str, p(1), n_E_theory, err_pct);
        end
    end
end
fprintf('======================================================================\n');