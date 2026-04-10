%% Rayleigh–Taylor Instability: Scheme Comparison vs Viscous Theory
%  Reads monitor/simulation from each scheme, compares amplitude growth
%  against the exact viscous growth rate from Chandrasekhar Eq. 113.
%
%  Theory: A(t) = A0 * exp(n * t)  in the LINEAR regime
%  n = 7.38089 s^-1 (from solve.py, Chandrasekhar Eq. 113)
%
%  Valid comparison range: A(t) << lambda = 1.0
%    With A0=0.01, linear regime holds while A < ~0.1
%    => t < ln(10)/7.38 ~ 0.31 s
%
%  Usage: Run from the RayleighTaylor/ directory.
clear; clc; close all;

%% --------
% ';lkhgfdsa Physical parameters --------
rho_l = 1.225;       % heavy (liquid, VF=1)
rho_g = 0.1694;      % light (gas,    VF=0)
mu_l    = 0.00313;
mu_g    = 0.000432834;% dynamic viscosity (both phases)
g_val = 9.81;        % gravity magnitude
A0    = 0.012;        % perturbation amplitude
Lx    = 1.0;         % domain width
sigma = 0.0;         % surface tension
deltax=1.0/256.0;
% viscoutype="arith_stats_512";
% viscoutype = "stats_512_inviscid";
% viscoutype = "arith_stats_256_sharp";
viscoutype = "arith_stats_256_sharp_same";

% viscoutype="harm";
At = (rho_l - rho_g) / (rho_l + rho_g);   % Atwood number
k  = 2*pi / Lx;                            % wavenumber
nubar=(mu_l+mu_g)/(rho_l+rho_g);
% Hardcoded exact viscous growth rate (from solve.py, Chandrasekhar Eq 113)
n_viscous = 6.34882;
% n_viscous=6.73591;
% n_viscous=sqrt(k * g_val * At);
% n_viscous=-nubar*k^2 + sqrt(k * g_val * At +nubar^2*k^4 -sigma/(rho_l+rho_g)*k^3);
% Inviscid growth rate (for reference)
n_inviscid = sqrt(k * g_val * At);

fprintf('Atwood number At   = %.4f\n', At);
fprintf('Wavenumber k       = %.4f (2*pi)\n', k);
fprintf('Inviscid  n        = %.4f s^-1\n', n_inviscid);
fprintf('Viscous   n (Eq113)= %.5f s^-1\n', n_viscous);

%% -------- Linear regime bounds --------
% Linear theory valid while A(t) << wavelength lambda = 1
% Conservative bound: A(t) < 0.1 => t < ln(0.1/A0)/n
% t_linear_end = log(0.1 / A0) / n_viscous;
% t_linear_end = log(1/(k*A0)) / n_viscous;
t_linear_end = acosh(1/(k*A0)) / n_viscous;
fprintf('Linear regime valid until t ~ %.3f s\n', t_linear_end);
fprintf('  (where amplitude reaches ~0.1, i.e. 10%% of wavelength)\n\n');

%% -------- Formatting --------
LW1 = 1.5;   % theory line width
LW2 = 1.7;   % simulation line width
FS  = 15;    % font size

%% -------- Scheme definitions --------
scheme_dirs   = {strcat('result_',viscoutype,'/1'), strcat('result_',viscoutype,'/2'), strcat('result_',viscoutype,'/3')};
scheme_labels = {'Default NGA2 (SG)', 'KE Conservative (SG)', 'SL Momentum (CG)'};
n_schemes     = numel(scheme_dirs);
blue=[71,135,224]/255;
red =[218,62,32]/255;
green= [82,147,47]/255;

%% -------- Read monitor files --------
time_all  = cell(n_schemes, 1);
spike_all = cell(n_schemes, 1);
bub_all   = cell(n_schemes, 1);
amp_all   = cell(n_schemes, 1);

for s = 1:n_schemes
    mfile = fullfile(scheme_dirs{s}, 'simulation');
    if ~isfile(mfile)
        warning('Monitor file not found: %s', mfile);
        continue
    end
    
    raw = readmatrix(mfile, 'FileType', 'text', 'NumHeaderLines', 2);
    
    time_all{s}  = raw(2:end, 2);    % Time
    spike_all{s} = raw(2:end, 15);   % Spike Y
    bub_all{s}   = raw(2:end, 16);   % Bubble Y
    amp_all{s}   = raw(2:end, 17);   % Amplitude
    
    fprintf('Loaded %d timesteps for: %s\n', numel(time_all{s}), scheme_labels{s});
end

%% -------- Theory curve --------
t_th = linspace(0, 0.5, 2000);
A_theory = A0 * cosh(n_viscous * t_th);

%% -------- Figure 1: Log amplitude — per-scheme shifted theory --------
colors = vertcat(vertcat(blue,red),green);
theory_color = [0.4 0.4 0.4];  % dark gray for theory lines

figure('Color','w', 'Position', [100 100 700 500])
hold on

for s = 1:n_schemes
    if isempty(amp_all{s}); continue; end
    
    % Compute ln(A) from simulation
    log_amp = log(amp_all{s});
    
    % Plot simulation
    plot(time_all{s}, log_amp, '-', 'Color', colors(s,:), ...
         'LineWidth', LW2, 'DisplayName', scheme_labels{s})
end

% Draw per-scheme theory lines anchored at t_ref=0.2
t_ref = 0.3;
for s = 1:n_schemes
    if isempty(amp_all{s}); continue; end
    [~, idx] = min(abs(time_all{s} - t_ref));
    logA_ref = log(amp_all{s}(idx));
    theory_line = logA_ref + n_viscous * (t_th - t_ref);
    hv = 'off'; if s == 1; hv = 'on'; end
    plot(t_th, theory_line, 'k--', 'LineWidth', LW1, ...
         'HandleVisibility', hv, ...
         'DisplayName', sprintf('Theory $n=%.3f$ s$^{-1}$', n_viscous))
end

% Mark end of linear regime
xline(t_linear_end, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
      'DisplayName', 'Linear regime limit')
hold on
plot(t_th,log(A_theory),'LineWidth',2.0,'Color','k','LineStyle',':','DisplayName', '$A_0 \cosh(\omega_{\text{theory}} t)$')

xlabel('$t$ (s)', 'Interpreter', 'latex')
ylabel('$\ln(A)$', 'Interpreter', 'latex')
title('RT Instability: Log Amplitude Growth', 'Interpreter', 'latex')
legend('Location', 'northwest', 'Interpreter', 'latex')
set(gca, 'FontSize', FS, 'FontName', 'Times New Roman', 'LineWidth', 1.7)
grid on
xlim([0,0.55])
% xlim([0 min(0.5, max(cellfun(@max, time_all(~cellfun(@isempty, time_all)))))])

%% -------- Growth rate fitting --------
fprintf('\n===== Growth Rate Fitting (linear regime only) =====\n');
fprintf('%-25s | %-10s | %-10s | %-8s\n', 'Scheme', 'n_fit', 'n_theory', 'Err (%)');
fprintf('---------------------------------------------------\n');

t_fit_start = 0.2;   % skip initial transient
t_fit_end   = 0.45;
% t_fit_end   = t_linear_end;

for s = 1:n_schemes
    if isempty(amp_all{s}); continue; end
    
    log_amp = log(amp_all{s});
    mask = (time_all{s} >= t_fit_start) & (time_all{s} <= t_fit_end) & ~isnan(log_amp) & ~isinf(log_amp);
    
    if sum(mask) < 5; continue; end
    
    p = polyfit(time_all{s}(mask), log_amp(mask), 1);
    n_fit = p(1);
    err_pct = abs(n_fit - n_viscous) / n_viscous * 100;
    
    fprintf('%-25s | %.4f    | %.4f    | %.2f%%\n', ...
        scheme_labels{s}, n_fit, n_viscous, err_pct);
end
fprintf('===================================================\n');

% %% ======================================================================
% %  PART 2: Physical Quantity Evolution from output.csv
% %  Follows the pattern of Phaseinversion/phaseinversion.m
% %% ======================================================================
% 
% %% -------- Read output.csv files --------
% % output.csv columns (up to 24):
% %  1: time       2: timestep
% %  3: KE         4: KE_l       5: KE_g
% %  6: rhoU       7: rhoV       8: rhoW
% %  9: rhoU_l    10: rhoV_l    11: rhoW_l
% % 12: rhoU_g    13: rhoV_g    14: rhoW_g
% % 15: PE_l      16: PE_g
% % 17: EN_l      18: EN_g
% % 19: KE_x      20: KE_xl     21: KE_xg
% % 22: KE_y      23: KE_yl     24: KE_yg
% 
% csv_time   = cell(n_schemes, 1);
% csv_KE     = cell(n_schemes, 1);
% csv_KE_l   = cell(n_schemes, 1);
% csv_KE_g   = cell(n_schemes, 1);
% csv_rhoU   = cell(n_schemes, 1);
% csv_rhoV   = cell(n_schemes, 1);
% csv_rhoW   = cell(n_schemes, 1);
% csv_rhoU_l = cell(n_schemes, 1);
% csv_rhoV_l = cell(n_schemes, 1);
% csv_rhoW_l = cell(n_schemes, 1);
% csv_rhoU_g = cell(n_schemes, 1);
% csv_rhoV_g = cell(n_schemes, 1);
% csv_rhoW_g = cell(n_schemes, 1);
% csv_PE_l   = cell(n_schemes, 1);
% csv_PE_g   = cell(n_schemes, 1);
% csv_EN_l   = cell(n_schemes, 1);
% csv_EN_g   = cell(n_schemes, 1);
% csv_KE_x   = cell(n_schemes, 1);
% csv_KE_xl  = cell(n_schemes, 1);
% csv_KE_xg  = cell(n_schemes, 1);
% csv_KE_y   = cell(n_schemes, 1);
% csv_KE_yl  = cell(n_schemes, 1);
% csv_KE_yg  = cell(n_schemes, 1);
% 
% for s = 1:n_schemes
%     csvfile = fullfile(scheme_dirs{s}, 'output.csv');
%     if ~isfile(csvfile)
%         warning('output.csv not found: %s', csvfile);
%         continue
%     end
% 
%     raw = readmatrix(csvfile, 'FileType', 'text', 'NumHeaderLines', 1);
%     ncols = size(raw, 2);
% 
%     csv_time{s}   = raw(:, 1);
%     csv_KE{s}     = raw(:, 3);
%     csv_KE_l{s}   = raw(:, 4);
%     csv_KE_g{s}   = raw(:, 5);
%     csv_rhoU{s}   = raw(:, 6);
%     csv_rhoV{s}   = raw(:, 7);
%     csv_rhoW{s}   = raw(:, 8);
%     csv_rhoU_l{s} = raw(:, 9);
%     csv_rhoV_l{s} = raw(:, 10);
%     csv_rhoW_l{s} = raw(:, 11);
%     csv_rhoU_g{s} = raw(:, 12);
%     csv_rhoV_g{s} = raw(:, 13);
%     csv_rhoW_g{s} = raw(:, 14);
%     csv_PE_l{s}   = raw(:, 15);
%     csv_PE_g{s}   = raw(:, 16);
% 
%     if ncols >= 18
%         csv_EN_l{s} = raw(:, 17);
%         csv_EN_g{s} = raw(:, 18);
%     end
% 
%     if ncols >= 24
%         csv_KE_x{s}  = raw(:, 19);
%         csv_KE_xl{s} = raw(:, 20);
%         csv_KE_xg{s} = raw(:, 21);
%         csv_KE_y{s}  = raw(:, 22);
%         csv_KE_yl{s} = raw(:, 23);
%         csv_KE_yg{s} = raw(:, 24);
%     end
% 
%     fprintf('Loaded %d timesteps from output.csv for: %s (%d columns)\n', ...
%         size(raw,1), scheme_labels{s}, ncols);
% end
% 
% % Determine common x-axis limit from output.csv data
% t_max_csv = max(cellfun(@(x) max(x, [], 'omitnan'), ...
%     csv_time(~cellfun(@isempty, csv_time))));
% t_xlim = min(0.5, t_max_csv);
% 
% %% -------- Figure 4: Log Total Kinetic Energy --------
% figure('Color','w', 'Position', [100 100 700 500])
% hold on
% 
% for s = 1:n_schemes
%     if isempty(csv_KE{s}); continue; end
%     plot(csv_time{s}, log(csv_KE{s}), '-', 'Color', colors(s,:), ...
%          'LineWidth', LW2, 'DisplayName', scheme_labels{s})
% end
% 
% % Per-scheme straight-line theory with slope 2*n_viscous anchored at t_ref=0.2
% % (KE ~ exp(2*n*t) in the linear regime)
% t_ref = 0.3;
% for s = 1:n_schemes
%     if isempty(csv_KE{s}); continue; end
%     [~, idx_s] = min(abs(csv_time{s} - t_ref));
%     logKE_ref = log(csv_KE{s}(idx_s));
%     theory_line = logKE_ref + 2*n_viscous * (t_th - t_ref);
%     hv = 'off'; if s == 1; hv = 'on'; end
%     plot(t_th, theory_line, 'k--', 'LineWidth', LW1, ...
%          'HandleVisibility', hv, ...
%          'DisplayName', sprintf('Theory $2n=%.3f$ s$^{-1}$', 2*n_viscous))
% end
% 
% xline(t_linear_end, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
%       'DisplayName', 'Linear regime limit')
% 
% xlabel('$t$ (s)', 'Interpreter', 'latex')
% ylabel('$\ln(\mathrm{KE})$', 'Interpreter', 'latex')
% title('RT Instability: Log Total Kinetic Energy', 'Interpreter', 'latex')
% legend('Location', 'northwest', 'Interpreter', 'latex')
% set(gca, 'FontSize', FS, 'FontName', 'Times New Roman', 'LineWidth', 1.7)
% grid on
% xlim([0 t_xlim])
% 
% %% -------- KE Growth rate fitting --------
% n_KE_theory = 2 * n_viscous;
% fprintf('\n===== KE Growth Rate Fitting (linear regime only) =====\n');
% fprintf('%-25s | %-10s | %-10s | %-8s\n', 'Scheme', 'n_fit', 'n_theory', 'Err (%)');
% fprintf('------------------------------------------------------\n');
% 
% for s = 1:n_schemes
%     if isempty(csv_KE{s}); continue; end
% 
%     log_KE = log(csv_KE{s});
%     mask = (csv_time{s} >= t_fit_start) & (csv_time{s} <= t_fit_end) ...
%          & ~isnan(log_KE) & ~isinf(log_KE);
% 
%     if sum(mask) < 5; continue; end
% 
%     p = polyfit(csv_time{s}(mask), log_KE(mask), 1);
%     n_fit_KE = p(1);
%     err_pct = abs(n_fit_KE - n_KE_theory) / n_KE_theory * 100;
% 
%     fprintf('%-25s | %.4f    | %.4f    | %.2f%%\n', ...
%         scheme_labels{s}, n_fit_KE, n_KE_theory, err_pct);
% end
% fprintf('======================================================\n');
% 
% %% -------- Figure 5: KE Liquid / KE Gas Ratio --------
% figure('Color','w', 'Position', [100 100 700 500])
% hold on
% 
% for s = 1:n_schemes
%     if isempty(csv_KE_l{s}) || isempty(csv_KE_g{s}); continue; end
%     ratio = csv_KE_l{s} ./ csv_KE_g{s};
%     plot(csv_time{s}, ratio, '-', 'Color', colors(s,:), ...
%          'LineWidth', LW2, 'DisplayName', scheme_labels{s})
% end
% 
% % Theoretical ratio = rho_l / rho_g
% yline(rho_l / rho_g, 'k--', 'LineWidth', LW1, ...
%       'DisplayName', sprintf('$\\rho_l / \\rho_g = %.2f$', rho_l/rho_g))
% 
% xline(t_linear_end, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
%       'DisplayName', 'Linear regime limit')
% 
% xlabel('$t$ (s)', 'Interpreter', 'latex')
% ylabel('$\mathrm{KE}_l / \mathrm{KE}_g$', 'Interpreter', 'latex')
% title('RT Instability: KE Liquid/Gas Ratio', 'Interpreter', 'latex')
% legend('Location', 'best', 'Interpreter', 'latex')
% set(gca, 'FontSize', FS, 'FontName', 'Times New Roman', 'LineWidth', 1.7)
% grid on
% xlim([0 t_xlim])
% 
% %% -------- Figure 7: Per-phase Vertical Momentum --------
% figure('Color','w', 'Position', [100 100 700 500])
% hold on
% 
% for s = 1:n_schemes
%     if isempty(csv_rhoV_l{s}); continue; end
%     plot(csv_time{s}, csv_rhoV_l{s}, '-', 'Color', colors(s,:), ...
%          'LineWidth', LW2, 'DisplayName', [scheme_labels{s} ' (liquid)'])
%     plot(csv_time{s}, csv_rhoV_g{s}, '--', 'Color', colors(s,:), ...
%          'LineWidth', LW2, 'DisplayName', [scheme_labels{s} ' (gas)'])
% end
% 
% xline(t_linear_end, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
%       'DisplayName', 'Linear regime limit')
% 
% xlabel('$t$ (s)', 'Interpreter', 'latex')
% ylabel('$\int \rho_\phi V_\phi \, dV$', 'Interpreter', 'latex')
% title('RT Instability: Per-Phase Vertical Momentum', 'Interpreter', 'latex')
% legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 10)
% set(gca, 'FontSize', FS, 'FontName', 'Times New Roman', 'LineWidth', 1.7)
% grid on
% xlim([0 t_xlim])
% 
% %% -------- Figure 8: Vertical Momentum Ratio |rhoV_l| / |rhoV_g| --------
% figure('Color','w', 'Position', [100 100 700 500])
% hold on
% 
% for s = 1:n_schemes
%     if isempty(csv_rhoV_l{s}) || isempty(csv_rhoV_g{s}); continue; end
%     ratio = abs(csv_rhoV_l{s}) ./ abs(csv_rhoV_g{s});
%     plot(csv_time{s}, ratio, '-', 'Color', colors(s,:), ...
%          'LineWidth', LW2, 'DisplayName', scheme_labels{s})
% end
% 
% % Theoretical ratio = rho_l / rho_g
% yline(rho_l / rho_g, 'k--', 'LineWidth', LW1, ...
%       'DisplayName', sprintf('$\\rho_l / \\rho_g = %.2f$', rho_l/rho_g))
% 
% xline(t_linear_end, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
%       'DisplayName', 'Linear regime limit')
% 
% xlabel('$t$ (s)', 'Interpreter', 'latex')
% ylabel('$|\rho V|_l / |\rho V|_g$', 'Interpreter', 'latex')
% title('RT Instability: Vertical Momentum Ratio', 'Interpreter', 'latex')
% legend('Location', 'best', 'Interpreter', 'latex')
% set(gca, 'FontSize', FS, 'FontName', 'Times New Roman', 'LineWidth', 1.7)
% grid on
% xlim([0 t_xlim])
% 
% 
% %% -------- Figure 8: Per-phase Horizontal Momentum --------
% figure('Color','w', 'Position', [100 100 700 500])
% hold on
% 
% for s = 1:n_schemes
%     if isempty(csv_rhoU_l{s}); continue; end
%     plot(csv_time{s}, csv_rhoU_l{s}, '-', 'Color', colors(s,:), ...
%          'LineWidth', LW2, 'DisplayName', [scheme_labels{s} ' (liquid)'])
%     plot(csv_time{s}, csv_rhoU_g{s}, '--', 'Color', colors(s,:), ...
%          'LineWidth', LW2, 'DisplayName', [scheme_labels{s} ' (gas)'])
% end
% 
% xline(t_linear_end, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
%       'DisplayName', 'Linear regime limit')
% 
% xlabel('$t$ (s)', 'Interpreter', 'latex')
% ylabel('$\int \rho_\phi U_\phi \, dV$', 'Interpreter', 'latex')
% title('RT Instability: Per-Phase Horizontal Momentum', 'Interpreter', 'latex')
% legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 10)
% set(gca, 'FontSize', FS, 'FontName', 'Times New Roman', 'LineWidth', 1.7)
% grid on
% xlim([0 t_xlim])
% 
% %% -------- Figure 9: Horizontal Momentum Ratio |rhoU_l| / |rhoU_g| --------
% figure('Color','w', 'Position', [100 100 700 500])
% hold on
% 
% for s = 1:n_schemes
%     if isempty(csv_rhoU_l{s}) || isempty(csv_rhoU_g{s}); continue; end
%     ratio = abs(csv_rhoU_l{s}) ./ abs(csv_rhoU_g{s});
%     plot(csv_time{s}, ratio, '-', 'Color', colors(s,:), ...
%          'LineWidth', LW2, 'DisplayName', scheme_labels{s})
% end
% 
% % Theoretical ratio = rho_l / rho_g
% yline(rho_l / rho_g, 'k--', 'LineWidth', LW1, ...
%       'DisplayName', sprintf('$\\rho_l / \\rho_g = %.2f$', rho_l/rho_g))
% 
% xline(t_linear_end, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
%       'DisplayName', 'Linear regime limit')
% 
% xlabel('$t$ (s)', 'Interpreter', 'latex')
% ylabel('$|\rho U|_l / |\rho U|_g$', 'Interpreter', 'latex')
% title('RT Instability: Horizontal Momentum Ratio', 'Interpreter', 'latex')
% legend('Location', 'best', 'Interpreter', 'latex')
% set(gca, 'FontSize', FS, 'FontName', 'Times New Roman', 'LineWidth', 1.7)
% grid on
% xlim([0 t_xlim])
% 
% 
% %% -------- Figure 11: Total Potential Energy (log) --------
% figure('Color','w', 'Position', [100 100 700 500])
% hold on
% 
% for s = 1:n_schemes
%     if isempty(csv_PE_l{s}); continue; end
%     plot(csv_time{s}, log(csv_PE_l{s}(1) + csv_PE_g{s}(1)-csv_PE_l{s} -csv_PE_g{s}), '-', 'Color', colors(s,:), ...
%          'LineWidth', LW2, 'DisplayName', scheme_labels{s})
% end
% 
% % Per-scheme straight-line theory with slope 2*n_viscous anchored at t_ref=0.2
% % (PE(0)-PE(t) ~ exp(2*n*t), so log slope = 2*n)
% t_ref = 0.3;
% for s = 1:n_schemes
%     if isempty(csv_PE_l{s}); continue; end
%     [~, idx_s] = min(abs(csv_time{s} - t_ref));
%     logPE_ref = log(csv_PE_l{s}(1) + csv_PE_g{s}(1) - csv_PE_l{s}(idx_s) - csv_PE_g{s}(idx_s));
%     theory_line = logPE_ref + 2*n_viscous * (t_th - t_ref);
%     hv = 'off'; if s == 1; hv = 'on'; end
%     plot(t_th, theory_line, 'k--', 'LineWidth', LW1, ...
%          'HandleVisibility', hv, ...
%          'DisplayName', sprintf('Theory $2n=%.3f$ s$^{-1}$', 2*n_viscous))
% end
% 
% xline(t_linear_end, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
%       'DisplayName', 'Linear regime limit')
% 
% xlabel('$t$ (s)', 'Interpreter', 'latex')
% ylabel('$\ln(PE_l + PE_g)$', 'Interpreter', 'latex')
% title('RT Instability: Log Total Potential Energy', 'Interpreter', 'latex')
% legend('Location', 'best', 'Interpreter', 'latex')
% set(gca, 'FontSize', FS, 'FontName', 'Times New Roman', 'LineWidth', 1.7)
% grid on
% xlim([0 t_xlim])
% ylim([-16,-8])
% 
% %% -------- PE Growth rate fitting --------
% % PE(0)-PE(t) grows ~ exp(2*n*t) in the linear regime, so slope of
% % log(PE(0)-PE(t)) should be 2*n_viscous
% t_fit_start = 0.2;   % skip initial transient
% t_fit_end   = 0.45;
% n_PE_theory = 2 * n_viscous;
% fprintf('\n===== PE Growth Rate Fitting (linear regime only) =====\n');
% fprintf('%-25s | %-10s | %-10s | %-8s\n', 'Scheme', 'n_fit', 'n_theory', 'Err (%)');
% fprintf('------------------------------------------------------\n');
% 
% for s = 1:n_schemes
%     if isempty(csv_PE_l{s}); continue; end
% 
%     PE_drop = csv_PE_l{s}(1) + csv_PE_g{s}(1) - csv_PE_l{s} - csv_PE_g{s};
%     log_PE_drop = log(PE_drop);
%     mask = (csv_time{s} >= t_fit_start) & (csv_time{s} <= t_fit_end) ...
%          & ~isnan(log_PE_drop) & ~isinf(log_PE_drop);
% 
%     if sum(mask) < 5; continue; end
% 
%     p = polyfit(csv_time{s}(mask), log_PE_drop(mask), 1);
%     n_fit_PE = p(1);
%     err_pct = abs(n_fit_PE - n_PE_theory) / n_PE_theory * 100;
% 
%     fprintf('%-25s | %.4f    | %.4f    | %.2f%%\n', ...
%         scheme_labels{s}, n_fit_PE, n_PE_theory, err_pct);
% end
% fprintf('======================================================\n');
% 
% %% -------- Figure: KE_y / KE_x Ratio (Total) --------
% figure('Color','w', 'Position', [100 100 700 500])
% hold on
% 
% for s = 1:n_schemes
%     if isempty(csv_KE_x{s}); continue; end
%     plot(csv_time{s}, csv_KE_y{s} ./ csv_KE_x{s}, '-', 'Color', colors(s,:), ...
%          'LineWidth', LW2, 'DisplayName', scheme_labels{s})
% end
% 
% xline(t_linear_end, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
%       'DisplayName', 'Linear regime limit')
% 
% xlabel('$t$ (s)', 'Interpreter', 'latex')
% ylabel('$\mathrm{KE}_y / \mathrm{KE}_x$', 'Interpreter', 'latex')
% title('RT Instability: Total KE$_y$ / KE$_x$ Ratio', 'Interpreter', 'latex')
% legend('Location', 'best', 'Interpreter', 'latex')
% set(gca, 'FontSize', FS, 'FontName', 'Times New Roman', 'LineWidth', 1.7)
% grid on
% xlim([0 t_xlim])
% 
% %% -------- Figure: KE_y / KE_x Ratio (Liquid) --------
% figure('Color','w', 'Position', [100 100 700 500])
% hold on
% 
% for s = 1:n_schemes
%     if isempty(csv_KE_xl{s}); continue; end
%     plot(csv_time{s}, csv_KE_yl{s} ./ csv_KE_xl{s}, '-', 'Color', colors(s,:), ...
%          'LineWidth', LW2, 'DisplayName', scheme_labels{s})
% end
% 
% xline(t_linear_end, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
%       'DisplayName', 'Linear regime limit')
% 
% xlabel('$t$ (s)', 'Interpreter', 'latex')
% ylabel('$\mathrm{KE}_{y,l} / \mathrm{KE}_{x,l}$', 'Interpreter', 'latex')
% title('RT Instability: Liquid KE$_y$ / KE$_x$ Ratio', 'Interpreter', 'latex')
% legend('Location', 'best', 'Interpreter', 'latex')
% set(gca, 'FontSize', FS, 'FontName', 'Times New Roman', 'LineWidth', 1.7)
% grid on
% xlim([0 t_xlim])
% 
% %% -------- Figure: KE_y / KE_x Ratio (Gas) --------
% figure('Color','w', 'Position', [100 100 700 500])
% hold on
% 
% for s = 1:n_schemes
%     if isempty(csv_KE_xg{s}); continue; end
%     plot(csv_time{s}, csv_KE_yg{s} ./ csv_KE_xg{s}, '-', 'Color', colors(s,:), ...
%          'LineWidth', LW2, 'DisplayName', scheme_labels{s})
% end
% 
% xline(t_linear_end, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
%       'DisplayName', 'Linear regime limit')
% 
% xlabel('$t$ (s)', 'Interpreter', 'latex')
% ylabel('$\mathrm{KE}_{y,g} / \mathrm{KE}_{x,g}$', 'Interpreter', 'latex')
% title('RT Instability: Gas KE$_y$ / KE$_x$ Ratio', 'Interpreter', 'latex')
% legend('Location', 'best', 'Interpreter', 'latex')
% set(gca, 'FontSize', FS, 'FontName', 'Times New Roman', 'LineWidth', 1.7)
% grid on
% xlim([0 t_xlim])
% 
% %% -------- Figure: Liquid/Gas KE Ratio in x-direction --------
% figure('Color','w', 'Position', [100 100 700 500])
% hold on
% 
% for s = 1:n_schemes
%     if isempty(csv_KE_xl{s}); continue; end
%     plot(csv_time{s}, csv_KE_xl{s} ./ csv_KE_xg{s}, '-', 'Color', colors(s,:), ...
%          'LineWidth', LW2, 'DisplayName', scheme_labels{s})
% end
% 
% yline(rho_l / rho_g, 'k--', 'LineWidth', LW1, ...
%       'DisplayName', sprintf('$\\rho_l / \\rho_g = %.2f$', rho_l/rho_g))
% 
% xline(t_linear_end, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
%       'DisplayName', 'Linear regime limit')
% 
% xlabel('$t$ (s)', 'Interpreter', 'latex')
% ylabel('$\mathrm{KE}_{x,l} / \mathrm{KE}_{x,g}$', 'Interpreter', 'latex')
% title('RT Instability: Liquid/Gas KE Ratio ($x$-direction)', 'Interpreter', 'latex')
% legend('Location', 'best', 'Interpreter', 'latex')
% set(gca, 'FontSize', FS, 'FontName', 'Times New Roman', 'LineWidth', 1.7)
% grid on
% xlim([0 t_xlim])
% 
% %% -------- Figure: Liquid/Gas KE Ratio in y-direction --------
% figure('Color','w', 'Position', [100 100 700 500])
% hold on
% 
% for s = 1:n_schemes
%     if isempty(csv_KE_yl{s}); continue; end
%     plot(csv_time{s}, csv_KE_yl{s} ./ csv_KE_yg{s}, '-', 'Color', colors(s,:), ...
%          'LineWidth', LW2, 'DisplayName', scheme_labels{s})
% end
% 
% yline(rho_l / rho_g, 'k--', 'LineWidth', LW1, ...
%       'DisplayName', sprintf('$\\rho_l / \\rho_g = %.2f$', rho_l/rho_g))
% 
% xline(t_linear_end, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
%       'DisplayName', 'Linear regime limit')
% 
% xlabel('$t$ (s)', 'Interpreter', 'latex')
% ylabel('$\mathrm{KE}_{y,l} / \mathrm{KE}_{y,g}$', 'Interpreter', 'latex')
% title('RT Instability: Liquid/Gas KE Ratio ($y$-direction)', 'Interpreter', 'latex')
% legend('Location', 'best', 'Interpreter', 'latex')
% set(gca, 'FontSize', FS, 'FontName', 'Times New Roman', 'LineWidth', 1.7)
% grid on
% xlim([0 t_xlim])
% 
% %% -------- Figure: Liquid KE_x / KE (total) --------
% figure('Color','w', 'Position', [100 100 700 500])
% hold on
% 
% for s = 1:n_schemes
%     if isempty(csv_KE_xl{s}) || isempty(csv_KE{s}); continue; end
%     plot(csv_time{s}, csv_KE_xl{s} ./ csv_KE{s}, '-', 'Color', colors(s,:), ...
%          'LineWidth', LW2, 'DisplayName', scheme_labels{s})
% end
% 
% xline(t_linear_end, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
%       'DisplayName', 'Linear regime limit')
% 
% xlabel('$t$ (s)', 'Interpreter', 'latex')
% ylabel('$\mathrm{KE}_{x,l} / \mathrm{KE}$', 'Interpreter', 'latex')
% title('RT Instability: Liquid KE$_x$ / Total KE', 'Interpreter', 'latex')
% legend('Location', 'best', 'Interpreter', 'latex')
% set(gca, 'FontSize', FS, 'FontName', 'Times New Roman', 'LineWidth', 1.7)
% grid on
% xlim([0 t_xlim])
% 
% %% -------- Figure: Liquid KE_y / KE (total) --------
% figure('Color','w', 'Position', [100 100 700 500])
% hold on
% 
% for s = 1:n_schemes
%     if isempty(csv_KE_yl{s}) || isempty(csv_KE{s}); continue; end
%     plot(csv_time{s}, csv_KE_yl{s} ./ csv_KE{s}, '-', 'Color', colors(s,:), ...
%          'LineWidth', LW2, 'DisplayName', scheme_labels{s})
% end
% 
% xline(t_linear_end, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
%       'DisplayName', 'Linear regime limit')
% 
% xlabel('$t$ (s)', 'Interpreter', 'latex')
% ylabel('$\mathrm{KE}_{y,l} / \mathrm{KE}$', 'Interpreter', 'latex')
% title('RT Instability: Liquid KE$_y$ / Total KE', 'Interpreter', 'latex')
% legend('Location', 'best', 'Interpreter', 'latex')
% set(gca, 'FontSize', FS, 'FontName', 'Times New Roman', 'LineWidth', 1.7)
% grid on
% xlim([0 t_xlim])
% 
% %% -------- Figure: Gas KE_x / KE (total) --------
% figure('Color','w', 'Position', [100 100 700 500])
% hold on
% 
% for s = 1:n_schemes
%     if isempty(csv_KE_xg{s}) || isempty(csv_KE{s}); continue; end
%     plot(csv_time{s}, csv_KE_xg{s} ./ csv_KE{s}, '-', 'Color', colors(s,:), ...
%          'LineWidth', LW2, 'DisplayName', scheme_labels{s})
% end
% 
% xline(t_linear_end, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
%       'DisplayName', 'Linear regime limit')
% 
% xlabel('$t$ (s)', 'Interpreter', 'latex')
% ylabel('$\mathrm{KE}_{x,g} / \mathrm{KE}$', 'Interpreter', 'latex')
% title('RT Instability: Gas KE$_x$ / Total KE', 'Interpreter', 'latex')
% legend('Location', 'best', 'Interpreter', 'latex')
% set(gca, 'FontSize', FS, 'FontName', 'Times New Roman', 'LineWidth', 1.7)
% grid on
% xlim([0 t_xlim])
% 
% %% -------- Figure: Gas KE_y / KE (total) --------
% figure('Color','w', 'Position', [100 100 700 500])
% hold on
% 
% for s = 1:n_schemes
%     if isempty(csv_KE_yg{s}) || isempty(csv_KE{s}); continue; end
%     plot(csv_time{s}, csv_KE_yg{s} ./ csv_KE{s}, '-', 'Color', colors(s,:), ...
%          'LineWidth', LW2, 'DisplayName', scheme_labels{s})
% end
% 
% xline(t_linear_end, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
%       'DisplayName', 'Linear regime limit')
% 
% xlabel('$t$ (s)', 'Interpreter', 'latex')
% ylabel('$\mathrm{KE}_{y,g} / \mathrm{KE}$', 'Interpreter', 'latex')
% title('RT Instability: Gas KE$_y$ / Total KE', 'Interpreter', 'latex')
% legend('Location', 'best', 'Interpreter', 'latex')
% set(gca, 'FontSize', FS, 'FontName', 'Times New Roman', 'LineWidth', 1.7)
% grid on
% xlim([0 t_xlim])
% 
% % figure('Color','w', 'Position', [100 100 700 500])
% % hold on
% % 
% % for s = 1:n_schemes
% %     if isempty(csv_EN_l{s}); continue; end
% %     plot(csv_time{s}, log(csv_EN_l{s} + csv_EN_g{s}), '-', 'Color', colors(s,:), ...
% %          'LineWidth', LW2, 'DisplayName', scheme_labels{s})
% % end
% % 
% % % Per-scheme straight-line theory with slope 2*n_viscous anchored at t_ref
% % % (Enstrophy ~ |grad u|^2 ~ exp(2*n*t) in the linear regime)
% % t_ref = 0.2;
% % for s = 1:n_schemes
% %     if isempty(csv_EN_l{s}); continue; end
% %     [~, idx_s] = min(abs(csv_time{s} - t_ref));
% %     logEN_ref = log(csv_EN_l{s}(idx_s) + csv_EN_g{s}(idx_s));
% %     theory_line = logEN_ref + 2*n_viscous * (t_th - t_ref);
% %     hv = 'off'; if s == 1; hv = 'on'; end
% %     plot(t_th, theory_line, 'k--', 'LineWidth', LW1, ...
% %          'HandleVisibility', hv, ...
% %          'DisplayName', sprintf('Theory $2n=%.3f$ s$^{-1}$', 2*n_viscous))
% % end
% % 
% % xline(t_linear_end, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
% %       'DisplayName', 'Linear regime limit')
% % 
% % xlabel('$t$ (s)', 'Interpreter', 'latex')
% % ylabel('$\ln(\Omega_l + \Omega_g)$', 'Interpreter', 'latex')
% % title('RT Instability: Log Total Enstrophy', 'Interpreter', 'latex')
% % legend('Location', 'northwest', 'Interpreter', 'latex')
% % set(gca, 'FontSize', FS, 'FontName', 'Times New Roman', 'LineWidth', 1.7)
% % grid on
% % xlim([0 t_xlim])
% 
% % %% -------- Enstrophy Growth rate fitting --------
% % n_EN_theory = 2 * n_viscous;
% % fprintf('\n===== Enstrophy Growth Rate Fitting (linear regime only) =====\n');
% % fprintf('%-25s | %-10s | %-10s | %-8s\n', 'Scheme', 'n_fit', 'n_theory', 'Err (%)');
% % fprintf('--------------------------------------------------------------\n');
% % 
% % for s = 1:n_schemes
% %     if isempty(csv_EN_l{s}); continue; end
% % 
% %     log_EN = log(csv_EN_l{s} + csv_EN_g{s});
% %     mask = (csv_time{s} >= t_fit_start) & (csv_time{s} <= t_fit_end) ...
% %          & ~isnan(log_EN) & ~isinf(log_EN);
% % 
% %     if sum(mask) < 5; continue; end
% % 
% %     p = polyfit(csv_time{s}(mask), log_EN(mask), 1);
% %     n_fit_EN = p(1);
% %     err_pct = abs(n_fit_EN - n_EN_theory) / n_EN_theory * 100;
% % 
% %     fprintf('%-25s | %.4f    | %.4f    | %.2f%%\n', ...
% %         scheme_labels{s}, n_fit_EN, n_EN_theory, err_pct);
% % end
% % fprintf('==============================================================\n');
% % % 
% % %% -------- Figure 14: Total Mechanical Energy (KE + PE) --------
% % figure('Color','w', 'Position', [100 100 700 500])
% % hold on
% % 
% % for s = 1:n_schemes
% %     if isempty(csv_KE{s}) || isempty(csv_PE_l{s}); continue; end
% % 
% %     % Total mechanical energy = KE + PE_l + PE_g
% %     tot_E = csv_KE{s} + csv_PE_l{s} + csv_PE_g{s};
% % 
% %     % Normalize as drift from initial value
% %     tot_E0 = tot_E(1);
% %     if abs(tot_E0) > eps
% %         tot_E_norm = (tot_E - tot_E0) / abs(tot_E0);
% %     else
% %         tot_E_norm = tot_E - tot_E0;
% %     end
% % 
% %     plot(csv_time{s}, tot_E_norm, '-', 'Color', colors(s,:), ...
% %          'LineWidth', LW2, 'DisplayName', scheme_labels{s})
% % end
% % 
% % xline(t_linear_end, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
% %       'DisplayName', 'Linear regime limit')
% % 
% % xlabel('$t$ (s)', 'Interpreter', 'latex')
% % ylabel('$(E_{tot} - E_{tot,0}) / |E_{tot,0}|$', 'Interpreter', 'latex')
% % title('RT Instability: Total Mechanical Energy Drift (KE + PE)', 'Interpreter', 'latex')
% % legend('Location', 'best', 'Interpreter', 'latex')
% % set(gca, 'FontSize', FS, 'FontName', 'Times New Roman', 'LineWidth', 1.7)
% % grid on
% % xlim([0 t_xlim])
% % 
