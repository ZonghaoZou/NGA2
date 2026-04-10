%% Rayleigh–Taylor Instability: Nonlinear Regime Analysis
%  Compares bubble/spike terminal velocities against Goncharov (2002)
%  "Analytical Model of Nonlinear, Single-Mode, Classical Rayleigh–Taylor
%   Instability at Arbitrary Atwood Numbers"
%
%  Goncharov terminal bubble velocity (2D):
%    V_b = sqrt( 2 * At * g / ((1 + At) * C_g * k) )
%    where C_g = 3 for 2D, C_g = 1 for 3D
%
%  Usage: Run from the RayleighTaylor/ directory.
clear; clc; close all;

%% -------- Physical parameters --------
rho_l = 1.225;       % heavy (liquid, VF=1)
rho_g = 0.1694;      % light (gas,    VF=0)
mu_l    = 0.00313;     % dynamic viscosity (both phases)
mu_g    = 0.000432834;     % dynamic viscosity (both phases)
g_val = 9.81;        % gravity magnitude
A0    = 0.012;       % perturbation amplitude
Lx    = 1.0;         % domain width
sigma = 0.0;         % surface tension
viscoutype = "arith_stats_256_sharp_same";

At = (rho_l - rho_g) / (rho_l + rho_g);   % Atwood number
k  = 2*pi / Lx;                            % wavenumber

fprintf('Atwood number At = %.4f\n', At);
fprintf('Wavenumber k     = %.4f (2*pi)\n', k);
blue=[71,135,224]/255;
red =[218,62,32]/255;
green= [82,147,47]/255;

%% -------- Goncharov (2002) terminal velocities --------
% 2D geometry: C_g = 3
C_g = 3;
nuh=mu_l/rho_l;
% V_b_Goncharov = sqrt(2 * At * g_val / ((1 + At) * C_g * k));
V_b_Sohn = -2*k*nuh/3 + sqrt(2 * At * g_val / ((1 + At) * C_g * k)  +(4*k^2*nuh^2)/9 -k*sigma/(9*rho_l));

% Spike terminal velocity (Goncharov Eq. for spike):
%   V_s = sqrt( 2 * At * g / ((1 - At) * C_g * k) )
V_s_Goncharov = sqrt(2 * At * g_val / ((1 - At) * C_g * k));

fprintf('Goncharov bubble velocity V_b = %.4f m/s\n', V_b_Goncharov);
fprintf('Goncharov spike  velocity V_s = %.4f m/s\n', V_s_Goncharov);

%% -------- Formatting --------
LW1 = 1.5;   % theory line width
LW2 = 1.7;   % simulation line width
FS  = 15;    % font size

%% -------- Scheme definitions --------
scheme_dirs   = {strcat('result_',viscoutype,'/1'), strcat('result_',viscoutype,'/2'), strcat('result_',viscoutype,'/3')};
scheme_labels = {'Default NGA2 (SG)', 'KE Conservative (SG)', 'SL Momentum (CG)'};
n_schemes     = numel(scheme_dirs);
colors = vertcat(vertcat(blue,red),green);

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

%% -------- Compute velocities via finite differencing --------
bub_vel = cell(n_schemes, 1);
spike_vel = cell(n_schemes, 1);
time_vel = cell(n_schemes, 1);

for s = 1:n_schemes
    if isempty(bub_all{s}); continue; end
    
    t = time_all{s};
    dt = diff(t);
    
    % Central differencing (interior), forward/backward at ends
    db = diff(bub_all{s});
    ds = diff(spike_all{s});
    
    bub_vel{s}   = db ./ dt;         % bubble velocity (positive = rising)
    spike_vel{s} = -ds ./ dt;        % spike velocity (positive = falling, negate)
    time_vel{s}  = (t(1:end-1) + t(2:end)) / 2;  % midpoint times
    
    % Smooth with moving average for cleaner signal
    N_smooth = 20;
    if length(bub_vel{s}) > 2*N_smooth
        bub_vel{s}   = movmean(bub_vel{s}, N_smooth);
        spike_vel{s} = movmean(spike_vel{s}, N_smooth);
    end
end

%% -------- Figure 1: Bubble and Spike positions --------
figure('Color','w', 'Position', [100 100 700 500])
hold on

for s = 1:n_schemes
    if isempty(bub_all{s}); continue; end
    plot(time_all{s}, bub_all{s}, '-', 'Color', colors(s,:), ...
         'LineWidth', LW2, 'DisplayName', [scheme_labels{s} ' (bubble)'])
    plot(time_all{s}, spike_all{s}, '--', 'Color', colors(s,:), ...
         'LineWidth', LW2, 'DisplayName', [scheme_labels{s} ' (spike)'])
end

xlabel('$t$ (s)', 'Interpreter', 'latex')
ylabel('$y$ position', 'Interpreter', 'latex')
title('RT Instability: Bubble and Spike Positions', 'Interpreter', 'latex')
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 10)
set(gca, 'FontSize', FS, 'FontName', 'Times New Roman', 'LineWidth', 1.7)
grid on

%% -------- Figure 2: Bubble velocity vs Goncharov --------
figure('Color','w', 'Position', [100 100 700 500])
hold on

for s = 1:n_schemes
    if isempty(bub_vel{s}); continue; end
    plot(time_vel{s}, bub_vel{s}, '-', 'Color', colors(s,:), ...
         'LineWidth', LW2, 'DisplayName', scheme_labels{s})
end

% Goncharov terminal velocity
yline(V_b_Goncharov, 'k--', 'LineWidth', LW1, ...
      'DisplayName', sprintf('Goncharov $V_b = %.4f$ m/s', V_b_Goncharov))

xlabel('$t$ (s)', 'Interpreter', 'latex')
ylabel('Bubble velocity $\dot{y}_b$ (m/s)', 'Interpreter', 'latex')
ylim([0,0.8])
title('RT Instability: Bubble Velocity vs Goncharov (2002)', 'Interpreter', 'latex')
legend('Location', 'best', 'Interpreter', 'latex')
set(gca, 'FontSize', FS, 'FontName', 'Times New Roman', 'LineWidth', 1.7)
grid on

%% -------- Figure 3: Spike velocity vs Goncharov --------
figure('Color','w', 'Position', [100 100 700 500])
hold on

for s = 1:n_schemes
    if isempty(spike_vel{s}); continue; end
    plot(time_vel{s}, spike_vel{s}, '-', 'Color', colors(s,:), ...
         'LineWidth', LW2, 'DisplayName', scheme_labels{s})
end

xlabel('$t$ (s)', 'Interpreter', 'latex')
ylabel('Spike velocity $|\dot{y}_s|$ (m/s)', 'Interpreter', 'latex')
title('RT Instability: Spike Velocity', 'Interpreter', 'latex')
legend('Location', 'best', 'Interpreter', 'latex')
set(gca, 'FontSize', FS, 'FontName', 'Times New Roman', 'LineWidth', 1.7)
grid on

%% -------- Terminal velocity summary --------
fprintf('\n===== Bubble Terminal Velocity Comparison =====\n');
fprintf('%-25s | %-12s | %-12s | %-8s\n', 'Scheme', 'V_b (sim)', 'V_b (Gonch)', 'Err (%)');
fprintf('------------------------------------------------------\n');

% Estimate terminal velocity as average over late-time plateau
for s = 1:n_schemes
    if isempty(bub_vel{s}); continue; end
    
    % Use last 20% of data for terminal velocity estimate
    N = length(bub_vel{s});
    idx_late = round(0.8*N):N;
    
    Vb_sim = mean(bub_vel{s}(idx_late));
    err_pct = abs(Vb_sim - V_b_Goncharov) / V_b_Goncharov * 100;
    
    fprintf('%-25s | %.4f      | %.4f      | %.2f%%\n', ...
        scheme_labels{s}, Vb_sim, V_b_Goncharov, err_pct);
end
fprintf('======================================================\n');
