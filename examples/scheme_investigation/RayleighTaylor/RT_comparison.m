%% Rayleigh–Taylor Instability: Nonlinear Regime Analysis
%  Compares bubble terminal velocities against:
%   1) Goncharov (2002) - Inviscid
%   2) Sohn (2004) - Viscous (Based on heavy fluid viscosity)
%  Generates separate plots for each gas viscosity case.
%
%  Usage: Run from the RayleighTaylor/ directory.
clear; clc; close all;

%% -------- Physical parameters --------
rho_l = 1.225;       % heavy (liquid, VF=1)
rho_g = 0.1694;      % light (gas,    VF=0)
mu_l  = 0.00313;     % dynamic viscosity (heavy fluid)
g_val = 9.81;        % gravity magnitude
Lx    = 1.0;         % domain width
sigma = 0.0;         % surface tension

At = (rho_l - rho_g) / (rho_l + rho_g);    % Atwood number
k  = 2*pi / Lx;                            % wavenumber
C_g = 3;                                   % 2D geometry coefficient

fprintf('Atwood number At = %.4f\n', At);
fprintf('Wavenumber k     = %.4f (2*pi)\n\n', k);

%% -------- Case Definitions --------
% Define the two simulations (differing only by gas viscosity)
cases(1).dir   = 'arith_stats_256_same';
cases(1).mu_g  = 0.00313;

cases(2).dir   = 'arith_stats_256';
cases(2).mu_g  = 0.000432834;

%% -------- Theoretical Terminal Velocities --------
% 1) Goncharov terminal bubble velocity (Inviscid)
V_b_Goncharov = sqrt(2 * At * g_val / ((1 + At) * C_g * k));

% 2) Sohn terminal bubble velocity (Viscous - depends ONLY on mu_l)
nuh = mu_l / rho_l; % Kinematic viscosity of the heavy fluid
V_b_Sohn = -2*k*nuh/3 + sqrt(2 * At * g_val / ((1 + At) * C_g * k) + (4*k^2*nuh^2)/9 - k*sigma/(9*rho_l));

%% -------- Formatting & Scheme Definitions --------
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

%% -------- Read files, Compute, and Plot --------
% Array to store terminal velocities for the summary table
Vb_sim = zeros(length(cases), n_schemes);

for c = 1:length(cases)
    
    % Create a separate figure for each case
    figure('Color','w', 'Position', [100+c*50, 100+c*50, 800, 600])
    hold on
    
    title_suffix = sprintf('($\\mu_g=%g$)', cases(c).mu_g);
    
    for s = 1:n_schemes
        % Construct directory path
        mfile = fullfile(sprintf('result_%s', cases(c).dir), num2str(s), 'simulation');
        
        if ~isfile(mfile)
            warning('Monitor file not found: %s', mfile);
            continue
        end
        
        raw = readmatrix(mfile, 'FileType', 'text', 'NumHeaderLines', 2);
        
        t   = raw(2:end, 2);    % Time
        bub = raw(2:end, 16);   % Bubble Y position
        
        if isempty(bub); continue; end
        
        % Compute velocities via finite differencing
        dt = diff(t);
        db = diff(bub);
        
        bub_vel  = db ./ dt;         
        time_vel = (t(1:end-1) + t(2:end)) / 2;  
        
        % Smooth with moving average for cleaner signal
        N_smooth = 200;
        if length(bub_vel) > 2*N_smooth
            bub_vel = movmean(bub_vel, N_smooth);
        end
        
        % Plot simulation curve (using solid lines for all schemes)
        disp_name = sprintf('%s ($\\mu_g=%g$)', scheme_labels{s}, cases(c).mu_g);
        plot(time_vel, bub_vel, '-', 'Color', colors(s,:), ...
             'LineWidth', LW2, 'DisplayName', disp_name)
         
        % Estimate terminal velocity (average over last 20%)
        N_pts = length(bub_vel);
        idx_late = round(0.8 * N_pts):N_pts;
        Vb_sim(c, s) = mean(bub_vel(idx_late));
    end
    
    %% -------- Add Theoretical Lines to Current Figure --------
    % Goncharov 
    yline(V_b_Goncharov, 'k--', 'LineWidth', 2.0, ...
          'DisplayName', sprintf('Goncharov (Inviscid) $V_b = %.4f$', V_b_Goncharov));
          
    % Sohn 
    yline(V_b_Sohn, 'k-.', 'LineWidth', LW1, ...
          'DisplayName', sprintf('Sohn (Viscous, $\\mu_l$) $V_b = %.4f$', V_b_Sohn));
          
    % Finalize Figure
    xlabel('$t$ (s)', 'Interpreter', 'latex')
    ylabel('Bubble velocity $\dot{y}_b$ (m/s)', 'Interpreter', 'latex')
    title(sprintf('RT Instability: Bubble Velocity %s', title_suffix), 'Interpreter', 'latex')
    legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 10, 'NumColumns', 2)
    set(gca, 'FontSize', FS, 'FontName', 'Times New Roman', 'LineWidth', 1.7)
    grid on
    ylim([0.6, 0.7]) 
    xlim([0.6, 2.1])
end

%% -------- Terminal velocity summary --------
fprintf('===== Bubble Terminal Velocity Comparison =====\n');
fprintf('%-38s | %-10s | %-12s | %-12s | %-8s\n', 'Scheme & Case', 'V_b (sim)', 'V_b (Gonch)', 'V_b (Sohn)', 'Err (Sohn)');
fprintf('------------------------------------------------------------------------------------------\n');
for c = 1:length(cases)
    for s = 1:n_schemes
        if Vb_sim(c,s) == 0; continue; end
        
        err_pct = abs(Vb_sim(c,s) - V_b_Sohn) / V_b_Sohn * 100;
        label_str = sprintf('%s (mu_g=%.4g)', scheme_labels{s}, cases(c).mu_g);
        
        fprintf('%-38s | %.4f     | %.4f      | %.4f      | %.2f%%\n', ...
            label_str, Vb_sim(c,s), V_b_Goncharov, V_b_Sohn, err_pct);
    end
end
fprintf('==========================================================================================\n');