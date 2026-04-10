%% Oscillating Droplet: Scheme Comparison vs. Lamb Theory (3D)
%  Reads monitor/simulation from each scheme folder and compares 
%  oscillation of semi-axis against the analytical result.
%
%  Model: XMAX(t) = R + A0*exp(-gamma*t)*cos(omega*t)
%         where omega is the 3D Lamb frequency (mode l=2)
clear
close all

%% ===== Hardcoded parameters (from input file) =====
R       = 2.0;        % Unperturbed radius (D/2 = 4/2)
lambda  = 0.01;       % Perturbation amplitude
sigma   = 1.0;        % Surface tension
rho_l   = 1.0;        % Liquid density
rho_g   = 0.01;       % Gas density
mu_l    = 1.0e-2;     % Liquid dynamic viscosity
mu_g    = 1.0e-4;     % Gas dynamic viscosity
l_mode  = 2;          % Mode number

%% ===== Theoretical values =====
gamma=2*R*((l_mode+1)*rho_l + l_mode*rho_g)*(sqrt(mu_l*rho_l)+sqrt(mu_g*rho_g));
% gamma_th = (2*l_mode+1)*(l_mode-1) * mu_l / (rho_l * R^2);
% Lamb angular frequency (3D, mode l):
omega = sqrt( l_mode*(l_mode-1)*(l_mode+1)*(l_mode+2) * sigma ...
               / ( (l_mode+1)*rho_l + l_mode*rho_g ) / R^3 );
omega_th = omega- ((2*l_mode+1)^2)*(sqrt(omega*rho_l*rho_g*mu_g*mu_l)/(sqrt(2)*gamma));
% Viscous damping rate (Lamb approximation):
gamma_th = (2*l_mode+1)^2*sqrt(omega*rho_l*rho_g*mu_g*mu_l)/(sqrt(2)*gamma) - (2*l_mode+1)^4*rho_l*rho_g*mu_g*mu_l/gamma^2 + (2*l_mode+1)*(2*(l_mode-1)*(l_mode+1)*mu_l^2*rho_l+2*l_mode*(l_mode+2)*mu_g^2*rho_g+mu_g*mu_l*((l_mode+2)*rho_l-(l_mode-1)*rho_g))/(R*gamma*(sqrt(mu_l*rho_l)+sqrt(mu_g*rho_g)));
% Initial amplitude
A0 = lambda * R;

fprintf('===== Lamb Theory (l=%d, 3D) =====\n', l_mode);
fprintf('  omega = %.6f rad/s\n', omega_th);
fprintf('  gamma = %.6f 1/s\n', gamma_th);
fprintf('  T     = %.4f\n', 2*pi/omega_th);
fprintf('  A0    = %.4f\n', A0);
fprintf('====================================\n\n');

%% ===== Theoretical curve =====
t_th = 0:0.1:20;
Aex_x = R + A0*exp(-gamma_th*t_th).*cos(omega_th*t_th);
Aex_y = R - (A0/2)*exp(-gamma_th*t_th).*cos(omega_th*t_th);

%% ===== Read simulation data =====
% Define your base schemes and resolutions
base_schemes = { 1, 'Default NGA2 (SG)'; ...
                 2, 'KE Conservative (SG)'; ...
                 3, 'SL Momentum (CG)' };

resolutions  = [32, 64, 128, 256]; % Add 256 here when ready

% We will dynamically build the list of paths to load
scheme_dirs   = {};
scheme_labels = {};

for i = 1:size(base_schemes, 1)
    id   = base_schemes{i, 1};
    name = base_schemes{i, 2};
    
    for res = resolutions
        % Construct path: result/1/simulation_32/monitor/simulation
        path_check = fullfile('result', num2str(id), sprintf('simulation_%d', res));
        
        if isfile(path_check)
            scheme_dirs{end+1}   = path_check; %#ok<SAGROW>
            scheme_labels{end+1} = sprintf('%s N=%d', name, res); %#ok<SAGROW>
        end
    end
end

n_schemes = numel(scheme_dirs);

% Storage
time_all  = cell(n_schemes, 1);
xmax_all  = cell(n_schemes, 1);
ymax_all  = cell(n_schemes, 1);
tke_all   = cell(n_schemes, 1);

for s = 1:n_schemes
    fpath = scheme_dirs{s};
    
    % Read with textscan (matching reference style)
    fileID = fopen(fpath, 'r');
    raw = readmatrix(fpath, 'FileType', 'text', 'NumHeaderLines', 2);
    fclose(fileID);
    
    start_index = 1;
    time_all{s}  = raw(start_index:end-1, 2);
    tke_all{s}   = raw(start_index+1:end, 15);
    xmax_all{s}  = abs(raw(start_index+1:end, 17));
    ymax_all{s}  = abs(raw(start_index+1:end, 18));
    
    fprintf('Loaded %d timesteps for: %s\n', numel(time_all{s}), scheme_labels{s});
end

%% ===== Figure 1: X semi-axis oscillation (All combined) =====
LW1 = 1.0;
LW2 = 2.0;
% figure('Color','w')
% hold on
% for s = 1:n_schemes
%     if ~isempty(xmax_all{s})
%         plot(time_all{s}, xmax_all{s}, 'LineWidth', LW2, 'LineStyle', '-')
%     end
% end
% plot(t_th, Aex_x, 'LineWidth', LW1, 'Marker', '*', 'Color', 'k')
% xlim([0, 20])
% ylim([R - 1.5*A0, R + 1.5*A0])
% xlabel('$t$', 'Interpreter', 'latex')
% ylabel('$x_{\max}$', 'Interpreter', 'latex')
% legend([scheme_labels, {'Lamb Theory'}], 'Location', 'southeast', 'Interpreter', 'latex');
% set(gca, 'Fontsize', 15, 'fontname', 'Times New Roman', 'LineWidth', 1.7)
% title('X Semi-Axis Oscillation (All Data)', 'Interpreter', 'latex')

%% ===== Figure 2: Y semi-axis oscillation (All combined) =====
% figure('Color','w')
% hold on
% for s = 1:n_schemes
%     if ~isempty(ymax_all{s})
%         plot(time_all{s}, ymax_all{s}, 'LineWidth', LW2, 'LineStyle', '-')
%     end
% end
% plot(t_th, Aex_y, 'LineWidth', LW1, 'Marker', '*', 'Color', 'k')
% xlim([0, 20])
% ylim([R - 1.0*A0, R + 1.0*A0])
% xlabel('$t$', 'Interpreter', 'latex')
% ylabel('$y_{\max}$', 'Interpreter', 'latex')
% legend([scheme_labels, {'Lamb Theory'}], 'Location', 'southeast', 'Interpreter', 'latex');
% set(gca, 'Fontsize', 15, 'fontname', 'Times New Roman', 'LineWidth', 1.7)
% title('Y Semi-Axis Oscillation', 'Interpreter', 'latex')

%% ===== Figure 3: TKE =====
% figure('Color','w')
% hold on
% for s = 1:n_schemes
%     if ~isempty(tke_all{s})
%         plot(time_all{s}, tke_all{s}, 'LineWidth', LW2, 'LineStyle', '-')
%     end
% end
% xlabel('$t$', 'Interpreter', 'latex')
% ylabel('Total Kinetic Energy', 'Interpreter', 'latex')
% legend(scheme_labels, 'Location', 'best', 'Interpreter', 'latex');
% set(gca, 'Fontsize', 15, 'fontname', 'Times New Roman', 'LineWidth', 1.7)
% title('Total Kinetic Energy', 'Interpreter', 'latex')

%% ===== Curve fitting & Table Generation =====
fprintf('\n===== Curve Fitting (lsqcurvefit) =====\n');
model_fun = @(params, x) params(1) + params(2) * exp(-params(3) * x) .* cos(params(4)*x);
initial_guess = [R, A0, gamma_th, omega_th];
lb = [-Inf, -Inf, 0, 0];
ub = [Inf, Inf, Inf, Inf];
opts = optimset('Display', 'off');

% Initialize column arrays
tbl_names = {};
tbl_omega_fit = []; tbl_omega_th = []; tbl_omega_err = [];
tbl_gamma_fit = []; tbl_gamma_th = []; tbl_gamma_err = [];

for s = 1:n_schemes
    if isempty(xmax_all{s}); continue; end
    fitted_params = lsqcurvefit(model_fun, initial_guess, time_all{s}, xmax_all{s}, lb, ub, opts);
    
    omega_fit = fitted_params(4);
    gamma_fit = fitted_params(3);
    R_fit     = fitted_params(1);
    A_fit     = fitted_params(2);
    
    % Store for table
    tbl_names{end+1,1} = scheme_labels{s}; %#ok<*SAGROW>
    
    % Omega (Frequency) Data
    tbl_omega_th(end+1,1)  = omega_th;       % Theoretical
    tbl_omega_fit(end+1,1) = omega_fit;      % Fitted
    tbl_omega_err(end+1,1) = abs(omega_fit/omega_th - 1) * 100; % % Error
    
    % Gamma (Damping) Data
    tbl_gamma_th(end+1,1)  = gamma_th;       % Theoretical
    tbl_gamma_fit(end+1,1) = gamma_fit;      % Fitted
    tbl_gamma_err(end+1,1) = abs(gamma_fit/gamma_th - 1) * 100; % % Error
end

% Create Table with Theoretical Columns included
VarNames = {'Scheme', ...
            'Omega_Theory', 'Omega_Fit', 'Omega_Err_Pct', ...
            'Gamma_Theory', 'Gamma_Fit', 'Gamma_Err_Pct'};
        
ResultsTable = table(tbl_names, ...
                     tbl_omega_th, tbl_omega_fit, tbl_omega_err, ...
                     tbl_gamma_th, tbl_gamma_fit, tbl_gamma_err, ...
                     'VariableNames', VarNames);
                 
fprintf('\n===== GRID CONVERGENCE ERROR SUMMARY =====\n');
disp(ResultsTable)
fprintf('==========================================\n');
%% ===== Figure 4 (Split): Shifted Simulation vs Lamb Theory =====
% We generate 3 separate figures (one per scheme) to show grid convergence cleanly.

% Define consistent colors for resolutions so N=32 is always the same color
% 1=Blue (32), 2=Red (64), 3=Yellow/Orange (128), 4=Purple (256)
res_colors = lines(length(resolutions)); 

for i = 1:size(base_schemes, 1)
    current_scheme_id   = base_schemes{i, 1};
    current_scheme_name = base_schemes{i, 2};
    
    figure('Color','w', 'Name', ['Shifted: ' current_scheme_name])
    hold on
    
    % 1. Plot Theory First
    plot(t_th, Aex_x, 'k--', 'LineWidth', LW1, 'DisplayName', 'Lamb Theory')
    
    % 2. Find all loaded data that belongs to this scheme
    for s = 1:n_schemes
        if contains(scheme_labels{s}, current_scheme_name)
            
            % Identify which resolution this is to assign the correct color
            color_idx = 1; % Default
            for r = 1:length(resolutions)
                if contains(scheme_labels{s}, ['N=' num2str(resolutions(r))])
                    color_idx = r;
                    break;
                end
            end
            
            % Re-calculate shift parameters
            if isempty(xmax_all{s}); continue; end
            fitted_params = lsqcurvefit(model_fun, initial_guess, time_all{s}, xmax_all{s}, lb, ub, opts);
            R_fit = fitted_params(1);
            A_fit = fitted_params(2);
            
            % Shift
            xmax_shifted = R + (xmax_all{s} - R_fit) * (A0 / A_fit);
            
            % Plot
            plot(time_all{s}, xmax_shifted, '-', 'Color', res_colors(color_idx,:), ...
                 'LineWidth', LW2, 'DisplayName', scheme_labels{s})
        end
    end
    
    % Formatting
    xlim([0, 20])
    ylim([R - 1.5*A0, R + 1.5*A0])
    xlabel('$t$', 'Interpreter', 'latex')
    ylabel('$x_{\max}$ (shifted)', 'Interpreter', 'latex')
    legend('Location', 'southeast', 'Interpreter', 'latex');
    set(gca, 'Fontsize', 15, 'fontname', 'Times New Roman', 'LineWidth', 1.7)
    title(['Shifted: ' current_scheme_name], 'Interpreter', 'latex')
    hold off
end

%% ===== Figure 5: Fitted vs Data (Verification - All Combined) =====
% figure('Color','w', 'Name', 'Curve Fitting Verification')
% hold on
% colors = lines(n_schemes);
% for s = 1:n_schemes
%     if isempty(xmax_all{s}); continue; end
% 
%     % Get fit again
%     fitted_params = lsqcurvefit(model_fun, initial_guess, time_all{s}, xmax_all{s}, lb, ub, opts);
% 
%     plot(time_all{s}, xmax_all{s}, 'o', 'Color', colors(s,:), ...
%          'MarkerSize', 3, 'DisplayName', [scheme_labels{s}, ' (data)'])
%     plot(time_all{s}, model_fun(fitted_params, time_all{s}), '-', ...
%          'Color', colors(s,:), 'LineWidth', LW2, ...
%          'DisplayName', [scheme_labels{s}, ' (fit)'])
% end
% plot(t_th, Aex_x, 'k--', 'LineWidth', LW1, 'DisplayName', 'Lamb Theory')
% xlim([0, 20])
% ylim([R - 1.5*A0, R + 1.5*A0])
% xlabel('$t$', 'Interpreter', 'latex')
% ylabel('$x_{\max}$', 'Interpreter', 'latex')
% set(gca, 'Fontsize', 15, 'fontname', 'Times New Roman', 'LineWidth', 1.7)
% title('Curve Fitting Verification', 'Interpreter', 'latex')