clear
format longG

%% -------- Configuration --------
% Scheme directories (matching run_all.sh output structure)
foldername1 = "1DefaultNGA2_SG";
foldername2 = "2SpatialTemporalKEcons_SG";
foldername3 = "3SLmomentum_CG";

% Mesh resolutions
nx_list = [16 32 64 128 256];
n_res = numel(nx_list);

% Physical parameters
D = 1.0e-4;
L = 2.5e-4;
dx_D = sqrt(2) * (L./nx_list) / D;   % CΔ = sqrt(2)*Δ/D  (same x-axis as spurious current)

%% -------- Colors & Style (matching spurious current) --------
blue   = [71,135,224]/255;
red    = [218,62,32]/255;
green  = [82,147,47]/255;
yellow = [241,178,56]/255;
LW1 = 1.0;
LW2 = 2.5;
MS  = 10;

%% -------- Read CSV data for all three schemes --------
schemes = {foldername1, foldername2, foldername3};
scheme_labels = {"Default NGA2", "KE Conservative", "SL Momentum"};
scheme_colors = {blue, red, green};
n_schemes = numel(schemes);

ca_rms = nan(n_schemes, n_res);
ca_max = nan(n_schemes, n_res);
ca_mean = nan(n_schemes, n_res);

for s = 1:n_schemes
    for r = 1:n_res
        nx = nx_list(r);
        csvfile = fullfile('result', schemes{s}, sprintf('%d.csv', nx));
        if isfile(csvfile)
            raw = readmatrix(csvfile);
            % Take max across runs (last row if multiple appended)
            ca_rms(s, r)  = max(raw(:, 1));
            ca_max(s, r)  = max(raw(:, 2));
            if size(raw, 2) >= 3
                ca_mean(s, r) = max(raw(:, 3));
            end
            fprintf('Loaded: %s, nx=%d -> Ca_rms=%.4e, Ca_max=%.4e, Ca_mean=%.4e\n', ...
                scheme_labels{s}, nx, ca_rms(s,r), ca_max(s,r), ca_mean(s,r));
        else
            fprintf('Missing: %s\n', csvfile);
        end
    end
end

%% -------- Reference slopes --------
cdx_ref = dx_D;  % Use same x-values as data
% Adjust scaling factors after seeing data; these are starting points
sy1 = 3e-3;
sy2 = 2e-5;
firstorder  = sy1 * (dx_D / dx_D(1));
secondorder = sy2 * (dx_D / dx_D(1)).^2;

%% -------- Figure 1: All schemes combined (Ca_rms & Ca_max vs CΔ) --------
figure('Color','w')
hold on

for s = 1:n_schemes
    valid = ~isnan(ca_rms(s,:));
    if any(valid)
        loglog(dx_D(valid), ca_rms(s, valid), ...
            'LineWidth', LW2, 'LineStyle', '--', 'Color', scheme_colors{s}, ...
            'Marker', "square", 'MarkerSize', MS, 'MarkerFaceColor', scheme_colors{s}, ...
            'DisplayName', sprintf('%s (rms)', scheme_labels{s}))
    end
    valid = ~isnan(ca_max(s,:));
    if any(valid)
        loglog(dx_D(valid), ca_max(s, valid), ...
            'LineWidth', LW2, 'LineStyle', '-', 'Color', scheme_colors{s}, ...
            'Marker', 'o', 'MarkerSize', MS, 'MarkerFaceColor', scheme_colors{s}, ...
            'DisplayName', sprintf('%s (max)', scheme_labels{s}))
    end
end

% Reference slopes
loglog(cdx_ref, firstorder, 'LineWidth', 3, 'LineStyle', ':', 'Color', 'k', 'DisplayName', '$O(\Delta x)$')
loglog(cdx_ref, secondorder, 'LineWidth', 3, 'LineStyle', '--', 'Color', 'k', 'DisplayName', '$O(\Delta x^2)$')

set(gca, 'Fontsize', 20)
set(gca, 'fontname', 'Times New Roman')
set(gca, 'LineWidth', 2.0)
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel("$C\Delta = \sqrt{2}\Delta/D$", 'Interpreter', 'latex')
ylabel("$Ca$", 'Interpreter', 'latex')
legend('Location', 'best', 'Interpreter', 'latex')

%% -------- Figure 2: Ca_mean vs CΔ --------
figure('Color','w')
hold on

for s = 1:n_schemes
    valid = ~isnan(ca_mean(s,:));
    if any(valid)
        loglog(dx_D(valid), ca_mean(s, valid), ...
            'LineWidth', LW2, 'LineStyle', '-', 'Color', scheme_colors{s}, ...
            'Marker', '^', 'MarkerSize', MS, 'MarkerFaceColor', scheme_colors{s}, ...
            'DisplayName', scheme_labels{s})
    end
end

loglog(cdx_ref, firstorder, 'LineWidth', 3, 'LineStyle', ':', 'Color', 'k', 'DisplayName', '$O(\Delta x)$')
loglog(cdx_ref, secondorder, 'LineWidth', 3, 'LineStyle', '--', 'Color', 'k', 'DisplayName', '$O(\Delta x^2)$')

set(gca, 'Fontsize', 20)
set(gca, 'fontname', 'Times New Roman')
set(gca, 'LineWidth', 2.0)
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel("$C\Delta = \sqrt{2}\Delta/D$", 'Interpreter', 'latex')
ylabel("$Ca_{\mathrm{mean}}$", 'Interpreter', 'latex')
legend('Location', 'best', 'Interpreter', 'latex')

%% -------- Print summary --------
for s = 1:n_schemes
    fprintf('\n=== %s ===\n', scheme_labels{s});
    disp('Ca rms')
    disp(ca_rms(s,:))
    disp('Ca max')
    disp(ca_max(s,:))
    disp('Ca mean')
    disp(ca_mean(s,:))
end
