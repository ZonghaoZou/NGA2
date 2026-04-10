% 1. Initialization and Parameters
filename = '1/monitor/simulation';

% Define your boundary layer thickness (Update this if your dg changes!)
dg = 2.5e-3; 
Ug=10;
threshold_multiplier = 0.3; % 10% of dg
amp_threshold = threshold_multiplier * dg;
mksize=2;
skipcell=1200;
% Read the file, skipping the header lines
data = readmatrix(filename, 'NumHeaderLines', 2);

% Extract Time (Column 2) and Amplitude (Column 15)
% Note: Using (:, 2) assuming readmatrix stripped all text headers
raw_time = data(:, 2);
raw_amp  = data(:, 15);
raw_amp_height  = data(:, 16);

% Filter out any non-positive amplitudes to prevent log(0) = -Inf
valid_idx = raw_amp > 0;
time = raw_time(valid_idx);
amp  = raw_amp(valid_idx);
N = length(time);

if isempty(time)
    error('No valid data points with amplitude > 0 found.');
end

log_beta = log(amp);
log_threshold = log(amp_threshold);

% =========================================================
% 2. Apply Filters (The New OLS Approach)
% =========================================================

% Filter 1: Start-up skip (skip the first 10 valid points)
start_idx = skipcell; 
if N < skipcell-1
    start_idx = 2; % Fallback for very short simulation runs
end

% Filter 2: Non-linear cutoff (stop when amplitude > 0.1 * dg)
end_idx = N;
for i = start_idx:N
    if amp(i) > amp_threshold || raw_amp_height (i) > amp_threshold
        end_idx = i - 1;
        break;
    end
end

if end_idx < start_idx
    error('Wave exceeded 0.1 dg before the start-up transient finished.');
end

% =========================================================
% 3. Extract the Pure Linear Band & Calculate OLS Growth Rate
% =========================================================
t_linear = time(start_idx:end_idx);
logb_linear = log_beta(start_idx:end_idx);

mean_t = mean(t_linear);
mean_logb = mean(logb_linear);

% Analytical OLS Slope formula
num = sum((t_linear - mean_t) .* (logb_linear - mean_logb));
den = sum((t_linear - mean_t).^2);

if den > 0
    growth_rate = num / den;
    intercept = mean_logb - growth_rate * mean_t;
else
    growth_rate = 0;
    intercept = 0;
end

% Print Terminal Diagnostics
fprintf('--- OLS Fit Diagnostics ---\n');
fprintf('Total data points: %d\n', N);
fprintf('Start-up skipped points: %d\n', start_idx - 1);
fprintf('Linear regime ended at point: %d (Amplitude = %.4e)\n', end_idx, amp(end_idx));
fprintf('Points used for OLS fit: %d\n', length(t_linear));
fprintf('Extracted Growth Rate (Slope): %e\n', growth_rate*dg/Ug);
fprintf('---------------------------\n');

% =========================================================
% 4. Visualization for Debugging
% =========================================================
figure('Name', 'Linear Stability Regime Debugger', 'Color', 'w', 'Position', [100, 100, 900, 600]);

% Plot the raw, unfiltered data as a faint gray background line
plot(time, log_beta, '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 1, 'HandleVisibility', 'off'); hold on;

% Highlight Phase 1: The ignored start-up noise (Yellow)
if start_idx > 1
    plot(time(1:start_idx-1), log_beta(1:start_idx-1), 'o', ...
        'Color', [0.9290 0.6940 0.1250], 'MarkerSize', mksize, 'MarkerFaceColor', [0.9290 0.6940 0.1250], ...
        'DisplayName', 'Start-Up Noise (Ignored)');
end

% Highlight Phase 3: The non-linear tail (Red)
if end_idx < N
    plot(time(end_idx+1:end), log_beta(end_idx+1:end), 'ro', ...
        'MarkerFaceColor', 'r', 'MarkerSize', mksize, 'DisplayName', 'Non-Linear Tail (Ignored)');
end
% Highlight Phase 2: The pure linear band used for the math (Blue)
plot(t_linear, logb_linear, 'bo', 'MarkerFaceColor', 'b', ...
    'MarkerSize', mksize, 'DisplayName', 'Pure Linear Band (Fitted)');

% Plot the fitted regression line extending across the whole plot
% This lets you visually see where the non-linear tail starts to deviate from the perfect line
plot(time, intercept + growth_rate * time, 'b--', 'LineWidth', 1.0, ...
    'DisplayName', 'OLS Fitted Line');

% Draw the hard amplitude threshold
yline(log_threshold, 'g-', sprintf('0.1 \\delta_g Threshold (amp = %.2e)', amp_threshold), ...
    'LabelHorizontalAlignment', 'left', 'LineWidth', 2, 'DisplayName', 'Cutoff Threshold');

% Formatting
xlabel('Time [s]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('log(Amplitude)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Theoretical Growth Rate Extraction (Slope = %.5e)', growth_rate), 'FontSize', 14);
legend('Location', 'northwest', 'FontSize', 11);
grid on;
ax = gca;
ax.FontSize = 11;