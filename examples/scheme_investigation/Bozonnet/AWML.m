% =========================================================================
% Air-Water Mixing Layer Post-Processing
% Benchmarking against Bozonnet et al. (2022) Case B1
% =========================================================================
clear; clc; close all;

%% 1. Physical & Numerical Parameters (Update these for your specific case)
rho_l  = 1000.0;         % Liquid density (kg/m^3)
rho_g  = 1.2;            % Gas density (kg/m^3)
mu_l   = 1e-3;
mu_g   = 1.8e-5;
H_g    = 0.01;
H_l    = 0.01;
U_l    = 0.50;           % Liquid velocity (m/s)
U_g    = 27.0;           % Gas velocity (m/s)
sigma  = 0.072;
gravity= 9.81;

deficit= 1.0;
delta_l= 5e-4;
delta_g= 6*H_g/sqrt(rho_g*U_g*H_g/mu_g);
n      = 2;
dx     = delta_g/2;
ncell  = 240;

fs_sample = 10000;      % Nominal sampling freq (Hz). Actual fs computed from data.
freq_loc=67;    %44,67,89
f_search_min = 10;          % Min frequency for peak search (Hz) - excludes drift
f_search_max = 200;         % Max frequency for peak search (Hz)

% Fit regions (normalized by delta_g) based on the paper's observations
x_exp_min = 2;         % Start of exponential growth region
x_exp_max = 25;         % End of exponential growth region
x_lin_min = 75;         % Start of linear self-similar region
x_lin_max = 175;        % End of linear self-similar region

% type="LVIRA";
type="PLICNET";
% Calculate Theoretical Dimotakis Velocity (U_D)
discard_timelist=[0.5,0.5,0.5];
U_D = (sqrt(rho_l)*U_l + sqrt(rho_g)*U_g) / (sqrt(rho_l) + sqrt(rho_g));
for caseuse=1:3
    
    %% 2. Load Data
    fprintf('Loading data...\n');
    % Assuming CSV format: [time, h(1), h(2), ..., h(end)]
    % Note: Using readmatrix skips headers automatically if present
    data = readmatrix(strcat('results/',type,'/',int2str(caseuse),'/interface_height.csv')); 
    
    t = data(:, 1);
    h = data(:, 2:end);
    [n_time, n_space] = size(h);
    
    % Reconstruct the spatial axis
    x = (0:n_space-1) * dx+dx/2;
    x_norm = x / delta_g;
    
   %% 2b. Convergence Check: Time-Averaging at Specific Downstream Locations

    % 1. Define target locations and find their exact indices in your grid
    target_x = [44, 67, 89];
    loc_indices = zeros(1, length(target_x));
    
    for i = 1:length(target_x)
        [~, loc_indices(i)] = min(abs(x_norm - target_x(i)));
    end
    
    % 2. Extract the time-history of h at these specific locations
    h_time_series = h(:, loc_indices); % Matrix of size [num_timesteps, 3]
    
    % 3. Calculate the Cumulative Moving Average (CMA) for each location
    % cumulative sum divided by the number of steps up to that point
    time_steps = (1:length(t))';
    h_cma = cumsum(h_time_series, 1) ./ time_steps;
    
    % 4. Plotting
    figure('Position', [100, 50, 1000, 400]);
    colors = lines(3); % Get distinct colors for the 3 locations
    
    for i = 1:3
        % Plot instantaneous fluctuating height (light/transparent)
        plot(t, h_time_series(:, i), '-', 'Color', [colors(i,:), 0.3], 'LineWidth', 0.5); hold on;
        
        % Plot the Cumulative Moving Average (Bold)
        plot(t, h_cma(:, i), '-', 'Color', colors(i,:), 'LineWidth', 2, ...
             'DisplayName', sprintf('CMA at x_{norm} = %d', target_x(i)));
    end
    
    grid on;
    xlabel('t (s)'); 
    ylabel('Interface Height, h (m)');
    title('Stationarity Check: Cumulative Time-Average at Fixed Locations');
    legend('show', 'Location', 'best');
    
    % 5. Quantitative Check (Your "Half vs 75%" idea)
    % Let's assume you visually determine from the plot that the first 20% 
    % of your simulation is pure transient startup. We discard it.
    discard_idx = floor(0.50 * length(t)); 
    
    % Take the remaining "stationary" data and split it in half
    stationary_h = h_time_series(discard_idx:end, :);
    mid_idx = floor(size(stationary_h, 1) / 2);
    
    h_first_half = mean(stationary_h(1:mid_idx, :), 1);
    h_second_half = mean(stationary_h(mid_idx+1:end, :), 1);
    h_full_window = mean(stationary_h, 1);
    
    % Print results to console to check convergence
    fprintf('\n--- Stationarity Check (Post-Transient) ---\n');
    for i = 1:3
        fprintf('Location x_norm = %d:\n', target_x(i));
        fprintf('  Mean (First Half):  %.6f\n', h_first_half(i));
        fprintf('  Mean (Second Half): %.6f\n', h_second_half(i));
        fprintf('  Percent Difference: %.3f%%\n', ...
            abs(h_first_half(i) - h_second_half(i)) / abs(h_full_window(i)) * 100);
    end
    discard_time=discard_timelist(caseuse);
    [~, discard_idx] = min(abs(t - discard_time));

    t = t(discard_idx:end);
    h = h(discard_idx:end, :);

    %% 3. Calculate Wave Amplitude A(x)
    fprintf('Calculating wave amplitude envelope...\n');
    % The paper defines amplitude by building a histogram and excluding the 
    % lowest and highest 0.5%. We can do this efficiently using quantiles.
    A = diff(quantile(h, [0.005, 0.995], 1)); 

    %% 4. Robust Frequency Extraction (handles non-uniform dt)
    % -----------------------------------------------------------------
    % Strategy: the simulation uses adaptive CFL-based time stepping, so
    % the raw time series is NOT uniformly sampled. We:
    %   (a) Interpolate to a uniform grid using the actual timestamps.
    %   (b) Apply three independent spectral methods and take the median
    %       as a robust consensus estimate.
    % Target: ~33 Hz (Bozonnet et al. 2022, Case B1)
    % -----------------------------------------------------------------
    fprintf('Extracting dominant frequency (robust multi-method)...\n');

    % --- Probe at x/delta_g = freq_loc (67, per the paper) ---
    [~, probe_idx] = min(abs(x_norm - freq_loc));
    h_raw = h(:, probe_idx);
    t_probe = t;

    % --- Handle non-uniform time stepping ---
    dt_actual = diff(t_probe);
    dt_mean   = mean(dt_actual);
    dt_std    = std(dt_actual);
    fprintf('   Actual mean dt = %.6e s  (std = %.2e, jitter = %.2f%%)\n', ...
            dt_mean, dt_std, 100*dt_std/dt_mean);

    % Interpolate onto a perfectly uniform grid (same number of points)
    t_uniform = linspace(t_probe(1), t_probe(end), length(t_probe));
    h_uniform = interp1(t_probe, h_raw, t_uniform, 'pchip');
    fs_uniform = 1 / mean(diff(t_uniform));  % true sampling freq

    % Detrend: linear detrend + subtract moving average to remove
    % exponential amplitude growth that creates spurious low-freq peaks
    h_detrend = detrend(h_uniform(:));  % remove linear trend first
    % Subtract a slow moving average (window ~ 3 expected wave periods)
    ma_win = min(round(3 * fs_uniform / 33), floor(length(h_detrend)/2));
    if mod(ma_win, 2) == 0; ma_win = ma_win + 1; end  % ensure odd
    h_smooth = movmean(h_detrend, ma_win);
    h_detrend = h_detrend - h_smooth;
    n_pts = length(h_detrend);

    fprintf('   Uniform fs = %.1f Hz,  N = %d pts,  T = %.4f s\n', ...
            fs_uniform, n_pts, t_uniform(end)-t_uniform(1));

    % === Method 1: Welch's PSD (non-parametric, gold standard) ==========
    % Segment length ~ N/2 gives ~5 Hz resolution for 0.4s data while
    % still allowing ~3-4 segments with 75% overlap for variance reduction.
    seg_len     = min(floor(n_pts / 2), 2048);
    noverlap    = floor(seg_len * 0.75);
    nfft_welch  = max(4096, 2^nextpow2(seg_len * 4)); % zero-pad for smooth curve
    [pxx_welch, f_welch] = pwelch(h_detrend, hanning(seg_len), ...
                                  noverlap, nfft_welch, fs_uniform);
    % Restrict peak search to physical frequency band
    band_w = (f_welch >= f_search_min) & (f_welch <= f_search_max);
    [~, idx_w]    = max(pxx_welch .* band_w);
    f_welch_peak  = f_welch(idx_w);

    % === Method 2: Zero-padded FFT (direct, simple) ====================
    % 8x zero-padding interpolates the DFT for a smoother peak location.
    nfft_fft  = 2^nextpow2(n_pts * 8);
    win       = hanning(n_pts);
    H_fft     = abs(fft(h_detrend .* win, nfft_fft));
    f_fft_ax  = (0:nfft_fft-1) * fs_uniform / nfft_fft;
    half_n    = floor(nfft_fft / 2);
    H_fft     = H_fft(1:half_n);
    f_fft_ax  = f_fft_ax(1:half_n);
    band_f = (f_fft_ax >= f_search_min) & (f_fft_ax <= f_search_max);
    [~, idx_f]   = max(H_fft' .* band_f);
    f_fft_peak   = f_fft_ax(idx_f);

    % === Consensus: average of Welch + FFT ===============================
    f_peak   = 0.5 * (f_welch_peak + f_fft_peak);
    f_spread = abs(f_welch_peak - f_fft_peak);

    fprintf('   Welch  peak: %6.2f Hz\n', f_welch_peak);
    fprintf('   FFT    peak: %6.2f Hz\n', f_fft_peak);
    fprintf('   --> Consensus (avg): %.2f Hz  (spread: %.2f Hz)\n', f_peak, f_spread);
    if f_spread > 3
        fprintf('   *** WARNING: Welch and FFT disagree by >3 Hz — inspect spectra ***\n');
    end
    %% 5. Exponential Spatial Growth Rate (k_i)
    % A = A_0 * exp(k_i * x)  ==>  ln(A) = k_i * x + ln(A_0)
    idx_exp = (x_norm >= x_exp_min) & (x_norm <= x_exp_max);

    % Linear fit on semi-log data
    P_exp = polyfit(x(idx_exp), log(A(idx_exp)), 1);
    k_i   = P_exp(1);
    A_0   = exp(P_exp(2));

    fprintf('--> Exponential Growth Rate (k_i): %.2f m^-1\n', k_i);

    %% 6. Self-Similar Linear Growth (C_0)
    % A = C_0 * sqrt(rho_g/rho_l) * U_g * (x / U_D)
    % dA/dx = C_0 * sqrt(rho_g/rho_l) * (U_g / U_D)
    idx_lin = (x_norm >= x_lin_min) & (x_norm <= x_lin_max);

    % Linear fit on linear data
    P_lin = polyfit(x(idx_lin), A(idx_lin), 1);
    slope_lin = P_lin(1);

    % Solve for C_0
    C_0 = slope_lin * (U_D / U_g) * sqrt(rho_l / rho_g);

    fprintf('--> Self-Similar Growth Constant (C_0): %.4f\n', C_0);

    %% 7. Interfacial Wave Speed (c)
    fprintf('Estimating numerical wave speed...\n');
    % Cross-correlate two adjacent probes in the linear region to find phase lag
    [~, idx1] = min(abs(x_norm - 80));
    [~, idx2] = min(abs(x_norm - 85)); % Probe 5*delta_g downstream

    [r, lags] = xcorr(h(:, idx2) - mean(h(:, idx2)), h(:, idx1) - mean(h(:, idx1)));
    [~, max_lag_idx] = max(r);
    time_lag = lags(max_lag_idx) * (1 / fs_uniform);

    dist = x(idx2) - x(idx1);
    c_num = dist / time_lag;

    fprintf('--> Theoretical U_D: %.2f m/s\n', U_D);
    fprintf('--> Numerical Wave Speed: %.2f m/s\n', c_num);

    %% 8. Visualizations
    figure('Position', [100, 100, 1200, 800]);

    % Plot 1: Amplitude vs. Downstream Distance (Log-Linear)
    subplot(2, 2, [1, 2]);
    semilogy(x_norm, A, 'b-', 'LineWidth', 1.5); hold on;
    semilogy(x_norm(idx_exp), exp(polyval(P_exp, x(idx_exp))), 'r--', 'LineWidth', 2);
    % Plot the linear fit model mapped back to the log plot
    semilogy(x_norm(idx_lin), polyval(P_lin, x(idx_lin)), 'k--', 'LineWidth', 2);
    grid on;
    xlabel('x / \delta_g');
    ylabel('A / \delta_g (Normalized by \delta_g for viewing)');
    title('Spatial Evolution of Wave Amplitude');
    legend('Simulation Data', 'Exponential Fit', 'Linear Fit', 'Location', 'SouthEast');
    xlim([0 max(x_norm)]);

    % Plot 2: Frequency Spectra (Welch + FFT)
    subplot(2, 2, 3);
    plot(f_welch, pxx_welch / max(pxx_welch), 'b-', 'LineWidth', 1.5); hold on;
    plot(f_fft_ax, H_fft / max(H_fft), 'k--', 'LineWidth', 1.0);
    xline(f_peak, 'm-', 'LineWidth', 1.5);  % consensus line
    grid on;
    xlim([0 150]);
    xlabel('f (Hz)');
    ylabel('Normalized Power');
    legend('Welch', 'FFT', sprintf('Consensus %.1f Hz', f_peak), ...
           'Location', 'NorthEast');
    title(sprintf('Frequency Spectra at x/\\delta_g = %d', freq_loc));

    % Plot 3: Hovmöller (Space-Time) Diagram
    subplot(2, 2, 4);
    % To keep the plot manageable, only plot the last 0.2 seconds
    t_plot_idx = t > (max(t) - 0.5); 
    imagesc(x_norm, t(t_plot_idx), h(t_plot_idx, :));
    colormap(gray); % Matches Figure 11 in the paper
    set(gca, 'YDir', 'normal');
    xlabel('x / \delta_g');
    ylabel('t (s)');
    title('Hovmöller Diagram (Interface Height)');
    % Overlay U_D trajectory
    hold on;
    x_traj = linspace(50, 150, 100); % Plot line between x/dg = 50 and 150
    t_start = t(find(t_plot_idx, 1)) + 0.05; % Arbitrary start time for the line
    t_traj = t_start + (x_traj * delta_g) / U_D;
    plot(x_traj, t_traj, 'm-', 'LineWidth', 2);
    legend('U_D Trajectory');
end