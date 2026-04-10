%% getting from KE conserving case
clear; close all
format longG

meshsizelist=[64,128,256];
caselist=[1,2,3];

% Initialize Data Storage
all_data = cell(length(meshsizelist), length(caselist));

% --- READ SIMULATION DATA ---
for i=1:length(meshsizelist)
    mesh=int2str(meshsizelist(i));
    for j=1:length(caselist)
        folder=int2str(caselist(j));
        % Using fullfile for safe path construction
        filename=fullfile("result", folder, "output_" + mesh + ".csv");
        
        if isfile(filename)
            startRow = 3;
            formatSpec = '%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%f%[^\n\r]';
            fileID = fopen(filename,'r');
            dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
            fclose(fileID);       
            all_data{i,j} = [dataArray{1:end-1}];
        else
            % warning('File not found: %s', filename);
        end
    end 
end 

%% 1. Setup Constants and Normalization
H=0.1; sigma=0.45; rho1=900; rho2=1000; mu1=0.1; mu2=0.001; g=9.81;

% Derived Scales
Ug     = sqrt(g*H/2)*(rho2-rho1)/rho1;
t_c    = H/(2*Ug);
Ep1_c  = 0.1035; % Hardcoded from user snippet
Ep2_c  = 0.3755; % Hardcoded from user snippet
Ek1_c  = rho1*Ug^2*H^3/16; % Updated to include H^3
Ek2_c  = rho2*Ug^2*H^3/16; % Updated to include H^3
EN1_c  = 0.0733;
EN2_c  = 1.3759;

%% 2. Setup Visualization Styles
% Colors for Cases (j index) - Swapped Convention
my_colors = [71,135,224; ...   % Blue
             218,62,32; ...    % Red
             82,147,47; ...    % Green
             241,178,56] / 255; % Yellow
% Line Styles for Mesh Sizes (i index) - Swapped Convention
my_styles = ["-", "--", ":", "-."];

%% 3. Initialize Figures (8 Total)
f1 = figure('Name', 'Liquid Kinetic Energy (Log-Log)'); hold on; grid on;
f2 = figure('Name', 'Gas Kinetic Energy (Semi-Log)');   hold on; grid on;
f3 = figure('Name', 'Liquid Potential Energy');         hold on; grid on;
f4 = figure('Name', 'Gas Potential Energy');            hold on; grid on;
f5 = figure('Name', 'Liquid Enstrophy (Linear)');       hold on; grid on;
f6 = figure('Name', 'Gas Enstrophy');                   hold on; grid on;
f7 = figure('Name', 'Total Mechanical Energy');         hold on; grid on;
f8 = figure('Name', 'Liquid Enstrophy (Log-Log)');      hold on; grid on;

%% 4. Loop and Plot Simulation Data
for i = 1:length(meshsizelist)
    for j = 1:length(caselist)
        
        data = all_data{i,j};
        if isempty(data), continue; end
        
        % --- Select Style (SWAPPED) ---
        % Case (j) determines Color
        c_idx = mod(j-1, 4) + 1; 
        
        % Mesh (i) determines Line Style
        s_idx = mod(i-1, 4) + 1; 
        
        curr_color = my_colors(c_idx, :);
        curr_style = my_styles(s_idx);
        disp_name  = sprintf('M:%d C:%d', meshsizelist(i), caselist(j));
        
        % --- Extract and Normalize ---
        t = data(:,1) / t_c;
        
        ke_l = data(:,4) / Ek2_c;
        ke_g = data(:,5) / Ek1_c;
        pe_l = data(:,15) / Ep2_c;
        pe_g = data(:,16) / Ep1_c;
        en_l = data(:,17) / EN2_c;
        en_g = data(:,18) / EN1_c;
        
        % --- Total Energy Calculation (Exact User Logic) ---
        % Sum raw columns: KE_liq + KE_gas + PE_liq + PE_gas
        % tot = data(:,4) + data(:,5) + data(:,15) + data(:,16);
        tot = data(:,3)+ data(:,15) + data(:,16);
        tot0 = tot(1);
        % tot = tot /tot0;
        % % Normalize drift by potential energy scale
        tot = -(tot -tot0) ./ (Ep1_c + Ep2_c-tot0)+1;
        
        % --- Plotting ---
        figure(f1); plot(t, ke_l, 'Color', curr_color, 'LineStyle', curr_style, 'LineWidth', 1.5, 'DisplayName', disp_name);
        figure(f2); plot(t, ke_g, 'Color', curr_color, 'LineStyle', curr_style, 'LineWidth', 1.5, 'DisplayName', disp_name);
        figure(f3); plot(t, pe_l, 'Color', curr_color, 'LineStyle', curr_style, 'LineWidth', 1.5, 'DisplayName', disp_name);
        figure(f4); plot(t, pe_g, 'Color', curr_color, 'LineStyle', curr_style, 'LineWidth', 1.5, 'DisplayName', disp_name);
        figure(f5); plot(t, en_l, 'Color', curr_color, 'LineStyle', curr_style, 'LineWidth', 1.5, 'DisplayName', disp_name);
        figure(f6); plot(t, en_g, 'Color', curr_color, 'LineStyle', curr_style, 'LineWidth', 1.5, 'DisplayName', disp_name);
        figure(f7); plot(t, tot,  'Color', curr_color, 'LineStyle', curr_style, 'LineWidth', 1.5, 'DisplayName', disp_name);
        figure(f8); plot(t, en_l, 'Color', curr_color, 'LineStyle', curr_style, 'LineWidth', 1.5, 'DisplayName', disp_name);
        
    end
end

%% 4.5. Loop and Plot Reference Data (DyJeatresult)
% Mapping: fluid1 = Gas, fluid2 = Liquid
ref_dir = "DyJeatresult/";

for j = 1:length(caselist)
    caseNum = caselist(j);
    
    % --- RESTRICTION: Only plot Case 1 from Experiment ---
    if caseNum ~= 1
        continue;
    end
    
    % Construct Filenames
    fn_ke_g = fullfile(ref_dir, sprintf("KEfluid1_case%d.csv", caseNum)); 
    fn_ke_l = fullfile(ref_dir, sprintf("KEfluid2_case%d.csv", caseNum)); 
    fn_pe_g = fullfile(ref_dir, sprintf("PEfluid1_case%d.csv", caseNum)); 
    fn_pe_l = fullfile(ref_dir, sprintf("PEfluid2_case%d.csv", caseNum)); 
    fn_en_g = fullfile(ref_dir, sprintf("ENfluid1_case%d.csv", caseNum)); 
    fn_en_l = fullfile(ref_dir, sprintf("ENfluid2_case%d.csv", caseNum)); 
    fn_tot  = fullfile(ref_dir, sprintf("TotalEnergy_case%d.csv", caseNum)); 
    
    
    % --- STYLE UPDATE: Solid Line, Black, Fixed Legend "DyJeat512" ---
    ref_style = '-';      % Solid line
    ref_color = 'k';      % Black color
    ref_marker = 'none'; 
    ref_legend = "DyJeat512"; 
    
    % Helper function to read and plot
    plot_ref = @(fname, fig_handle) ...
        plotting_helper(fname, fig_handle, ref_legend, ref_color, ref_style, ref_marker);
    
    % --- Plot Individual Quantities ---
    plot_ref(fn_ke_g, f2); % Gas KE
    plot_ref(fn_ke_l, f1); % Liq KE
    plot_ref(fn_pe_g, f4); % Gas PE
    plot_ref(fn_pe_l, f3); % Liq PE
    plot_ref(fn_en_g, f6); % Gas Ens
    plot_ref(fn_en_l, f5); % Liq Ens
    plot_ref(fn_en_l, f8); % Liq Ens (Log-Log)
    
    % For Total Energy, plot the reference file directly.
    % (Assumes file is already scaled to drift as per instructions)
    plot_ref(fn_tot,  f7); 
end

%% 5. Final Formatting
figs = [f1, f2, f3, f4, f5, f6, f7, f8];
titles = [...
    "Liquid KE (Log-Log)", ...
    "Gas KE (Semi-Log)", ...
    "Liquid PE", ...
    "Gas PE", ...
    "Liquid Enstrophy (Linear)", ...
    "Gas Enstrophy", ...
    "Total Mechanical Energy", ...
    "Liquid Enstrophy (Log-Log)"];
y_labels = [...
    "KE / Ek2_c", ...
    "KE / Ek1_c", ...
    "PE / Ep2_c", ...
    "PE / Ep1_c", ...
    "Ens / EN2_c", ...
    "Ens / EN1_c", ...
    "Energy Sum", ...
    "Ens / EN2_c"];
for k = 1:8
    figure(figs(k));
    
    % --- 1. Set Time Scale (X-Axis) ---
    % Log-Log plots: f1 (Liq KE) and f8 (Liq Enstrophy Log-Log)
    if k == 1 || k == 8
        set(gca, 'XScale', 'log');
        xlim([0.1, 20]); % Log scale range
    elseif k==7
        set(gca, 'XScale', 'linear');
        xlim([0, 10]);   % Linear scale range
    else
        set(gca, 'XScale', 'linear');
        xlim([0, 20]);   % Linear scale range
    end
    
    % --- 2. Set Value Scale (Y-Axis) ---
    % Log Y plots: f1 (Liq KE), f2 (Gas KE), f8 (Liq Enstrophy Log-Log)
    if k == 1 || k == 2 || k == 8
        set(gca, 'YScale', 'log');
    else
        set(gca, 'YScale', 'linear');
    end
    
    % --- 3. Labels ---
    xlabel('Time (t/t_c)');
    ylabel(y_labels(k));
    title(titles(k));
    legend('show', 'Location', 'bestoutside');
    grid on;
end

%% --- Helper Function ---
function data = plotting_helper(filename, fig_h, disp_name, col, sty, mrk)
    data = [];
    if isfile(filename)
        try
            % Read CSV (Assumes Col 1 = Time, Col 2 = Value)
            raw = readmatrix(filename);
            if size(raw,2) >= 2
                data = raw;
                figure(fig_h);
                plot(raw(:,1), raw(:,2), 'Color', col, 'LineStyle', sty, ...
                     'Marker', mrk, 'LineWidth', 1.5, 'DisplayName', disp_name);
            end
        catch
            warning('Could not read %s', filename);
        end
    end
end