% Enhanced script to compute derivatives for multiple fields
% This script computes d(heat_release)/d(field) for multiple fields
% where numerator is always heat release and denominator changes for each field

clear; clc;
close all;
addpath('~/MATLAB');
addpath('~/satyam_files/CH4_jet_PF/2025_Runs/derivatives_files/functions');

%% User Input - Work Directory
work_dir = '/work/home/satyam/satyam_files/CH4_jet_PF/2025_Runs/derivatives_files';

%% Configuration
data_dir = '/work/home/satyam/satyam_files/CH4_jet_PF/2025_Runs/c_cond_stats/C_cond_fields_800_10D_n';

D = 2e-3;  % Diameter scale
save_results_flag = false;
z_ref_CH4 = 8.5*D;z_ref = 8*D;
% Create output directories
sensitivities_dir = fullfile(work_dir, 'sensitivities_10D');
figures_dir = fullfile(sensitivities_dir, 'figures');

if ~exist(sensitivities_dir, 'dir')
    mkdir(sensitivities_dir);
    fprintf('Created sensitivities directory: %s\n', sensitivities_dir);
end

if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
    fprintf('Created figures directory: %s\n', figures_dir);
end

% Define fields to process (heat release is always the numerator)
field_configs = {
    % {field_name, smooth_suffix, latex_label, plot_title, derivative_label}
    'Temperature', '_smooth_bc_smooth', 'T', '$T$', '$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial T}$';
    'density', '_smooth_bc_smooth', '\rho', '$\rho$', '$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial \rho}$';
%     'CH2O', '_smooth', 'Y_{CH_2O}', '$Y_{CH2O}$', '$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial Y_{CH2O}}$';
%     'CH3', '_smooth', 'Y_{CH_3}', '$Y_{CH3}$', '$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial Y_{CH3}}$';
    'CH4', '_smooth_bc_smooth', 'Y_{CH_4}', '$Y_{CH4}$', '$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial Y_{CH4}}$';
%     'CO', '_smooth', 'Y_{CO}', '$Y_{CO}$', '$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial Y_{CO}}$';
    'CO2', '_smooth_bc_smooth', 'Y_{CO_2}', '$Y_{CO2}$', '$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial Y_{CO2}}$';
%     'H', '_smooth', 'Y_{H}', '$Y_{H}$', '$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial Y_{H}}$';
%     'H2', '_smooth', 'Y_{H_2}', '$Y_{H2}$', '$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial Y_{H2}}$';
    'H2O', '_smooth_bc_smooth', 'Y_{H_2O}', '$Y_{H2O}$', '$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial Y_{H2O}}$';
%     'HO2', '_smooth', 'Y_{HO_2}', '$Y_{HO2}$', '$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial Y_{HO2}}$';
%     'N2', '_smooth', 'Y_{N_2}', '$Y_{N2}$', '$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial Y_{N2}}$';
%     'O', '_smooth', 'Y_{O}', '$Y_{O}$', '$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial Y_{O}}$';
    'O2', '_smooth_bc_smooth', 'Y_{O_2}', '$Y_{O2}$', '$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial Y_{O2}}$';
%     'OH', '_smooth', 'Y_{OH}', '$Y_{OH}$', '$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial Y_{OH}}$';
};
% field_configs = {
%     % {field_name, smooth_suffix, latex_label, plot_title, short_name}
%     'Temperature', '_smooth_bc_smooth', 'T', '$T$', 'T';
%     'density', '_smooth_bc_smooth', '\rho', '$\rho$', 'rho';
%     'CH4', '_smooth_bc_smooth', 'Y_{CH_4}', '$Y_{CH4}$', 'CH4';
%     'O2', '_smooth_bc_smooth', 'Y_{O_2}', '$Y_{O2}$', 'O2';
%     'CO2', '_smooth_bc_smooth', 'Y_{CO_2}', '$Y_{CO2}$', 'CO2';
%     'H2O', '_smooth_bc_smooth', 'Y_{H_2O}', '$Y_{H2O}$', 'H2O';
% };

%% Load heat release data (numerator for all calculations)
fprintf('Loading heat release data...\n');
heat_release_field = 'Heatrelease_smooth';
try
    heatRelease = load(sprintf('%s/%s.mat', data_dir, heat_release_field));
    fprintf('? Heat release data loaded successfully\n');
catch ME
    fprintf('? Error loading heat release data: %s\n', ME.message);
    return;
end
[C_MAT, Z_MAT] = get_CZ_coord_data(data_dir);
% Compute heat release derivative (numerator for all calculations)
fprintf('Computing heat release derivative...\n');
dq_dc = compute_dfdr(heatRelease.DF, C_MAT);
% dq_dc = apply_downstream_replacement(dq_dc,Z_MAT,z_ref);

%% Process each field
fprintf('\nProcessing %d fields...\n', size(field_configs, 1));
results = struct();
successful_fields = {};
successful_fields_fig = {};
fieldname = 'hrr';
results.(fieldname) = struct();
results.(fieldname).field_data = heatRelease.DF;
results.(fieldname).field_derivative = dq_dc;
results.(fieldname).C_MAT = C_MAT;
results.(fieldname).Z_MAT = Z_MAT;
results.(fieldname).latex_label = '\dot{\omega}_{T}';
results.(fieldname).derivative_label = '\dot{\omega}_{T}';

for i = 1:size(field_configs, 1)
    field_name = field_configs{i, 1};
    smooth_suffix = field_configs{i, 2};
    latex_label = field_configs{i, 3};
    plot_title = field_configs{i, 4};
    derivative_label = field_configs{i, 5};
    
    fprintf('\n--- Processing Field %d/%d: %s ---\n', i, size(field_configs, 1), field_name);
    
    % Construct field filename
    field_filename = sprintf('%s%s', field_name, smooth_suffix);
    field_path = sprintf('%s/%s.mat', data_dir, field_filename);
    
    % Check if field exists
    if ~exist(field_path, 'file')
        fprintf('? Field file not found: %s\n', field_path);
        % Try without smooth suffix
        field_filename_alt = field_name;
        field_path_alt = sprintf('%s/%s.mat', data_dir, field_filename_alt);
        if exist(field_path_alt, 'file')
            fprintf('  Found alternative: %s\n', field_filename_alt);
            field_path = field_path_alt;
            field_filename = field_filename_alt;
        else
            fprintf('  Skipping field: %s\n', field_name);
            continue;
        end
    end
    
    try
        % Load field data
        fprintf('  Loading field data...\n');
        field_data = load(field_path);
        
        % Compute field derivative
        fprintf('  Computing derivative d%s/dc...\n', field_name);
        df_dc = compute_dfdr(field_data.DF, C_MAT);
        if strcmp(field_name,'CH4')
            df_dc = apply_downstream_replacement(df_dc,Z_MAT,z_ref_CH4);
        else
            df_dc = apply_downstream_replacement(df_dc,Z_MAT,z_ref);
        end
        if strcmp(field_name,'CO2')
                df_dc = set_boundary_to_zero(df_dc, 'Boundaries', {'right'});
        end
        % Handle zero values in denominator
        zero_idx = find(abs(df_dc) < 1e-16);
        if ~isempty(zero_idx)
            fprintf('  Handling %d zero values in denominator...\n', length(zero_idx));
            df_dc(zero_idx) = 1e-8;
        end
        
        % Compute final derivative: d(heat_release)/d(field)
        fprintf('  Computing final derivative: d(heat_release)/d(%s)...\n', field_name);
        derivative_result = dq_dc ./ df_dc;
        if strcmp(field_name,'CH4')
            derivative_result = set_sensitivities_bc(derivative_result, Z_MAT, z_ref_CH4);
        else
            derivative_result = set_sensitivities_bc(derivative_result, Z_MAT, z_ref);
        end
%         if strcmp(field_name,'CH4') 
%             derivative_result = thresholding_and_smoothing(derivative_result,5,0.01);
%         elseif strcmp(field_name,'N2')
%             derivative_result = thresholding_and_smoothing(derivative_result,3,0.001);
%         end
        % Store results
        results.(field_name) = struct();
        results.(field_name).field_data = field_data.DF;
        results.(field_name).field_derivative = df_dc;
        results.(field_name).final_derivative = derivative_result;
        results.(field_name).C_MAT = C_MAT;
        results.(field_name).Z_MAT = Z_MAT;
        results.(field_name).latex_label = latex_label;
        results.(field_name).plot_title = plot_title;
        results.(field_name).derivative_label = derivative_label;
        
        successful_fields_fig{end+1} = latex_label;
        successful_fields{end+1} = field_name;
        fprintf('  ? Successfully processed %s\n', field_name);
        
    catch ME
        fprintf('  ? Error processing %s: %s\n', field_name, ME.message);
        continue;
    end
end

%% Generate plots for all successful fields
if ~isempty(successful_fields)
    fprintf('\n=== Generating Plots ===\n');
    fprintf('Successfully processed %d fields: %s\n', length(successful_fields), strjoin(successful_fields, ', '));
    
    % Create plots for each field
    for i = 1:length(successful_fields)
        field_name = successful_fields{i};
        f_name = successful_fields_fig{i};
        plot_field_results(results.(field_name), results.hrr, f_name, D, i, figures_dir,C_MAT,Z_MAT);
        plot_surface_field_results(results.(field_name), results.hrr, f_name, D, i, figures_dir, C_MAT, Z_MAT)
    end
    
    % Create summary comparison plot
%     create_summary_plot(results, successful_fields,successful_fields_fig, heatRelease, D, figures_dir);
    
    % Save results (sensitivities and data)
    if save_results_flag
        fprintf("Saving Data \n");
        save_results(results, successful_fields, sensitivities_dir, work_dir);
    end
    
else
    fprintf('\n? No fields were successfully processed.\n');
end

fprintf('\n=== Processing Complete ===\n');

%% Function Definitions

function dfield_dc = compute_dfdr(field_data, C_MAT)
    % Compute derivative of field with respect to c using central differences
    % This function approximates the missing compute_dfdr function
    
    [rows, cols] = size(field_data);
    dfield_dc = zeros(size(field_data));
    
    % Get c spacing (assuming uniform grid)
    dc = C_MAT(1, 2) - C_MAT(1, 1);
    
    % Central differences for interior points
    for i = 1:rows
        for j = 2:cols-1
            dfield_dc(i, j) = (field_data(i, j+1) - field_data(i, j-1)) / (2 * dc);
        end
    end
    
    % Forward difference for first column
    for i = 1:rows
        dfield_dc(i, 1) = (field_data(i, 2) - field_data(i, 1)) / dc;
    end
    
    % Backward difference for last column
    for i = 1:rows
        dfield_dc(i, end) = (field_data(i, end) - field_data(i, end-1)) / dc;
    end
end

function plot_field_results(field_result, heatRelease, field_name, D, figure_offset, figures_dir, C_MAT, Z_MAT)
    % Plot original field, its derivative, and final result in a single
    % maximized figure
    
    base_figure = 100 + figure_offset * 10;
    
    % Create maximized figure with 3 subplots
    figure(base_figure);
    set(gcf, 'WindowState', 'maximized');  % Maximize the figure
    set(gcf, 'Position', get(0, 'Screensize'));  % Alternative maximization method
    
    % Subplot 1:heat release derivative
    subplot(1, 4, 1);
        myutils.plot_contourf(gcf, C_MAT, Z_MAT/D, heatRelease.field_derivative, ...
            '$c$', '$z/D$', sprintf('$\\frac{\\partial \\langle %s|c\\rangle}{\\partial c}$', heatRelease.latex_label));

    % Subplot 2: Original field
    subplot(1, 4, 2);

        myutils.plot_contourf(gcf, field_result.C_MAT, field_result.Z_MAT/D, field_result.field_data, ...
            '$c$', '$z/D$', sprintf('$\\langle %s|c\\rangle$ ', field_name));
    
    % Subplot 2: Field derivative
    subplot(1, 4, 3);
        myutils.plot_contourf(gcf, field_result.C_MAT, field_result.Z_MAT/D, field_result.field_derivative, ...
            '$c$', '$z/D$', sprintf('$\\frac{\\partial \\langle %s|c\\rangle}{\\partial c}$ ', field_name));
    
    % Subplot 3: Final derivative result (sensitivity)
    subplot(1, 4, 4);
        myutils.plot_contourf(gcf, field_result.C_MAT, field_result.Z_MAT/D, field_result.final_derivative, ...
            '$c$', '$z/D$', sprintf('$\\frac{\\partial \\langle \\dot{\\omega}_{T}|c\\rangle}{\\partial %s}$ ', field_name));
    
    % Add overall title
    sgtitle(sprintf('%s Field Analysis', field_name), 'FontSize', 16, 'FontWeight', 'bold');
    
    % Adjust subplot spacing
    set(gcf, 'Units', 'normalized');
    
    % Save the combined figure
    try
        % Save as PNG and FIG
        saveas(gcf, fullfile(figures_dir, sprintf('%s_combined_analysis.png', field_name)));
        saveas(gcf, fullfile(figures_dir, sprintf('%s_combined_analysis.fig', field_name)));
        
        % Also save individual plots for compatibility
        % Extract and save subplot 1 (original field)
%         figure(base_figure + 100);  % Use different figure numbers to avoid conflicts
%         try
%             myutils.plot_contourf(gcf, field_result.C_MAT, field_result.Z_MAT/D, field_result.field_data, ...
%                 '$c$', '$z/D$', sprintf('$\\langle %s|c\\rangle$', field_name));
%         catch
%             contourf(field_result.C_MAT, field_result.Z_MAT/D, field_result.field_data, 20);
%             colorbar;
%             xlabel('$c$', 'Interpreter', 'latex');
%             ylabel('$z/D$', 'Interpreter', 'latex');
%             title(sprintf('$\\langle %s|c\\rangle$', field_name), 'Interpreter', 'latex');
%         end
%         saveas(gcf, fullfile(figures_dir, sprintf('%s_original_field.png', field_name)));
%         close(gcf);
%         
%         % Extract and save subplot 2 (field derivative)
%         figure(base_figure + 200);
%         try
%             myutils.plot_contourf(gcf, field_result.C_MAT, field_result.Z_MAT/D, field_result.field_derivative, ...
%                 '$c$', '$z/D$', sprintf('$\\frac{\\partial \\langle %s|c\\rangle}{\\partial c}$', field_name));
%         catch
%             contourf(field_result.C_MAT, field_result.Z_MAT/D, field_result.field_derivative, 20);
%             colorbar;
%             xlabel('$c$', 'Interpreter', 'latex');
%             ylabel('$z/D$', 'Interpreter', 'latex');
%             title(sprintf('$\\frac{\\partial \\langle %s|c\\rangle}{\\partial c}$', field_name), 'Interpreter', 'latex');
%         end
%         saveas(gcf, fullfile(figures_dir, sprintf('%s_field_derivative.png', field_name)));
%         close(gcf);
%         
%         % Extract and save subplot 3 (sensitivity)
%         figure(base_figure + 300);
%         try
%             myutils.plot_contourf(gcf, field_result.C_MAT, field_result.Z_MAT/D, field_result.final_derivative, ...
%                 '$c$', '$z/D$', sprintf('$\\frac{\\partial \\langle \\dot{\\omega}_{T}|c\\rangle}{\\partial %s}$', field_name));
%         catch
%             contourf(field_result.C_MAT, field_result.Z_MAT/D, field_result.final_derivative, 20);
%             colorbar;
%             xlabel('$c$', 'Interpreter', 'latex');
%             ylabel('$z/D$', 'Interpreter', 'latex');
%             title(sprintf('$\\frac{\\partial \\langle \\dot{\\omega}_{T}|c\\rangle}{\\partial %s}$', field_name), 'Interpreter', 'latex');
%         end
%         saveas(gcf, fullfile(figures_dir, sprintf('%s_sensitivity.png', field_name)));
%         close(gcf);
        
        fprintf('  ? Combined and individual figures saved for %s\n', field_name);
    catch ME
        fprintf('  ? Error saving figures for %s: %s\n', field_name, ME.message);
    end
    
    fprintf('  ? Generated combined plot for %s (figure %d)\n', field_name, base_figure);
end

function create_summary_plot(results, successful_fields, successful_fields_fig, heatRelease, D, figures_dir)
    % Create a summary plot comparing all derivative results
    
    if length(successful_fields) < 2
        fprintf('  Need at least 2 fields for summary plot\n');
        return;
    end
    
    figure("WindowState","maximized");

    n_fields = length(successful_fields);
    subplot_rows = ceil(sqrt(n_fields));
    subplot_cols = ceil(n_fields / subplot_rows);
    
    for i = 1:n_fields
        field_name = successful_fields{i};
        field_result = results.(field_name);
        field_name_fig = successful_fields_fig{i};
        
        subplot(subplot_rows, subplot_cols, i);
            myutils.plot_contourf(gcf, field_result.C_MAT, field_result.Z_MAT/D, field_result.final_derivative, ...
                '$r/D$', '$c$',sprintf('$\\frac{\\partial \\dot{\\omega}_T}{\\partial %s}$', field_name_fig));
            pbaspect([1 2 1]);
    end
    
    sgtitle('Summary: Heat Release Derivatives w.r.t. Different Fields', ...
        'Interpreter', 'latex', 'FontSize', 7);
    
    % Save summary plot
    try
        saveas(gcf, fullfile(figures_dir, 'summary_sensitivities_comparison.png'));
        saveas(gcf, fullfile(figures_dir, 'summary_sensitivities_comparison.fig'));
        fprintf('  ? Summary plot saved\n');
    catch ME
        fprintf('  ? Error saving summary plot: %s\n', ME.message);
    end
    
    fprintf('  ? Generated summary comparison plot (figure 999)\n');
end

function save_results(results, successful_fields, sensitivities_dir, work_dir)

    % Save all results and sensitivities
    
    try
        
        % Save individual field sensitivities (final derivatives only)
        for i = 1:length(successful_fields)
            field_name = successful_fields{i};
            field_result = results.(field_name);
            
            % Create sensitivity data structure
            sensitivity_data = struct();
            sensitivity_data.sensitivity = field_result.final_derivative;  % The main result
            sensitivity_data.field_name = field_name;
            sensitivity_data.latex_label = field_result.latex_label;
            sensitivity_data.derivative_label = field_result.derivative_label;
            
            % Save individual sensitivity file
            sensitivity_file = fullfile(sensitivities_dir, sprintf('sensitivity_%s.mat', field_name));
            save(sensitivity_file, '-struct', 'sensitivity_data');
            
        end
       
        
        fprintf('  ? Results saved to sensitivities directory: %s\n', sensitivities_dir);
        fprintf('    - Individual sensitivities: sensitivity_<fieldname>.mat\n');
        
    catch ME
        fprintf('  ? Error saving results: %s\n', ME.message);
    end
end
