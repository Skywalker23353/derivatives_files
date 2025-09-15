% Enhanced script to compute derivatives for multiple species reaction rates
% This script computes d(species_reaction_rate)/d(field) for multiple species and fields
% where numerator is species reaction rate and denominator is various fields

clear; clc;
% close all;
addpath('~/MATLAB');
addpath('~/satyam_files/CH4_jet_PF/2025_Runs/derivatives_files/functions');

%% User Input - Work Directory
work_dir = '/work/home/satyam/satyam_files/CH4_jet_PF/2025_Runs/derivatives_files';

%% Configuration
data_dir = '/work/home/satyam/satyam_files/CH4_jet_PF/2025_Runs/c_cond_stats/C_cond_fields_800';

D = 2e-3;  % Diameter scale
save_results_flag = true;
generate_plots_flag = false;
threshold_and_smooth_results_flag = true;
smth_window = 3;

%% Load coordinate data (C_MAT and Z_MAT)
fprintf('Loading coordinate data from CZ_data.mat...\n');
try
    coord_data = load(fullfile(data_dir, 'CZ_data.mat'));
    C_MAT = coord_data.C_MAT;
    Z_MAT = coord_data.Z_MAT;
    fprintf(' Successfully loaded coordinate data\n');
catch ME
    fprintf(' Error loading coordinate data: %s\n', ME.message);
    fprintf('Please ensure CZ_data.mat exists in: %s\n', data_dir);
    return;
end

% Create output directories
sensitivities_dir = fullfile(work_dir, 'species_sensitivities');
figures_dir = fullfile(work_dir, 'species_figures');

if ~exist(sensitivities_dir, 'dir')
    mkdir(sensitivities_dir);
    fprintf('Created species sensitivities directory: %s\n', sensitivities_dir);
end

if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
    fprintf('Created species figures directory: %s\n', figures_dir);
end

% Define species to process (numerators)
species_configs = {
    % {species_name, smooth_suffix, file_prefix, latex_label, plot_title}
    'CH4', '_smooth', 'SYm_CH4', '\dot{\omega}_{CH_4}', 'CH_4';
    'O2', '_smooth', 'SYm_O2', '\dot{\omega}_{O_2}', 'O_2';
    'CO2', '_smooth', 'SYm_CO2', '\dot{\omega}_{CO_2}', 'CO_2';
    'H2O', '_smooth', 'SYm_H2O', '\dot{\omega}_{H_2O}', 'H_2O';
};

% Define fields to process (denominators)
field_configs = {
    % {field_name, smooth_suffix, latex_label, plot_title, short_name}
    'Temperature', '_smooth', 'T', '$T$', 'T';
    'density', '_smooth', '\rho', '$\rho$', 'rho';
    'CH4', '_smooth', 'Y_{CH_4}', '$Y_{CH4}$', 'CH4';
    'O2', '_smooth', 'Y_{O_2}', '$Y_{O2}$', 'O2';
    'CO2', '_smooth', 'Y_{CO_2}', '$Y_{CO2}$', 'CO2';
    'H2O', '_smooth', 'Y_{H_2O}', '$Y_{H2O}$', 'H2O';
};

%% Load all field data first (denominators)
fprintf('Loading field data (denominators)...\n');
field_data = struct();
successful_fields = {};

for i = 1:size(field_configs, 1)
    field_name = field_configs{i, 1};
    smooth_suffix = field_configs{i, 2};
    latex_label = field_configs{i, 3};
    plot_title = field_configs{i, 4};
    short_name = field_configs{i, 5};
    
    fprintf('  Loading %s...\n', field_name);
    
    % Construct field filename
    field_filename = sprintf('%s%s', field_name, smooth_suffix);
    field_path = sprintf('%s/%s.mat', data_dir, field_filename);
    
    % Check if field exists
    if ~exist(field_path, 'file')
        fprintf('    Field file not found: %s\n', field_path);
        % Try without smooth suffix
        field_filename_alt = field_name;
        field_path_alt = sprintf('%s/%s.mat', data_dir, field_filename_alt);
        if exist(field_path_alt, 'file')
            fprintf('    Found alternative: %s\n', field_filename_alt);
            field_path = field_path_alt;
        else
            fprintf('    Skipping field: %s\n', field_name);
            continue;
        end
    end
    
    try
        % Load field data
        loaded_data = load(field_path);
        
        % Store field information
        field_data.(field_name) = struct();
        field_data.(field_name).data = loaded_data.DF;
        field_data.(field_name).C_MAT = C_MAT;  % Use loaded coordinate data
        field_data.(field_name).Z_MAT = Z_MAT;  % Use loaded coordinate data
        field_data.(field_name).latex_label = latex_label;
        field_data.(field_name).plot_title = plot_title;
        field_data.(field_name).short_name = short_name;
        
        % Compute field derivative
        fprintf('    Computing derivative d%s/dc...\n', field_name);
        field_data.(field_name).derivative = compute_dfdr(loaded_data.DF, C_MAT);
        
        successful_fields{end+1} = field_name;
        fprintf('    ? Successfully loaded %s\n', field_name);
        
    catch ME
        fprintf('    ? Error loading %s: %s\n', field_name, ME.message);
        continue;
    end
end

%% Process each species (numerators)
fprintf('\nProcessing %d species...\n', size(species_configs, 1));
results = struct();
successful_species = {};

for i = 1:size(species_configs, 1)
    species_name = species_configs{i, 1};
    smooth_suffix = field_configs{i, 2};
    file_prefix = species_configs{i, 3};
    species_latex = species_configs{i, 4};
    species_title = species_configs{i, 5};
    
    fprintf('\n--- Processing Species %d/%d: %s ---\n', i, size(species_configs, 1), species_name);
    
    % Load species reaction rate data
%     species_filename = file_prefix;
    species_filename = sprintf('%s%s', file_prefix, smooth_suffix);
    species_path = sprintf('%s/%s.mat', data_dir, species_filename);
    
    if ~exist(species_path, 'file')
        fprintf('  Species file not found: %s\n', species_path);
        continue;
    end
    
    try
        % Load species data
        fprintf('  Loading species reaction rate data...\n');
        species_data = load(species_path);
        
        % Compute species derivative with respect to c
        fprintf('  Computing derivative d%s/dc...\n', species_name);
        dspecies_dc = compute_dfdr(species_data.DF, C_MAT);
        
        % Initialize species results
        results.(species_name) = struct();
        results.(species_name).species_data = species_data.DF;
        results.(species_name).species_derivative = dspecies_dc;
        results.(species_name).C_MAT = C_MAT;  % Use loaded coordinate data
        results.(species_name).Z_MAT = Z_MAT;  % Use loaded coordinate data
        results.(species_name).latex_label = species_latex;
        results.(species_name).plot_title = species_title;
        results.(species_name).field_derivatives = struct();
        
        % Process each field for this species
        successful_combinations = {};
        
        for j = 1:length(successful_fields)
            field_name = successful_fields{j};
            field_info = field_data.(field_name);
            
            fprintf('  Computing d%s/d%s...\n', species_name, field_name);
            
            try
                % Get field derivative
                df_dc = field_info.derivative;
                
                % Handle zero values in denominator
                zero_idx = find(abs(df_dc) < 1e-16);
                if ~isempty(zero_idx)
                    fprintf('    Handling %d zero values in denominator...\n', length(zero_idx));
                    df_dc(zero_idx) = 1e-8;
                end
                
                % Compute final derivative: d(species_reaction_rate)/d(field)
                final_derivative = dspecies_dc ./ df_dc;
                if threshold_and_smooth_results_flag
                    final_derivative = custom_thresholding_and_smoothing(final_derivative,sensitivities_dir,species_name,field_data.(field_name).short_name,smth_window);
                end
                % Store results
                results.(species_name).field_derivatives.(field_name) = struct();
                results.(species_name).field_derivatives.(field_name).field_data = field_info.data;
                results.(species_name).field_derivatives.(field_name).field_derivative = df_dc;
                results.(species_name).field_derivatives.(field_name).final_derivative = final_derivative;
                results.(species_name).field_derivatives.(field_name).latex_label = field_info.latex_label;
                results.(species_name).field_derivatives.(field_name).plot_title = field_info.plot_title;
                results.(species_name).field_derivatives.(field_name).short_name = field_info.short_name;
                
                successful_combinations{end+1} = sprintf('%s_%s', species_name, field_name);
                fprintf('     Successfully computed d%s/d%s\n', species_name, field_name);
                
            catch ME
                fprintf('     Error computing d%s/d%s: %s\n', species_name, field_name, ME.message);
                continue;
            end
        end
        
        if ~isempty(successful_combinations)
            successful_species{end+1} = species_name;
            fprintf('   Successfully processed %s with %d field combinations\n', species_name, length(successful_combinations));
        end
        
    catch ME
        fprintf('   Error processing %s: %s\n', species_name, ME.message);
        continue;
    end
end

%% Generate plots and save results
if ~isempty(successful_species)
    fprintf('\n=== Generating Plots and Saving Results ===\n');
    fprintf('Successfully processed %d species: %s\n', length(successful_species), strjoin(successful_species, ', '));
    
    % Generate plots for each species
    for i = 1:length(successful_species)
        species_name = successful_species{i};
        species_result = results.(species_name);
        
        if generate_plots_flag
            fprintf('\nGenerating plots for %s...\n', species_name);
        
            % Type 1 plots: Individual derivative analysis (3 subplots each)
            generate_type1_plots(species_result, species_name, successful_fields, D, figures_dir);
            
            % Type 2 plot: All field derivatives for this species (6 subplots)
            generate_type2_plot(species_result, species_name, successful_fields, D, figures_dir);
        else
            fprintf("Figures flag set to false");
        end

        % Save results
        if save_results_flag
            save_species_results(species_result, species_name, successful_fields, sensitivities_dir);
        end
    end
    
else
    fprintf('\n No species were successfully processed.\n');
end

fprintf('\n=== Processing Complete ===\n');

%% Function Definitions

function generate_type1_plots(species_result, species_name, successful_fields, D, figures_dir)
    % Generate Type 1 plots: Individual derivative analysis (3 subplots each)
    % This creates 24 figures total (6 fields × 4 species)
    
    field_names = fieldnames(species_result.field_derivatives);
    
    for i = 1:length(field_names)
        field_name = field_names{i};
        field_result = species_result.field_derivatives.(field_name);
        
        % Create figure
        figure_num = 1000 + i + length(successful_fields) * find(strcmp({'CH4', 'O2', 'CO2', 'H2O'}, species_name));
        figure(figure_num);
        set(gcf, 'WindowState', 'maximized');
        
        % Subplot 1: Species derivative w.r.t. c
        subplot(1, 3, 1);
        try
            myutils.plot_contourf(gcf, species_result.C_MAT, species_result.Z_MAT/D, species_result.species_derivative, ...
                '$c$', '$z/D$', sprintf('$\\frac{\\partial \\langle %s|c\\rangle}{\\partial c}$', species_result.latex_label));
        catch
            contourf(species_result.C_MAT, species_result.Z_MAT/D, species_result.species_derivative, 20);
            colorbar;
            xlabel('$c$', 'Interpreter', 'latex');
            ylabel('$z/D$', 'Interpreter', 'latex');
            title(sprintf('$\\frac{\\partial \\langle %s|c\\rangle}{\\partial c}$', species_result.latex_label), 'Interpreter', 'latex');
        end
        
        % Subplot 2: Field derivative w.r.t. c
        subplot(1, 3, 2);
        try
            myutils.plot_contourf(gcf, species_result.C_MAT, species_result.Z_MAT/D, field_result.field_derivative, ...
                '$c$', '$z/D$', sprintf('$\\frac{\\partial \\langle %s|c\\rangle}{\\partial c}$', field_result.latex_label));
        catch
            contourf(species_result.C_MAT, species_result.Z_MAT/D, field_result.field_derivative, 20);
            colorbar;
            xlabel('$c$', 'Interpreter', 'latex');
            ylabel('$z/D$', 'Interpreter', 'latex');
            title(sprintf('$\\frac{\\partial \\langle %s|c\\rangle}{\\partial c}$', field_result.latex_label), 'Interpreter', 'latex');
        end
        
        % Subplot 3: Final derivative (sensitivity)
        subplot(1, 3, 3);
        try
            myutils.plot_contourf(gcf, species_result.C_MAT, species_result.Z_MAT/D, field_result.final_derivative, ...
                '$c$', '$z/D$', sprintf('$\\frac{\\partial \\langle %s|c\\rangle}{\\partial %s}$', species_result.latex_label, field_result.latex_label));
        catch
            contourf(species_result.C_MAT, species_result.Z_MAT/D, field_result.final_derivative, 20);
            colorbar;
            xlabel('$c$', 'Interpreter', 'latex');
            ylabel('$z/D$', 'Interpreter', 'latex');
            title(sprintf('$\\frac{\\partial \\langle %s|c\\rangle}{\\partial %s}$', species_result.latex_label, field_result.latex_label), 'Interpreter', 'latex');
        end
        
        % Add overall title
        sgtitle(sprintf('%s Derivative Analysis w.r.t. %s', species_name, field_name), 'FontSize', 16, 'FontWeight', 'bold');
        
        % Save figure
        try
            saveas(gcf, fullfile(figures_dir, sprintf('Type1_%s_vs_%s_analysis.png', species_name, field_name)));
            saveas(gcf, fullfile(figures_dir, sprintf('Type1_%s_vs_%s_analysis.fig', species_name, field_name)));
            fprintf('   Type 1 figure saved: %s vs %s\n', species_name, field_name);
        catch ME
            fprintf('   Error saving Type 1 figure for %s vs %s: %s\n', species_name, field_name, ME.message);
        end
    end
end

function generate_type2_plot(species_result, species_name, successful_fields, D, figures_dir)
    % Generate Type 2 plot: All field derivatives for this species (6 subplots)
    % This creates 4 figures total (one per species)
    
    field_names = fieldnames(species_result.field_derivatives);
    
    % Create figure
    figure_num = 2000 + find(strcmp({'CH4', 'O2', 'CO2', 'H2O'}, species_name));
    figure(figure_num);
    set(gcf, 'WindowState', 'maximized');
    
    % Calculate subplot arrangement
    n_fields = length(field_names);
    subplot_rows = 2;
    subplot_cols = 3;
    
    for i = 1:length(field_names)
        field_name = field_names{i};
        field_result = species_result.field_derivatives.(field_name);
        
        subplot(subplot_rows, subplot_cols, i);
        try
            myutils.plot_contourf(gcf, species_result.C_MAT, species_result.Z_MAT/D, field_result.final_derivative, ...
                '$c$', '$z/D$', sprintf('$\\frac{\\partial %s}{\\partial %s}$', species_result.latex_label, field_result.latex_label));
        catch
            contourf(species_result.C_MAT, species_result.Z_MAT/D, field_result.final_derivative, 20);
            colorbar;
            xlabel('$c$', 'Interpreter', 'latex');
            ylabel('$z/D$', 'Interpreter', 'latex');
            title(sprintf('$\\frac{\\partial %s}{\\partial %s}$', species_result.latex_label, field_result.latex_label), 'Interpreter', 'latex');
        end
    end
    
    % Add overall title
    sgtitle(sprintf('%s Reaction Rate Derivatives w.r.t. All Fields', species_name), 'FontSize', 16, 'FontWeight', 'bold');
    
    % Save figure
    try
        saveas(gcf, fullfile(figures_dir, sprintf('Type2_%s_all_derivatives.png', species_name)));
        saveas(gcf, fullfile(figures_dir, sprintf('Type2_%s_all_derivatives.fig', species_name)));
        fprintf('   Type 2 figure saved: %s all derivatives\n', species_name);
    catch ME
        fprintf('   Error saving Type 2 figure for %s: %s\n', species_name, ME.message);
    end
end

function save_species_results(species_result, species_name, successful_fields, sensitivities_dir)
    % Save results for each species-field combination
    
    field_names = fieldnames(species_result.field_derivatives);
    
    for i = 1:length(field_names)
        field_name = field_names{i};
        field_result = species_result.field_derivatives.(field_name);
        short_name = field_result.short_name;
        
        % Create result data structure
        sensitivity_data = struct();
        sensitivity_data.sensitivity = field_result.final_derivative;
        sensitivity_data.species_name = species_name;
        sensitivity_data.field_name = field_name;
        sensitivity_data.species_latex = species_result.latex_label;
        sensitivity_data.field_latex = field_result.latex_label;
        
        % Save with naming convention: dCH4_dT.mat, dCH4_drho.mat, etc.
        output_filename = sprintf('d%s_d%s.mat', species_name, short_name);
        output_path = fullfile(sensitivities_dir, output_filename);
        
        try
            save(output_path, '-struct', 'sensitivity_data');
            fprintf('   Saved: %s\n', output_filename);
        catch ME
            fprintf('   Error saving %s: %s\n', output_filename, ME.message);
        end
    end
end
