% Enhanced script to compute derivatives for multiple species reaction rates
% This script computes d(species_reaction_rate)/d(field) for multiple species and fields
% where numerator is species reaction rate and denominator is various fields

clear; clc;
close all;
addpath('~/MATLAB');
addpath(fullfile(pwd, 'functions'));

%% User Input - Work Directory
work_dir = '/work/home/satyam/satyam_files/CH4_jet_PF/2025_Runs/derivatives_files';

%% Configuration
data_dir = '/work/home/satyam/satyam_files/CH4_jet_PF/2025_Runs/c_cond_stats/C_cond_fields_800';

D = 2e-3;  % Diameter scale
save_results_flag = true;
generate_plots_flag = false;
threshold_and_smooth_results_flag = true;
smth_window = 3;
plot_surface_figures_flag = true;  % Toggle for plotting surface figures
save_surface_figures_flag = true;  % Toggle for saving surface figures

%% Load coordinate data (C_MAT and Z_MAT)
[C_MAT, Z_MAT] = get_CZ_coord_data(data_dir);

% Check if coordinate data was loaded successfully
if isempty(C_MAT) || isempty(Z_MAT)
    fprintf('Failed to load coordinate data. Exiting...\n');
    return;
end

% Create output directories
dirs = create_output_dirs(work_dir);
sensitivities_dir = dirs.sensitivities_dir;
figures_dir = dirs.figures_dir;
surface_figures_dir = dirs.surface_figures_dir;

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
[fields_struct, successful_fields] = compute_field_derivatives(field_configs, data_dir, C_MAT, Z_MAT);

%% Compute species numerators (derivatives with respect to c)
[species_struct, successful_species] = compute_species_numerators(species_configs, field_configs, data_dir, C_MAT, Z_MAT);

%% Compute final derivatives/sensitivities
species_struct = compute_sensitivities(species_struct, fields_struct, successful_species, successful_fields, threshold_and_smooth_results_flag, sensitivities_dir, smth_window);

%% Generate plots and save results
if ~isempty(successful_species)
    fprintf('\n=== Generating Plots and Saving Results ===\n');
    fprintf('Successfully processed %d species: %s\n', length(successful_species), strjoin(successful_species, ', '));
    if plot_surface_figures_flag || save_surface_figures_flag
        plot_and_save_individual_field_surface_figures(fields_struct,D, surface_figures_dir, plot_surface_figures_flag, save_surface_figures_flag)
    end    
    % Generate plots for each species
    for i = 1:length(successful_species)
        species_name = successful_species{i};
        species_result = species_struct.(species_name);

        if generate_plots_flag
            fprintf('\nGenerating plots for %s...\n', species_name);
            generate_type1_plots(species_struct, fields_struct, species_name, successful_fields, D, figures_dir);
            generate_type2_plot(species_struct, fields_struct, species_name, successful_fields, D, figures_dir);
        else
            fprintf("Figures flag set to false");
        end

        % Plot and/or save individual surface figures based on flags
        if plot_surface_figures_flag || save_surface_figures_flag
            save_individual_surface_figures(species_struct, fields_struct, species_name, successful_fields, D, surface_figures_dir, plot_surface_figures_flag, save_surface_figures_flag,i);
        end
    end
     % Save results
        if save_results_flag
            save_species_results(species_struct, fields_struct, sensitivities_dir);
        end
else
    fprintf('\n No species were successfully processed. \n');
end
fprintf('\n=== Processing Complete ===\n');