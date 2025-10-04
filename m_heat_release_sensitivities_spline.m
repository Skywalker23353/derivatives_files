% Spline-based sensitivity computation script
% This script computes d(heat_release)/d(field) using spline fitting and analytical differentiation
% Results are evaluated directly on the RZ plane using ppval

clear; clc; close all;
addpath('~/MATLAB/');
addpath('~/satyam_files/CH4_jet_PF/2025_Runs/derivatives_files/Interpolate_derivs_on_R_Z_plane/functions');

%% Configuration
work_dir = '/work/home/satyam/satyam_files/CH4_jet_PF/2025_Runs/derivatives_files/Interpolate_derivs_on_R_Z_plane';
deriv_dir = '/work/home/satyam/satyam_files/CH4_jet_PF/2025_Runs/derivatives_files/sensitivities_10D';
data_dir = '/work/home/satyam/satyam_files/CH4_jet_PF/2025_Runs/c_cond_stats/C_cond_fields_800_10D_n';

write_to_h5_file_flag = false;
h5filename = 'Reactants_1';
h5_outdir = '/work/home/satyam/satyam_files/CH4_jet_PF/2025_Runs/LES_base_case_v6/filtering_run3/sensitivities/10D';

D = 2e-3;
window = 3; % Window size for smoothing
rmx = 5;
zmx = 10;
Yu = 0.222606; % Yb = 0.0423208;
Yb = 0.0411;
Min_c_limit = 1e-3; Max_c_limit = 0.98;
load("comb_interpolted_hrr_field_10D.mat","model_scaling_factor"); % hrr scaling factor
l_ref = 2e-3;
U_ref = 65;
V_ref = l_ref^3;
Cp_ref = 1100;
variable_ref_val_list = {'density',0.4237; 'Temperature',800;};
Q_bar = 163;

% Define fields to process (same as original script)
field_configs = {
    % {field_name, smooth_suffix, latex_label, plot_title, derivative_label}
    'Temperature', '_smooth_bc_smooth', 'T', '$T$', '$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial T}$';
    'density', '_smooth_bc_smooth', '\rho', '$\rho$', '$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial \rho}$';
    'CH4', '_smooth_bc_smooth', 'Y_{CH_4}', '$Y_{CH4}$', '$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial Y_{CH4}}$';
    'O2', '_smooth_bc_smooth', 'Y_{O_2}', '$Y_{O2}$', '$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial Y_{O2}}$';
    'CO2', '_smooth_bc_smooth', 'Y_{CO_2}', '$Y_{CO2}$', '$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial Y_{CO2}}$';
    'H2O', '_smooth_bc_smooth', 'Y_{H_2O}', '$Y_{H2O}$', '$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial Y_{H2O}}$';
};

%% Step 1: Load grids and setup
fprintf('=== Loading Grid Data and Setup ===\n');

% Load combustor and nozzle coordinate data
fprintf('Loading combustor domain data...\n');
case_dir = 'interp_16_unilat_structured_pl_combustor_domain';
comb_data = load(sprintf('%s/%s/interp_fields/YO2_field_structured_grid_16_planes_1011260202_avg.mat', work_dir, case_dir));

fprintf('Loading nozzle domain data...\n');
case_dir = 'interp_16_unilat_structured_pl_nozzle_domain_with_zero';
noz_data = load(sprintf('%s/%s/interp_fields/YO2_field_structured_grid_16_planes_1011260202_avg.mat', work_dir, case_dir));

% Load structured grid coordinates
fprintf('Loading structured grid coordinates...\n');
load('structured_grid_from_LES_grid_with_zero_16_unilat_planes.mat');

% Create C field using YO2 field
fprintf('Computing C fields...\n');
comb_data.C_field_MAT = (Yu - comb_data.O2_mean)/(Yu - Yb);
noz_data.C_field_MAT = (Yu - noz_data.O2_mean)/(Yu - Yb);

% Load reference C_MAT and Z_MAT for sensitivities
fprintf('Loading reference C and Z matrices...\n');
load(fullfile(data_dir, 'CZ_data.mat'), "C_MAT", "Z_MAT");

% Apply Z restriction logic
Z_idx_mx = find((Z_MAT)/D >= zmx, 1);
r_idx_mx = find(R1(1,:)/D >= rmx,1);
R1 = R1(1:Z_idx_mx, 1:r_idx_mx);
Z1 = Z1(1:Z_idx_mx, 1:r_idx_mx);
comb_data.C_field_MAT = comb_data.C_field_MAT(1:Z_idx_mx, 1:r_idx_mx);

% Sanity check for C field
comb_data.C_field_MAT(find(comb_data.C_field_MAT < Min_c_limit)) = 0;
comb_data.C_field_MAT(find(comb_data.C_field_MAT > Max_c_limit)) = 1;

fprintf('Setup complete. Grid size: %dx%d\n', size(Z1, 1), size(Z1, 2));

%% Step 2: Load conditional statistics data
fprintf('\n=== Loading Conditional Statistics Data ===\n');

% Load heat release data
fprintf('Loading heat release data...\n');
heat_release_field = 'Heatrelease_smooth';
try
    heatRelease = load(sprintf('%s/%s.mat', data_dir, heat_release_field));
    fprintf('✓ Heat release data loaded successfully\n');
catch ME
    fprintf('✗ Error loading heat release data: %s\n', ME.message);
    return;
end

[C_MAT_stats, Z_MAT_stats] = get_CZ_coord_data(data_dir);

%% Step 3: Initialize results structures
fprintf('\n=== Initializing Results Structures ===\n');

comb_sensitivities = struct();
noz_sensitivities = struct();

% Copy coordinate and C field data to results
comb_sensitivities.R1 = R1;
comb_sensitivities.Z1 = Z1;
comb_sensitivities.C_field_MAT = comb_data.C_field_MAT;
comb_sensitivities.Yu = Yu;
comb_sensitivities.Yb = Yb;

noz_sensitivities.R2 = R2;
noz_sensitivities.Z2 = Z2;
noz_sensitivities.C_field_MAT = noz_data.C_field_MAT;
noz_sensitivities.Yu = Yu;
noz_sensitivities.Yb = Yb;

%% Step 4: Process heat release data with splines
fprintf('\n=== Processing Heat Release with Splines ===\n');

% Initialize heat release derivative matrix on RZ grid
comb_hrr_derivative = zeros(size(R1));
noz_hrr_derivative = zeros(size(R2));

% Process each z-level
for i = 1:size(Z1, 1)
    if mod(i, 10) == 0 || i == 1 || i == size(Z1, 1)
        fprintf('Processing heat release Z level %d/%d...\n', i, size(Z1, 1));
    end

    % Get c values for this z-level from RZ grid
    c_values_rz = comb_data.C_field_MAT(i, :);

    % Fit spline to heat release data at this z-level
    try
        hrr_spline = spline(C_MAT_stats(1, :), heatRelease.DF(i, :));
        hrr_deriv_spline = fnder(hrr_spline, 1); % First derivative

        % Evaluate derivative at RZ c-values using ppval
        comb_hrr_derivative(i, :) = ppval(hrr_deriv_spline, c_values_rz);
    catch ME
        fprintf('Warning: Spline fitting failed for heat release at z-level %d: %s\n', i, ME.message);
    end
end

% Apply scaling to heat release derivative
comb_hrr_derivative = model_scaling_factor * comb_hrr_derivative * V_ref / Q_bar;

% Store heat release derivative
comb_sensitivities.heat_release_derivative = comb_hrr_derivative;

%% Step 5: Process each field with splines
fprintf('\n=== Processing Fields with Splines ===\n');

for field_idx = 1:size(field_configs, 1)
    field_name = field_configs{field_idx, 1};
    smooth_suffix = field_configs{field_idx, 2};
    clean_field_name = field_name;

    fprintf('\n--- Processing %d/%d: %s ---\n', field_idx, size(field_configs, 1), clean_field_name);

    % Construct field filename
    field_filename = sprintf('%s%s', field_name, smooth_suffix);
    field_path = sprintf('%s/%s.mat', data_dir, field_filename);

    % Check if field exists
    if ~exist(field_path, 'file')
        fprintf('✗ Field file not found: %s\n', field_path);
        % Try without smooth suffix
        field_filename_alt = field_name;
        field_path_alt = sprintf('%s/%s.mat', data_dir, field_filename_alt);
        if exist(field_path_alt, 'file')
            fprintf('  Found alternative: %s\n', field_filename_alt);
            field_path = field_path_alt;
        else
            fprintf('  Skipping field: %s\n', field_name);
            continue;
        end
    end

    try
        fprintf('  Loading field data...\n');
        field_data = load(field_path);

        % Initialize field derivative matrices on RZ grid
        comb_field_derivative = zeros(size(R1));
        noz_field_derivative = zeros(size(R2));

        %% Step 5a: Process combustor domain with splines
        fprintf('  Processing combustor domain with splines...\n');

        for i = 1:size(Z1, 1)
            if mod(i, 20) == 0 || i == 1 || i == size(Z1, 1)
                fprintf('    Processing Z level %d/%d...\n', i, size(Z1, 1));
            end

            % Get c values for this z-level from RZ grid
            c_values_rz = comb_data.C_field_MAT(i, :);

            % Fit spline to field data at this z-level
            try
                field_spline = spline(C_MAT_stats(1, :), field_data.DF(i, :));
                field_deriv_spline = fnder(field_spline, 1); % First derivative

                % Evaluate derivative at RZ c-values using ppval
                comb_field_derivative(i, :) = ppval(field_deriv_spline, c_values_rz);
            catch ME
                fprintf('    Warning: Spline fitting failed at z-level %d: %s\n', i, ME.message);
            end
        end

        % Apply scaling based on field type
        if strcmp(clean_field_name, "Temperature") || strcmp(clean_field_name, "density")
            variable_ref_val = variable_ref_val_list{strcmp(variable_ref_val_list(:,1), clean_field_name), 2};
            comb_field_derivative = model_scaling_factor * comb_field_derivative * variable_ref_val * V_ref / Q_bar;
        else
            comb_field_derivative = model_scaling_factor * comb_field_derivative * V_ref / Q_bar;
        end

        % Handle zero values in denominator
        zero_idx = find(abs(comb_field_derivative) < 1e-16);
        if ~isempty(zero_idx)
            fprintf('  Handling %d zero values in field derivative...\n', length(zero_idx));
            comb_field_derivative(zero_idx) = 1e-8;
        end

        % Compute sensitivity: d(heat_release)/d(field)
        comb_sensitivity = comb_hrr_derivative ./ comb_field_derivative;

        % % Apply boundary conditions
        % z_ref = 8*D; % Default z_ref
        % if strcmp(clean_field_name, 'CH4')
        %     z_ref = 8.5*D; % Special case for CH4
        % end
        % comb_sensitivity = set_sensitivities_bc(comb_sensitivity, Z1, z_ref);

        % Store result
        comb_sensitivities.(clean_field_name) = comb_sensitivity;

        %% Step 5b: Process nozzle domain
        fprintf('  Processing nozzle domain...\n');

        % Set boundary values from combustor
        nozzle_field = zeros(size(R2));
        nozzle_field(end,:) = comb_sensitivity(1, 1:size(R2, 2));

        % Apply smoothing to nozzle
        noz_sensitivities.(clean_field_name) = myutils.f_return_smooth_field(nozzle_field, window, 'col');

        fprintf('  ✓ Successfully processed %s\n', clean_field_name);

    catch ME
        fprintf('  ✗ Error processing %s: %s\n', field_name, ME.message);
        continue;
    end
end

%% Step 6: Save results
fprintf('\n=== Saving Results ===\n');

% Create output directory if it doesn't exist
output_dir = './spline_sensitivities_output';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('Created output directory: %s\n', output_dir);
end

% Save combustor sensitivities
comb_output_file = fullfile(output_dir, 'combustor_spline_sensitivities.mat');
save(comb_output_file, '-struct', 'comb_sensitivities');
fprintf('✓ Saved combustor sensitivities to: %s\n', comb_output_file);

% Save nozzle sensitivities
noz_output_file = fullfile(output_dir, 'nozzle_spline_sensitivities.mat');
save(noz_output_file, '-struct', 'noz_sensitivities');
fprintf('✓ Saved nozzle sensitivities to: %s\n', noz_output_file);

%% Step 7: Write to H5 file if requested
if write_to_h5_file_flag
    fprintf('\n=== Writing to H5 File ===\n');

    % Get successfully processed field names
    processed_fields = {};
    for i = 1:size(field_configs, 1)
        field_name = field_configs{i, 1};
        if isfield(comb_sensitivities, field_name)
            processed_fields{end+1} = field_name;
        end
    end

    N_fields = length(processed_fields);
    inp_fieldNames = processed_fields;

    try
        write_field_to_h5_file(inp_fieldNames, N_fields, comb_sensitivities, noz_sensitivities, h5_outdir, h5filename);
        fprintf('✓ Successfully wrote to H5 file: %s/%s.h5\n', h5_outdir, h5filename);
    catch ME
        fprintf('✗ Error writing to H5 file: %s\n', ME.message);
    end
end

%% Summary
fprintf('\n=== Processing Summary ===\n');
fprintf('Spline-based sensitivity computation completed.\n');
fprintf('Method: Analytical differentiation of fitted splines\n');
fprintf('Evaluation: Direct evaluation on RZ grid using ppval()\n');
fprintf('Combustor grid size: %dx%d\n', size(Z1, 1), size(Z1, 2));
fprintf('Nozzle grid size: %dx%d\n', size(Z2, 1), size(Z2, 2));
fprintf('Files saved to: %s\n', output_dir);

if write_to_h5_file_flag
    fprintf('H5 file saved to: %s/%s.h5\n', h5_outdir, h5filename);
end

fprintf('\n=== Processing Complete ===\n');
