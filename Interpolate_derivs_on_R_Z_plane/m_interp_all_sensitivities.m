clear; clc;close all;
addpath('~/MATLAB/');
addpath('~/satyam_files/CH4_jet_PF/2025_Runs/derivatives_files/Interpolate_derivs_on_R_Z_plane/functions');

%% Configuration
work_dir = '/work/home/satyam/satyam_files/CH4_jet_PF/2025_Runs/derivatives_files/Interpolate_derivs_on_R_Z_plane';
deriv_dir = '/work/home/satyam/satyam_files/CH4_jet_PF/2025_Runs/derivatives_files/sensitivities';

% parameters
write_to_h5_file_flag = false;
h5filename = 'Reactants';
h5_outdir = '/work/home/satyam/satyam_files/CH4_jet_PF/2025_Runs/LES_base_case_v6/filtering_run3/sensitivities';
D = 2e-3;
window = 3; % Window size for nozzle data smoothening (adjust as needed)
rmx = 5;
zmx=8.5;
Yu = 0.222606; %Yb = 0.0423208;
Yb = 0.0399;
load("comb_interpolted_hrr_field.mat","model_scaling_factor"); %hrr scaling factor
l_ref = 2e-3;
U_ref = 65;
V_ref = l_ref^3;
Cp_ref = 1100;
variable_ref_val_list = {'density',0.4237; 'Temperature',800;}; 
Q_bar = 163;

% Define all sensitivity fields to process
sensitivity_fields = {
    'sensitivity_Temperature';
    'sensitivity_density';
    'sensitivity_CH4';
    'sensitivity_O2';
    'sensitivity_CO2';
    'sensitivity_H2O';
    'sensitivity_N2';
};

%% Step 1-4: Load grids and setup (done once)
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
comb_data.C_field_MAT = (Yu - comb_data.O2_fmean)/(Yu - Yb);
noz_data.C_field_MAT = (Yu - noz_data.O2_fmean)/(Yu - Yb);

% Load reference C_MAT and Z_MAT for sensitivities
fprintf('Loading reference C and Z matrices...\n');
load('/work/home/satyam/satyam_files/CH4_jet_PF/2025_Runs/c_cond_stats/C_cond_fields_800/Heatrelease_smooth.mat', "C_MAT", "Z_MAT");

% Apply Z restriction logic
Z_idx_mx = find((Z_MAT)/D >= zmx, 1);
r_idx_mx = find(R1(1,:)/D >= rmx,1);
R1 = R1(1:Z_idx_mx, 1:r_idx_mx); 
Z1 = Z1(1:Z_idx_mx, 1:r_idx_mx);
comb_data.C_field_MAT = comb_data.C_field_MAT(1:Z_idx_mx, 1:r_idx_mx);

% Sanity check for C field
% comb_data.C_field_MAT(find(comb_data.C_field_MAT < 1e-5)) = 0;
comb_data.C_field_MAT(find(comb_data.C_field_MAT > 1)) = 1;

fprintf('Setup complete. Grid size: %dx%d\n', size(Z1, 1), size(Z1, 2));

%% Step 5-7: Process all sensitivity fields
fprintf('\n=== Processing Sensitivity Fields ===\n');

% Initialize results structures
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

% Process each sensitivity field
for field_idx = 1:length(sensitivity_fields)
    field_name = sensitivity_fields{field_idx};
    fprintf('\n--- Processing %d/%d: %s ---\n', field_idx, length(sensitivity_fields), field_name);
    
    % Load sensitivity data
    sensitivity_file = sprintf('%s/%s.mat', deriv_dir, field_name);
    
    if ~exist(sensitivity_file, 'file')
        fprintf('? Warning: Sensitivity file not found: %s\n', sensitivity_file);
        continue;
    end
    
    try
        fprintf('  Loading sensitivity data...\n');
        deriv_data = load(sensitivity_file);
        
        % Extract clean field name for storage (remove 'sensitivity_' prefix)
        clean_field_name = strrep(field_name, 'sensitivity_', '');
        
        %% Step 5: Interpolate sensitivity onto combustor grid
        fprintf('  Interpolating onto combustor grid...\n');
        Interp_deriv_field = zeros(size(comb_data.C_field_MAT));
        
        for i = 1:size(Z1, 1)
            if mod(i, 10) == 0 || i == 1 || i == size(Z1, 1)
                fprintf('    Processing Z level %d/%d...\n', i, size(Z1, 1));
            end
            Interp_deriv_field(i, :) = linear_interp_field(comb_data.C_field_MAT(i, :)', ...
                                                          C_MAT(1, :)', ...
                                                          deriv_data.sensitivity(i, :)');
        end
        
        % Store interpolated result in combustor structure
        if strcmp(clean_field_name,"Temperature") || strcmp(clean_field_name,"density")
            fprintf("Scaling %s\n",field_name);
            variable_ref_val = variable_ref_val_list{strcmp(variable_ref_val_list(:,1), clean_field_name), 2};
            Interp_deriv_field = model_scaling_factor * Interp_deriv_field * variable_ref_val * V_ref / Q_bar; % Scaling
        else
            fprintf("Scaling %s\n",field_name);
            Interp_deriv_field = model_scaling_factor * Interp_deriv_field * V_ref / Q_bar;
        end
        comb_sensitivities.(clean_field_name) = Interp_deriv_field;
        
        %% Step 6-7: Process nozzle grid (set boundary from combustor and smoothen)
        fprintf('  Processing nozzle grid (setting boundary from combustor and smoothening)...\n');
        
        % Start with original nozzle sensitivity data
        nozzle_field = zeros(size(noz_sensitivities.R2));
        
        % Set boundary values from combustor interpolated field
        % Take boundary values from combustor and apply to nozzle
        nozzle_field(end,:) = Interp_deriv_field(1,1:size(R2,2));
        
        % Apply smoothening to the nozzle field with updated boundary
        temp_f = myutils.f_return_smooth_field(nozzle_field, window, 'col');
        noz_sensitivities.(clean_field_name) = temp_f;
        clear Interp_deriv_field temp_f nozzle_field;
        
        fprintf('  ? Successfully processed %s\n', clean_field_name);
        
    catch ME
        fprintf('  ? Error processing %s: %s\n', field_name, ME.message);
        continue;
    end
end

%% Save results
fprintf('\n=== Saving Results ===\n');

% Create output directory if it doesn't exist
output_dir = './interpolated_sensitivities_output';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('Created output directory: %s\n', output_dir);
end

% Create plots directory
plots_dir = fullfile(output_dir, 'plots');
if ~exist(plots_dir, 'dir')
    mkdir(plots_dir);
    fprintf('Created plots directory: %s\n', plots_dir);
end

% Save combustor sensitivities
comb_output_file = fullfile(output_dir, 'combustor_interpolated_sensitivities.mat');
save(comb_output_file, '-struct', 'comb_sensitivities');
fprintf('? Saved combustor sensitivities to: %s\n', comb_output_file);

% Save nozzle sensitivities
noz_output_file = fullfile(output_dir, 'nozzle_sensitivities.mat');
save(noz_output_file, '-struct', 'noz_sensitivities');
fprintf('? Saved nozzle sensitivities to: %s\n', noz_output_file);

% % Save a summary info file
summary_info = struct();
summary_info.processed_fields = {};
summary_info.combustor_grid_size = size(Z1);
summary_info.nozzle_grid_size = size(Z2);
summary_info.processing_date = datetime('now');
summary_info.Z_restriction_value = 8.5;
summary_info.D = D;

% Collect successfully processed fields
field_names = fieldnames(comb_sensitivities);
excluded_fields = {'R1', 'Z1', 'C_field_MAT', 'Yu', 'Yb'};
summary_info.processed_fields = setdiff(field_names, excluded_fields);
% 
% summary_file = fullfile(output_dir, 'processing_summary.mat');
% save(summary_file, 'summary_info');
% fprintf('? Saved processing summary to: %s\n', summary_file);

%% Generate comparison plots
fprintf('\n=== Generating Visualization Plots ===\n');

% Define sensitivity labels
sensitivity_labels = struct();
sensitivity_labels.Temperature = '$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial T}$';
sensitivity_labels.density = '$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial \rho}$';
sensitivity_labels.CH4 = '$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial Y_{CH4}}$';
sensitivity_labels.O2 = '$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial Y_{O2}}$';
sensitivity_labels.CO2 = '$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial Y_{CO2}}$';
sensitivity_labels.H2O = '$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial Y_{H2O}}$';
sensitivity_labels.N2 = '$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial Y_{N2}}$';

for i = 1:length(summary_info.processed_fields)
    field_name = summary_info.processed_fields{i};
    figidx = 100 + i; % Use unique figure index for each field
    
    % Get appropriate label
    if isfield(sensitivity_labels, field_name)
        label = sensitivity_labels.(field_name);
    else
        label = sprintf('$\\frac{\\partial \\langle \\dot{\\omega}_{T}|c\\rangle}{\\partial %s}$', strrep(field_name, '_', '\\_'));
    end
    
    fprintf('Creating plots for %s (Figure %d)...\n', field_name, figidx);
    try
        plot_sensitivity_comparison(comb_sensitivities, noz_sensitivities, field_name, label, plots_dir, figidx);
        fprintf('  ? Plots saved for %s\n', field_name);
    catch ME
        fprintf('  ? Error creating plots for %s: %s\n', field_name, ME.message);
    end
end

%% Display summary
fprintf('\n=== Processing Summary ===\n');
fprintf('Successfully processed %d sensitivity fields:\n', length(summary_info.processed_fields));
for i = 1:length(summary_info.processed_fields)
    fprintf('  %d. %s\n', i, summary_info.processed_fields{i});
end
fprintf('\nFiles saved to: %s\n', output_dir);
fprintf('Combustor grid size: %dx%d\n', summary_info.combustor_grid_size(1), summary_info.combustor_grid_size(2));
fprintf('Nozzle grid size: %dx%d\n', summary_info.nozzle_grid_size(1), summary_info.nozzle_grid_size(2));

fprintf('\n=== Processing Complete ===\n');
%%
if write_to_h5_file_flag
    N_fields = length(sensitivity_fields);
    inp_fieldNames = strrep(sensitivity_fields,'sensitivity_', '');
    write_field_to_h5_file(inp_fieldNames,N_fields,comb_sensitivities,noz_sensitivities,h5_outdir,h5filename);
end
