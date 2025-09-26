clear;clc;close all;
addpath('~/MATLAB/');
addpath('~/satyam_files/CH4_jet_PF/2025_Runs/derivatives_files/Interpolate_derivs_on_R_Z_plane/functions');

%% Configuration
work_dir = '/work/home/satyam/satyam_files/CH4_jet_PF/2025_Runs/derivatives_files/Interpolate_derivs_on_R_Z_plane';
deriv_dir = '/work/home/satyam/satyam_files/CH4_jet_PF/2025_Runs/derivatives_files/sensitivities_10D';
species_deriv_dir = '/work/home/satyam/satyam_files/CH4_jet_PF/2025_Runs/derivatives_files/species_sensitivities_10D';
%% Configuration
% parameters
write_to_h5_file_flag = false;
save_results_flag = true;
h5filename = 'Reactants_1';
h5_outdir = '/work/home/satyam/satyam_files/CH4_jet_PF/2025_Runs/LES_base_case_v6/filtering_run3/src_sensitivities/10D';
D = 2e-3;
window = 3; % Window size for nozzle data smoothening (adjust as needed)
rmx = 5;
zmx=10;
Yu = 0.222606;Yb = 0.041;% Yb = 0.0423208;
Min_c_limit = 1e-3;Max_c_limit = 1.0;
% Yb = 0.039;
l_ref = 2e-3;
U_ref = 65;
V_ref = l_ref^3;
Cp_ref = 1100;
rho_ref = 0.4237;
T_ref = 800;
variable_ref_val_list = {'density',rho_ref; 'Temperature',T_ref;}; 
omega_dot_k_scaling = 1;
omega_dot_T_scaling = (rho_ref*Cp_ref*T_ref*U_ref)/l_ref;
load("comb_interpolted_hrr_field.mat","model_scaling_factor"); %hrr scaling factor
%%
% Define all sensitivity fields to process
sensitivity_fields = {
    'sensitivity_Temperature';
    'sensitivity_density';
    'sensitivity_CH4';
    'sensitivity_O2';
    'sensitivity_CO2';
    'sensitivity_H2O';
};
% Define chemical source sensitivity fields (species derivatives)
chem_src_sensitivity_fields = {
    % CH4 derivatives
    'dCH4_dT';      
    'dCH4_drho';    
    'dCH4_dCH4';    
    'dCH4_dO2';     
    'dCH4_dCO2';    
    'dCH4_dH2O';    
    
    % O2 derivatives  
    'dO2_dT';       
    'dO2_drho';     
    'dO2_dCH4';     
    'dO2_dO2';      
    'dO2_dCO2';     
    'dO2_dH2O';     
    
    % CO2 derivatives
    'dCO2_dT';      
    'dCO2_drho';    
    'dCO2_dCH4';    
    'dCO2_dO2';     
    'dCO2_dCO2';    
    'dCO2_dH2O';    
    
    % H2O derivatives
    'dH2O_dT';      
    'dH2O_drho';    
    'dH2O_dCH4';    
    'dH2O_dO2';     
    'dH2O_dCO2';    
    'dH2O_dH2O';    
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
comb_data.C_field_MAT = (Yu - comb_data.O2_mean)/(Yu - Yb);
noz_data.C_field_MAT = (Yu - noz_data.O2_mean)/(Yu - Yb);

% Load reference C_MAT and Z_MAT for sensitivities
fprintf('Loading reference C and Z matrices...\n');
load('/work/home/satyam/satyam_files/CH4_jet_PF/2025_Runs/c_cond_stats/C_cond_fields_800_10D/CZ_data.mat', "C_MAT", "Z_MAT");
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
        fprintf("Scaling %s\n",field_name);
        
        Interp_deriv_field = model_scaling_factor * Interp_deriv_field  / omega_dot_T_scaling; % Scaling
        if strcmp(clean_field_name,'Temperature')
            Interp_deriv_field = Interp_deriv_field.*T_ref;
        elseif strcmp(clean_field_name,'density')
            Interp_deriv_field = Interp_deriv_field.*rho_ref;
        end
        clean_field_name = sprintf('dw_T_d%s',clean_field_name);
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

%% Step 8: Process Chemical Source Sensitivity Fields
fprintf('\n=== Processing Chemical Source Sensitivity Fields ===\n');

% Process each chemical source sensitivity field
for field_idx = 1:length(chem_src_sensitivity_fields)
    field_name = chem_src_sensitivity_fields{field_idx};
    fprintf('\n--- Processing Chem Src Sensitivity %d/%d: %s ---\n', field_idx, length(chem_src_sensitivity_fields), field_name);
    
    % Load chemical source sensitivity data
    chem_sensitivity_file = sprintf('%s/%s.mat', species_deriv_dir, field_name);
    
    if ~exist(chem_sensitivity_file, 'file')
        fprintf('? Warning: Chemical source sensitivity file not found: %s\n', chem_sensitivity_file);
        continue;
    end
    
    try
        fprintf('  Loading chemical source sensitivity data...\n');
        chem_deriv_data = load(chem_sensitivity_file);
        
        %% Step 5: Interpolate chemical source sensitivity onto combustor grid
        fprintf('  Interpolating chemical source sensitivity onto combustor grid...\n');
        Interp_chem_deriv_field = zeros(size(comb_data.C_field_MAT));
        
        for i = 1:size(Z1, 1)
            if mod(i, 10) == 0 || i == 1 || i == size(Z1, 1)
                fprintf('    Processing Z level %d/%d...\n', i, size(Z1, 1));
            end
            Interp_chem_deriv_field(i, :) = linear_interp_field(comb_data.C_field_MAT(i, :)', ...
                                                              C_MAT(1, :)', ...
                                                              chem_deriv_data.sensitivity(i, :)');
        end
        
        % Store interpolated result in combustor structure
        % Apply same scaling as other sensitivities - no special scaling for species derivatives
        fprintf("Scaling chemical source sensitivity %s\n", field_name);
        Interp_chem_deriv_field = model_scaling_factor * Interp_chem_deriv_field / omega_dot_k_scaling;
%         fieldName = split(field_name, '_');
%         fieldName = fieldName{end};
%         if strcmp(fieldName,'dT')
%             Interp_chem_deriv_field = Interp_chem_deriv_field*T_ref;
%         elseif strcmp(field_name,'drho')
%             Interp_chem_deriv_field = Interp_chem_deriv_field*rho_ref;
%         end
        comb_sensitivities.(field_name) = Interp_chem_deriv_field;
        
        %% Step 6-7: Process nozzle grid (set boundary from combustor and smoothen)
        fprintf('  Processing nozzle grid for chemical source sensitivity...\n');
        
        % Start with original nozzle sensitivity data
        nozzle_chem_field = zeros(size(noz_sensitivities.R2));
        
        % Set boundary values from combustor interpolated field
        % Take boundary values from combustor and apply to nozzle
        nozzle_chem_field(end,:) = Interp_chem_deriv_field(1,1:size(R2,2));
        
        % Apply smoothening to the nozzle field with updated boundary
        temp_chem_f = myutils.f_return_smooth_field(nozzle_chem_field, window, 'col');
        noz_sensitivities.(field_name) = temp_chem_f;
        clear Interp_chem_deriv_field temp_chem_f nozzle_chem_field;
        
        fprintf('  ? Successfully processed chemical source sensitivity: %s\n', field_name);
        
    catch ME
        fprintf('  ? Error processing chemical source sensitivity %s: %s\n', field_name, ME.message);
        continue;
    end
end

%% Save results
fprintf('\n=== Saving Results ===\n');

% Create output directory if it doesn't exist
output_dir = './interpolated_src_terms_sensitivities_output';
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
if save_results_flag
    % Save combustor sensitivities
    comb_output_file = fullfile(output_dir, 'combustor_interpolated_sensitivities.mat');
    save(comb_output_file, '-struct', 'comb_sensitivities');
    fprintf('? Saved combustor sensitivities to: %s\n', comb_output_file);
    
    % Save nozzle sensitivities
    noz_output_file = fullfile(output_dir, 'nozzle_sensitivities.mat');
    save(noz_output_file, '-struct', 'noz_sensitivities');
    fprintf('? Saved nozzle sensitivities to: %s\n', noz_output_file);
end

% Save a summary info file
summary_info = struct();
summary_info.regular_sensitivity_fields = {};
summary_info.chem_src_sensitivity_fields = {};
summary_info.combustor_grid_size = size(Z1);
summary_info.nozzle_grid_size = size(Z2);
summary_info.processing_date = datetime('now');
summary_info.Z_restriction_value = 8.5;
summary_info.D = D;

% Collect successfully processed fields
field_names = fieldnames(comb_sensitivities);
excluded_fields = {'R1', 'Z1', 'C_field_MAT', 'Yu', 'Yb'};
all_processed_fields = setdiff(field_names, excluded_fields);

% Separate regular sensitivities from chemical source sensitivities
regular_processed = {};
chem_src_processed = {};

for i = 1:length(all_processed_fields)
    field = all_processed_fields{i};
    % Check if it's a chemical source sensitivity (starts with 'd' and contains '_d')
    if startsWith(field, 'dw_T')
        regular_processed{end+1} = field;
    elseif startsWith(field, 'd') && contains(field, '_d')
        chem_src_processed{end+1} = field;
    end
end

summary_info.regular_sensitivity_fields = regular_processed;
summary_info.chem_src_sensitivity_fields = chem_src_processed;
summary_info.total_processed_fields = length(all_processed_fields);
% 
% summary_file = fullfile(output_dir, 'processing_summary.mat');
% save(summary_file, 'summary_info');
% fprintf('? Saved processing summary to: %s\n', summary_file);

%% Generate comparison plots
fprintf('\n=== Generating Visualization Plots ===\n');

% Define sensitivity labels for regular sensitivities
sensitivity_labels = struct();
sensitivity_labels.dw_T_dTemperature = '$\frac{\partial  \dot{\omega}_{T}}{\partial T}$';
sensitivity_labels.dw_T_ddensity = '$\frac{\partial  \dot{\omega}_{T}}{\partial \rho}$';
sensitivity_labels.dw_T_dCH4 = '$\frac{\partial  \dot{\omega}_{T}}{\partial Y_{CH4}}$';
sensitivity_labels.dw_T_dO2 = '$\frac{\partial  \dot{\omega}_{T}}{\partial Y_{O2}}$';
sensitivity_labels.dw_T_dCO2 = '$\frac{\partial  \dot{\omega}_{T}}{\partial Y_{CO2}}$';
sensitivity_labels.dw_T_dH2O = '$\frac{\partial  \dot{\omega}_{T}}{\partial Y_{H2O}}$';
sensitivity_labels.dw_T_dN2 = '$\frac{\partial  \dot{\omega}_{T}}{\partial Y_{N2}}$';

% Define labels for chemical source sensitivities
chem_src_labels = struct();
% CH4 derivatives
chem_src_labels.dCH4_dT = '$\frac{\partial  \dot{\omega}_{CH_4}}{\partial T}$';
chem_src_labels.dCH4_drho = '$\frac{\partial  \dot{\omega}_{CH_4}}{\partial \rho}$';
chem_src_labels.dCH4_dCH4 = '$\frac{\partial  \dot{\omega}_{CH_4}}{\partial Y_{CH4}}$';
chem_src_labels.dCH4_dO2 = '$\frac{\partial  \dot{\omega}_{CH_4}}{\partial Y_{O2}}$';
chem_src_labels.dCH4_dCO2 = '$\frac{\partial  \dot{\omega}_{CH_4}}{\partial Y_{CO2}}$';
chem_src_labels.dCH4_dH2O = '$\frac{\partial  \dot{\omega}_{CH_4}}{\partial Y_{H2O}}$';

% O2 derivatives
chem_src_labels.dO2_dT = '$\frac{\partial  \dot{\omega}_{O_2}}{\partial T}$';
chem_src_labels.dO2_drho = '$\frac{\partial  \dot{\omega}_{O_2}}{\partial \rho}$';
chem_src_labels.dO2_dCH4 = '$\frac{\partial  \dot{\omega}_{O_2}}{\partial Y_{CH4}}$';
chem_src_labels.dO2_dO2 = '$\frac{\partial  \dot{\omega}_{O_2}}{\partial Y_{O2}}$';
chem_src_labels.dO2_dCO2 = '$\frac{\partial  \dot{\omega}_{O_2}}{\partial Y_{CO2}}$';
chem_src_labels.dO2_dH2O = '$\frac{\partial  \dot{\omega}_{O_2}}{\partial Y_{H2O}}$';

% CO2 derivatives
chem_src_labels.dCO2_dT = '$\frac{\partial  \dot{\omega}_{CO_2}}{\partial T}$';
chem_src_labels.dCO2_drho = '$\frac{\partial  \dot{\omega}_{CO_2}}{\partial \rho}$';
chem_src_labels.dCO2_dCH4 = '$\frac{\partial  \dot{\omega}_{CO_2}}{\partial Y_{CH4}}$';
chem_src_labels.dCO2_dO2 = '$\frac{\partial  \dot{\omega}_{CO_2}}{\partial Y_{O2}}$';
chem_src_labels.dCO2_dCO2 = '$\frac{\partial  \dot{\omega}_{CO_2}}{\partial Y_{CO2}}$';
chem_src_labels.dCO2_dH2O = '$\frac{\partial  \dot{\omega}_{CO_2}}{\partial Y_{H2O}}$';

% H2O derivatives
chem_src_labels.dH2O_dT = '$\frac{\partial  \dot{\omega}_{H_2O}}{\partial T}$';
chem_src_labels.dH2O_drho = '$\frac{\partial  \dot{\omega}_{H_2O}}{\partial \rho}$';
chem_src_labels.dH2O_dCH4 = '$\frac{\partial  \dot{\omega}_{H_2O}}{\partial Y_{CH4}}$';
chem_src_labels.dH2O_dO2 = '$\frac{\partial  \dot{\omega}_{H_2O}}{\partial Y_{O2}}$';
chem_src_labels.dH2O_dCO2 = '$\frac{\partial  \dot{\omega}_{H_2O}}{\partial Y_{CO2}}$';
chem_src_labels.dH2O_dH2O = '$\frac{\partial  \dot{\omega}_{H_2O}}{\partial Y_{H2O}}$';

% Generate plots for regular sensitivities
fprintf('Generating plots for regular sensitivities...\n');
for i = 1:length(summary_info.regular_sensitivity_fields)
    field_name = summary_info.regular_sensitivity_fields{i};
    figidx = 100 + i; % Use unique figure index for each field
    
    % Get appropriate label
    if isfield(sensitivity_labels, field_name)
        label = sensitivity_labels.(field_name);
    else
        label = sprintf('$\\frac{\\partial \ \\dot{\\omega}_{T}\}{\\partial %s}$', strrep(field_name, '_', '\\_'));
    end
    
    fprintf('Creating plots for regular sensitivity %s (Figure %d)...\n', field_name, figidx);
    try
        plot_sensitivity_comparison(comb_sensitivities, noz_sensitivities, field_name, label, plots_dir, figidx);
        fprintf('  ? Plots saved for regular sensitivity %s\n', field_name);
    catch ME
        fprintf('  ? Error creating plots for regular sensitivity %s: %s\n', field_name, ME.message);
    end
end

% Generate plots for chemical source sensitivities
fprintf('Generating plots for chemical source sensitivities...\n');
for i = 1:length(summary_info.chem_src_sensitivity_fields)
    field_name = summary_info.chem_src_sensitivity_fields{i};
    figidx = 500 + i; % Use unique figure index range for chemical source sensitivities
    
    % Get appropriate label
    if isfield(chem_src_labels, field_name)
        label = chem_src_labels.(field_name);
    else
        % Generic label for any missed chemical source sensitivity
        label = sprintf('$\\frac{\\partial \ \\dot{\\omega}_{species}\}{\\partial field}$ (%s)', strrep(field_name, '_', '\\_'));
    end
    
    fprintf('Creating plots for chemical source sensitivity %s (Figure %d)...\n', field_name, figidx);
    try
        plot_sensitivity_comparison(comb_sensitivities, noz_sensitivities, field_name, label, plots_dir, figidx);
        fprintf('  ? Plots saved for chemical source sensitivity %s\n', field_name);
    catch ME
        fprintf('  ? Error creating plots for chemical source sensitivity %s: %s\n', field_name, ME.message);
    end
end

%% Display summary
fprintf('\n=== Processing Summary ===\n');
fprintf('Successfully processed %d total fields:\n', summary_info.total_processed_fields);

fprintf('\nRegular sensitivity fields (%d):\n', length(summary_info.regular_sensitivity_fields));
for i = 1:length(summary_info.regular_sensitivity_fields)
    fprintf('  %d. %s\n', i, summary_info.regular_sensitivity_fields{i});
end

fprintf('\nChemical source sensitivity fields (%d):\n', length(summary_info.chem_src_sensitivity_fields));
for i = 1:length(summary_info.chem_src_sensitivity_fields)
    fprintf('  %d. %s\n', i, summary_info.chem_src_sensitivity_fields{i});
end

fprintf('\nFiles saved to: %s\n', output_dir);
fprintf('Combustor grid size: %dx%d\n', summary_info.combustor_grid_size(1), summary_info.combustor_grid_size(2));
fprintf('Nozzle grid size: %dx%d\n', summary_info.nozzle_grid_size(1), summary_info.nozzle_grid_size(2));

fprintf('\n=== Processing Complete ===\n');
%%
if write_to_h5_file_flag
    fprintf('\n=== Writing All Sensitivities to H5 File ===\n');
    
    % Combine all field names (regular + chemical source sensitivities)
    all_field_names = [summary_info.regular_sensitivity_fields, summary_info.chem_src_sensitivity_fields];
    N_fields = length(all_field_names);
    
    % Display what will be written
    fprintf('Writing all %d sensitivity fields together to H5 file: %s.h5\n', N_fields, h5filename);
    fprintf('Output directory: %s\n', h5_outdir);
    
    fprintf('\nFields being written (%d regular + %d chemical source = %d total):\n', ...
            length(summary_info.regular_sensitivity_fields), ...
            length(summary_info.chem_src_sensitivity_fields), ...
            N_fields);
    
    % Show regular sensitivities
    fprintf('  Regular sensitivities:\n');
    for i = 1:length(summary_info.regular_sensitivity_fields)
        fprintf('    %d. %s\n', i, summary_info.regular_sensitivity_fields{i});
    end
    
    % Show chemical source sensitivities
    fprintf('  Chemical source sensitivities:\n');
    for i = 1:length(summary_info.chem_src_sensitivity_fields)
        fprintf('    %d. %s\n', i, summary_info.chem_src_sensitivity_fields{i});
    end
    
    % Write all fields together in single operation
    fprintf('\nWriting all sensitivities together to H5 file...\n');
    try
        write_field_to_h5_file(all_field_names, N_fields, comb_sensitivities, noz_sensitivities, h5_outdir, h5filename);
        fprintf('? Successfully wrote all %d sensitivity fields to H5 file: %s.h5\n', N_fields, h5filename);
    catch ME
        fprintf('? Error writing to H5 file: %s\n', ME.message);
    end
    
    fprintf('=== H5 File Writing Complete ===\n');
end
%%
% myutils.plot_field(1,comb_sensitivities.R1,comb_sensitivities.Z1/D,comb_data.C_field_MAT,'$C$');