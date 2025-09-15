clear; clc;%close all;
addpath('~/MATLAB/');
addpath('~/satyam_files/CH4_jet_PF/2025_Runs/derivatives_files/Interpolate_derivs_on_R_Z_plane/functions');
%% Configuration
work_dir = '/work/home/satyam/satyam_files/CH4_jet_PF/2025_Runs/derivatives_files/Interpolate_derivs_on_R_Z_plane';
D = 2e-3;
global_hrr_mean = 163;
save_data = 0;

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
Yu = 0.222606; %Yb = 0.0423208;
Yb = 0.0399;
comb_data.C_field_MAT = (Yu - comb_data.O2_fmean)/(Yu - Yb);
noz_data.C_field_MAT = (Yu - noz_data.O2_fmean)/(Yu - Yb);

% Load hrr data
hrr = load('/work/home/satyam/satyam_files/CH4_jet_PF/2025_Runs/c_cond_stats/C_cond_fields_800/Heatrelease_smooth.mat');

% Apply Z restriction logic
Z_idx_mx = find((hrr.Z_MAT)/D >= 8.5, 1);
r_idx_mx = find((R1(1,:))/D >= 5, 1);
R1 = R1(1:Z_idx_mx, 1:r_idx_mx); 
Z1 = Z1(1:Z_idx_mx, 1:r_idx_mx);
comb_data.C_field_MAT = comb_data.C_field_MAT(1:Z_idx_mx, 1:r_idx_mx);

% Sanity check for C field
% comb_data.C_field_MAT(find(comb_data.C_field_MAT < 1e-5)) = 0;
comb_data.C_field_MAT(find(comb_data.C_field_MAT > 1)) = 1;

fprintf('Setup complete. Grid size: %dx%d\n', size(Z1, 1), size(Z1, 2));
%%
% Initialize results structures
comb = struct();
noz = struct();

% Copy coordinate and C field data to results
comb.R1 = R1;
comb.Z1 = Z1;
comb.C_field_MAT = comb_data.C_field_MAT;
comb.Yu = Yu;
comb.Yb = Yb;

noz.R2 = R2;
noz.Z2 = Z2;
noz.C_field_MAT = noz_data.C_field_MAT;
noz.Yu = Yu;
noz.Yb = Yb;
    

%% Step 5: Interpolate field onto combustor grid
fprintf('  Interpolating onto combustor grid...\n');
Interp_field = zeros(size(comb_data.C_field_MAT));

for i = 1:size(Z1, 1)
    if mod(i, 10) == 0 || i == 1 || i == size(Z1, 1)
        fprintf('    Processing Z level %d/%d...\n', i, size(Z1, 1));
    end
    Interp_field(i, :) = linear_interp_field(comb_data.C_field_MAT(i, :)', ...
                                                  hrr.C_MAT(1, :)', ...
                                                  hrr.DF(i, :)');
end

% Store interpolated result in combustor structure
comb.hrr_field = Interp_field;
    
%% Step 6-7: Process nozzle grid (set boundary from combustor and smoothen)
fprintf('  Processing nozzle grid (setting boundary from combustor and smoothening)...\n');

% Start with original nozzle sensitivity data
nozzle_field = zeros(size(noz.R2));

% Set boundary values from combustor interpolated field
% Take boundary values from combustor and apply to nozzle
nozzle_field(end,:) = Interp_field(1,1:size(R2,2));

% Apply smoothening to the nozzle field with updated boundary
temp_f = myutils.f_return_smooth_field(nozzle_field, 3, 'col');
noz.hrr_field = temp_f;
clear Interp_field temp_f nozzle_field;

fprintf('  ? Successfully processed hrr field\n');
    
%%
figure(1)
myutils.plot_field(1,comb.R1/D,comb.Z1/D,comb.hrr_field,'hrr');hold on;
myutils.plot_field(1,noz.R2/D,noz.Z2/D,noz.hrr_field,'hrr');
pbaspect([1 2 1]);
hold off;
%%
model_scaling_factor = model_factor(global_hrr_mean,comb.hrr_field,comb.R1,comb.Z1);
comb.model_scaling_factor = model_scaling_factor;
if save_data;save('comb_interpolted_hrr_field.mat','-struct','comb');end