clear;clc;
addpath('~/MATLAB/');
addpath('~/satyam_files/CH4_jet_PF/2025_Runs/derivatives_files/Interpolate_derivs_on_R_Z_plane/functions');
%%
work_dir = '/work/home/satyam/satyam_files/CH4_jet_PF/2025_Runs/derivatives_files/Interpolate_derivs_on_R_Z_plane';
%load deriv field
deriv_dir = '/work/home/satyam/satyam_files/CH4_jet_PF/2025_Runs/derivatives_files/sensitivities';
field = 'sensitivity_Temperature';
deriv_data = load(sprintf('%s/%s.mat',deriv_dir,field));
%load YO2 field
case_dir = 'interp_16_unilat_structured_pl_combustor_domain';
comb_data = load(sprintf('%s/%s/interp_fields/YO2_field_structured_grid_16_planes_1011260202_avg.mat',work_dir,case_dir));
case_dir = 'interp_16_unilat_structured_pl_nozzle_domain_with_zero';
noz_data = load(sprintf('%s/%s/interp_fields/YO2_field_structured_grid_16_planes_1011260202_avg.mat',work_dir,case_dir));
load('structured_grid_from_LES_grid_with_zero_16_unilat_planes.mat');
load('/work/home/satyam/satyam_files/CH4_jet_PF/2025_Runs/c_cond_stats/C_cond_fields_800/Heatrelease_smooth.mat',"C_MAT","Z_MAT");
D = 2e-3;
Z_idx_mx = find((Z_MAT)/D >= 8.5,1);
R1 = R1(1:Z_idx_mx,:); Z1 = Z1(1:Z_idx_mx,:);

%create C field using YO2 field
Yu = 0.222606; Yb = 0.0423208;
Y = [Yu,Yb]; 
comb_data.C_field_MAT = (Yu - comb_data.O2_fmean)/(Yu - Yb);
noz_data.C_field_MAT = (Yu - noz_data.O2_fmean)/(Yu - Yb);

comb_data.C_field_MAT = comb_data.C_field_MAT(1:Z_idx_mx,:);
% noz_data.C_field_MAT = (Yu - noz_data.O2_fmean)/(Yu - Yb);

% comb_data.C_field_MAT(find(comb_data.C_field_MAT < 1e-3)) = 0;
% comb_data.C_field_MAT(find(comb_data.C_field_MAT > 1)) = 1;
%%
Interp_deriv_field = zeros(size(comb_data.C_field_MAT));
for i = 1:size(Z1,1)
    Z_idx = i
    Interp_deriv_field(i,:) = linear_interp_field(comb_data.C_field_MAT(i,:),C_MAT(1,:),deriv_data.sensitivity(i,:));
end
comb_data.deriv_field = Interp_deriv_field;clear Interp_deriv_field;
Interp_deriv_field = zeros(size(R2));
Interp_deriv_field(end,:) = comb_data.deriv_field(1,1:size(R2,2));
window = 3;
temp_f = myutils.f_return_smooth_field(Interp_deriv_field, window, 'col');
noz_data.deriv_field = temp_f;clear Interp_deriv_field temp_f;
%%
figidx = 2;
figure(figidx)
myutils.plot_field(figidx,R1/D,Z1/D,comb_data.deriv_field,'$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial T}$');hold on ;
myutils.plot_field(figidx,R2/D,Z2/D,noz_data.deriv_field,'$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial T}$');
xlim([0 3]);
pbaspect([1 3 1]);
hold off;

%% write to h5 file

% zt1 = Z2(1,1)
% zt2 = Z2(end,1)
% 
% z = zt1:1e-5:zt2;
% 
% % c = 1 - (zt2 - z)/(zt2 - zt1);
% c = tanh(2*(z - zt1)/(zt2 - zt1));
% figure(101);
% plot(z,c,'+-k','DisplayName','c');
% grid on ;
