clear all;
clc;
addpath(genpath('/work/home/satyam/Multisolv_June_intel/tools/slicer'));
addpath("~/MATLAB/");
addpath('~/satyam_files/CH4_jet_PF/2025_Runs/derivatives_files/Interpolate_derivs_on_R_Z_plane/functions');
tic
%% USER INPUTS
Interp_field_names = {{'O2_mean'}};
LES_data_dir = "/work/home/satyam/satyam_files/CH4_jet_PF/2025_Runs/LES_base_case_v6/" + ...
    "filtering_run3/Final_baseflow_with_EBU_components";
work_dir = pwd;
grids = {'combustor_domain','nozzle_domain_with_zero'};
%output file index
Ns = 1011260202;
Ne = Ns;
nplanes= 16;
centerline_avg = 0;
%phasename
LES_OP_phasename = 'Reactants';
%LES mat file name

mat_file_name = sprintf("YO2_field_structured_grid_%d_planes_%d",nplanes,Ns);
load(sprintf('%s/structured_grid_from_LES_grid_with_zero_%d_unilat_planes.mat',work_dir,nplanes));

size_list = [size(R1);size(R2)];
%%
for ii = 1:length(grids)
    %path to folder with slicer outputs/interpolation coeffs
    interp_coeff_dir = sprintf('%s/interp_%d_unilat_structured_pl_%s',work_dir,nplanes,grids{ii});
    OP_dir = sprintf('%s/interp_%d_unilat_structured_pl_%s/interp_fields',work_dir,nplanes,grids{ii});
    mkdir (OP_dir);

    NI = size_list(ii,1);
    NJ = size_list(ii,2);
    f_create_2D_slice_fields(LES_data_dir,interp_coeff_dir,OP_dir,mat_file_name,...
        Interp_field_names,LES_OP_phasename,Ns,Ne,nplanes,NI,NJ,centerline_avg);%% Always check this function
    clear NI NJ;
end
%%
OP_dir = sprintf('%s/interp_%d_unilat_structured_pl_%s/interp_fields',work_dir,nplanes,grids{1});
comb = load(sprintf('%s/%s_avg.mat',OP_dir,mat_file_name));
O2_fmean = sum(comb.O2_fmean,3)/size(comb.O2_fmean,3);
% save(sprintf("%s/%s_avg.mat",OP_dir,mat_file_name),"O2_fmean");
clear O2_fmean;
%%

OP_dir = sprintf('%s/interp_%d_unilat_structured_pl_%s/interp_fields',work_dir,nplanes,grids{2});
noz = load(sprintf('%s/%s_avg.mat',OP_dir,mat_file_name));
O2_fmean = sum(noz.O2_fmean,3)/size(noz.O2_fmean,3);
% save(sprintf("%s/%s_avg.mat",OP_dir,mat_file_name),"O2_fmean");
%%

% figIdx = 10;
% figure(figIdx)
% % myutils.plot_field(figIdx,R1,Z1,sum(comb.O2_fmean,3)/size(comb.O2_fmean,3),'$Y_{O_2}$');
% myutils.plot_field(figIdx,R1,Z1,comb.O2_fmean,'$Y_{O_2}$');
% pbaspect([9 16 1]);
% hold on ;
% myutils.plot_field(figIdx,R2,Z2,noz.O2_fmean,'$Y_{O_2}$');
% pbaspect([9 16 1]);
% hold off ;
%%
