clear;clc;tic;
% %%
addpath('~/satyam_files/CH4_jet_PF/2025_Runs/derivatives_files/Interpolate_derivs_on_R_Z_plane/functions');
%%%%%%%
%To read the LES grid coordinates 
%Generate 16 unilateral grids 
%Write to dat files 
%%%%%%%
no_of_planes = 16;
case_path = '/work/home/satyam/satyam_files/LES_base_case_v4';
grid_name = 'burner';
src_blk_id = [13,12,6,1];% Initial plane LES source_block no's. i.e Initial plane lies in these blocks.
LES_grids_file_name = 'Reactants_grids_0.h5';
LES_grids_file_path = sprintf('%s/%s',case_path,LES_grids_file_name);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataset_name = 'x';
dataset_path = sprintf('/%s/source_blocks/%d/%s',grid_name,src_blk_id(1),dataset_name);
x1 = h5read(LES_grids_file_path,dataset_path);
x1= permute(x1,[2 1 3]);
dataset_path = sprintf('/%s/source_blocks/%d/%s',grid_name,src_blk_id(2),dataset_name);
x2 = h5read(LES_grids_file_path,dataset_path);
x2= permute(x2,[2 1 3]);
dataset_path = sprintf('/%s/source_blocks/%d/%s',grid_name,src_blk_id(3),dataset_name);
x3 = h5read(LES_grids_file_path,dataset_path);
x3= permute(x3,[2 1 3]);
dataset_path = sprintf('/%s/source_blocks/%d/%s',grid_name,src_blk_id(4),dataset_name);
x4 = h5read(LES_grids_file_path,dataset_path);
x4= permute(x4,[2 1 3]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataset_name = 'y';
dataset_path = sprintf('/%s/source_blocks/%d/%s',grid_name,src_blk_id(1),dataset_name);
y1 = h5read(LES_grids_file_path,dataset_path);
y1= permute(y1,[2 1 3]);
dataset_path = sprintf('/%s/source_blocks/%d/%s',grid_name,src_blk_id(2),dataset_name);
y2 = h5read(LES_grids_file_path,dataset_path);
y2= permute(y2,[2 1 3]);
dataset_path = sprintf('/%s/source_blocks/%d/%s',grid_name,src_blk_id(3),dataset_name);
y3 = h5read(LES_grids_file_path,dataset_path);
y3= permute(y3,[2 1 3]);
dataset_path = sprintf('/%s/source_blocks/%d/%s',grid_name,src_blk_id(4),dataset_name);
y4 = h5read(LES_grids_file_path,dataset_path);
y4= permute(y4,[2 1 3]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataset_name = 'z';
dataset_path = sprintf('/%s/source_blocks/%d/%s',grid_name,src_blk_id(1),dataset_name);
z1 = h5read(LES_grids_file_path,dataset_path);
z1= permute(z1,[2 1 3]);
dataset_path = sprintf('/%s/source_blocks/%d/%s',grid_name,src_blk_id(2),dataset_name);
z2 = h5read(LES_grids_file_path,dataset_path);
z2= permute(z2,[2 1 3]);
dataset_path = sprintf('/%s/source_blocks/%d/%s',grid_name,src_blk_id(3),dataset_name);
z3 = h5read(LES_grids_file_path,dataset_path);
z3= permute(z3,[2 1 3]);
dataset_path = sprintf('/%s/source_blocks/%d/%s',grid_name,src_blk_id(4),dataset_name);
z4 = h5read(LES_grids_file_path,dataset_path);
z4= permute(z4,[2 1 3]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Taking the first index of the blocks to form the LES domain plane
x1 = squeeze(x1(7,:,2:end));
y1 = squeeze(y1(7,:,2:end));
z1= squeeze(z1(7,:,2:end));
%r1 = sqrt(x1.^2 + y1.^2);
r1 = x1(1,:)';

x2 = squeeze(x2(7,:,2:end));
y2 = squeeze(y2(7,:,2:end));
z2= squeeze(z2(7,:,2:end));
% r2 = sqrt(x2.^2 + y2.^2);
r2 = x2(1,:)';
x3 = squeeze(x3(7,:,:));
y3 = squeeze(y3(7,:,:));
z3= squeeze(z3(7,:,:));
% r3 = sqrt(x3.^2 + y3.^2);
r3 = x3(1,:)';
%Taking the last index of block id : 3 to form the nozzle domain plane
x4 = squeeze(x4(7,:,:));
y4 = squeeze(y4(7,:,:));
z4= squeeze(z4(7,:,:));
% r4 = sqrt(x4.^2 + y4.^2);
r4 = x4(1,:)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Creating list of points form the inner circular edge to the centerline.
spacing = (r4(end-1) - r4(end));
%Nozzle domain(ommiting the boundary points from nozzle domain)
% c_r = (0:spacing:r4(1,1));
% [r5 , z5] = meshgrid(c_r, z4(1:(end-1),1));
r5 = (0:spacing:r4(end))';
r5 = sort(r5,'descend');
%nozzle domain w/o zero
%z5 = z4(1:(end-1),1)';
%nozzle domain with zero
z5 = z4(:,1)';
%LES domain(retaining the boundary points lying on the boundary of LES and
%nozzle domian)
% [r6, z6] = meshgrid(c_r,z3(:,1));
%r6 = (0:spacing:r4(end))';   
z6 = z3(:,1)';
%%%%%%%%%%%%%%
% Aggregating coordinates and writing to dat file
R_vec_uni_1 = [r1;r2;r3;r5];
R_vec_uni_2 = [r4;r5];

R_vec_uni_1 = sort(R_vec_uni_1,'ascend');
R_vec_uni_2 = sort(R_vec_uni_2,'ascend');

R_vec_uni_2(end) = 0.997*R_vec_uni_2(end);

[~,c] = min(abs(R_vec_uni_1-R_vec_uni_2(end)));
R_vec_uni_1(1:c,1) = R_vec_uni_2(:,1);
R_vec_uni_1(c+1) = [];

%%
R_vec_uni_1 = R_vec_uni_1';
R_vec_uni_2 = R_vec_uni_2';
% R_vec_bilat_1 = [-flip(R_vec_uni_1(1,2:end)),R_vec_uni_1(1,:)];
% R_vec_bilat_2 = [-flip(R_vec_uni_2(1,2:end)),R_vec_uni_2(1,:)];
%To check wether there are any repeated values in radial direction.
% if length(R_vec_bilat_1) == length(unique(R_vec_bilat_1))
%     fprintf('R_vec_1 has no repeated elements.\n');
% else
%     fprintf('R_vec_1 has repeated elements.\n');
% end
% 
% if length(R_vec_bilat_2) == length(unique(R_vec_bilat_2))
%     fprintf('R_vec_2 has no repeated elements.\n');
% else
%     fprintf('R_vec_2 has repeated elements.\n');
% end
%%
Z_vec_1 = z6;
Z_vec_2 = z5;

Z_vec_2(1) = 0.999*Z_vec_2(1);

% [R1,Z1] = meshgrid(R_vec_bilat_1, Z_vec_1);
% [R2,Z2] = meshgrid(R_vec_bilat_2, Z_vec_2);

[R1,Z1] = meshgrid(R_vec_uni_1, Z_vec_1);
[R2,Z2] = meshgrid(R_vec_uni_2, Z_vec_2);

%save('./structured_grid_from_LES_grid.mat',"R1","Z1","R2","Z2");
save(sprintf('./structured_grid_from_LES_grid_with_zero_%d_unilat_planes.mat',no_of_planes),"R1","Z1","R2","Z2");
% Writing the coordinates to dat files.


% coords_dir = 'coords_8_bilat_structured_pl_nozzle_domain';
% mkdir (coords_dir);
% write_structured_coords_bilat(R2,Z2,coords_dir,no_of_planes);
coords_dir = sprintf("coords_%d_unilat_structured_pl_nozzle_domain_with_zero",no_of_planes);
mkdir (coords_dir);
write_structured_coords(R2,Z2,coords_dir,no_of_planes);

coords_dir = sprintf("coords_%d_unilat_structured_pl_combustor_domain",no_of_planes);
mkdir (coords_dir);
write_structured_coords(R1,Z1,coords_dir,no_of_planes);

fprintf("LES coordinates read \n created mat file \n and written to dat files \n");
toc;
% clear all;
