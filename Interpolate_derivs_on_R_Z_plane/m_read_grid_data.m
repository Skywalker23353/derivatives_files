clear;clc;close all;

Phase_Name = 'Reactants';

% file_dir = ...
%     '/work/home/satyam/satyam_files/LES_base_case';

file_dir = '/work/home/satyam/satyam_files/Burner_6x_LSA';
src_Indx = 5501;
F_src_name = sprintf('%s/%s_%d.h5',file_dir,Phase_Name,src_Indx);

Grid_Fname = sprintf('%s/%s_grids_0.h5',file_dir,Phase_Name);
Grid_Info = h5info(Grid_Fname);
Num_grids = length(Grid_Info.Groups);
num_src_blks_grid = zeros(Num_grids,1);
for ii=1:Num_grids
    num_src_blks_grid(ii) = length(Grid_Info.Groups(ii).Groups.Groups);
    grid_name{ii,1} = Grid_Info.Groups(ii).Name;
end