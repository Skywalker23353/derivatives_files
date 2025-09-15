clear;
%%
tic
work_dir = pwd;
case_name = "16_unilat_structured_pl_combustor_domain";

N_planes = 16;

 for ii = 1:1:N_planes
    plane_id = ii;
    coord_file = load(sprintf("%s/coords_%s/coords_pl_%d.dat",work_dir,case_name,plane_id));
    interp_file = readtable(sprintf("%s/interp_%s/interp_coeff_%d.dat",work_dir,case_name,plane_id));
    [nrc,ncc] = size(coord_file);
    [nri,nci] = size(interp_file);
    a(ii,1) = nrc - nri;
    clear coord_file interp_file;
 end
max(abs(a))
toc

% O means all points are interpolated.
% NON zero output mean some points are not interpolated.