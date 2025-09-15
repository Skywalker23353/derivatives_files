function [Interp_plane_fields] = get_azi_interp_fields(LES_OP_dir,interp_file_dir,Interp_field_names,LES_OP_phasename,LES_OP_Indx,plane_id)
    tic;
    mapping_filename = sprintf('%s/interp_coeffs_%d.dat',...
        interp_file_dir,plane_id);

    [Interp_Fields] = ...
        interpsol_slicer(LES_OP_dir,LES_OP_phasename,LES_OP_Indx,mapping_filename,Interp_field_names);
        Interp_plane_fields = Interp_Fields;
    
    clear Interp_Fields;
    toc;
end


