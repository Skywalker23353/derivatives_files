function field_mat = get_azi_slice(LES_OP_dir,...
        interp_file_dir,Interp_field_name,LES_OP_phasename,idx,NI,NJ,plane_id)

[Interp_Azi_fields] = get_azi_interp_fields(LES_OP_dir,...
        interp_file_dir,Interp_field_name,LES_OP_phasename,idx,plane_id);
field_mat = get_field_mapped_to_grid(NI,NJ,Interp_Azi_fields);
end