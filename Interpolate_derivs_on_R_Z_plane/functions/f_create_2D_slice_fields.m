function [] = f_create_2D_slice_fields(LES_data_dir,Interp_coeff_dir,OP_dir, ...
    mat_file_name,field_names,LES_OP_phasename, ...
    Ns,Ne,nplanes,NI,NJ, ...
    centerline_avg)
    LES_field_mat = cell(1,nplanes);
    NJ_mid = (NJ+1)/2;
    tempst = struct();
    for ii = Ns:Ne
        for ij=1:length(field_names)
                field_name = field_names{ij};
                for jj = 1:nplanes
                    LES_field_mat{jj} = get_azi_slice(LES_data_dir,...
                        Interp_coeff_dir,field_name,LES_OP_phasename,ii,NI,NJ,jj);
                end
                mat_3D = cat(3,LES_field_mat{:}); 
                tempst.(field_name{1}) = mat_3D;
                if centerline_avg == 1
                 % This section is to remove centerline values and replace with avg values:
                    tempst.(field_name{1})(:,NJ_mid,:) = (tempst.(field_name{1})(:,(NJ_mid-1),:) + tempst.(field_name{1})(:,(NJ_mid+1),:))/2;
                end
                clear mat_3D ;
                clear LES_field_mat;
        end
    end
    save(sprintf('%s/%s.mat',OP_dir,mat_file_name),'-struct',"tempst");
    clear tempst;
end