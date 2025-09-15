function write_field_to_h5_file(fieldNames,N_fields,comb,noz,outdir,filename)

    h5_filename = sprintf('%s/%s_0.h5',outdir,filename);
    h5_grids_filename = sprintf('%s/%s_grids_0.h5',outdir,filename);
    phase_name = filename;
    grid_name1 = 'Combustor';
    grid_name2 = 'Nozzle';
    
    for i = 1:N_fields
        fieldName = fieldNames{i};
        fprintf("Writing %s \n",fieldName);
        
        comb_field_rh = comb.(fieldName);
        comb_field_lh = flip(comb_field_rh,2);

        noz_field_rh = noz.(fieldName);
        noz_field_lh = flip(noz_field_rh,2);
    
        comb_field = cat(2,comb_field_lh(:,1:end-1),comb_field_rh);
        noz_field = cat(2,noz_field_lh(:,1:end-1),noz_field_rh);
    
        %% % Increasing the radial size of nozzle domain back to original
        noz.R2(:,end) = noz.R2(:,end) + 3e-06;
        noz.R2(:,1) = noz.R2(:,1) - 3e-06;
        %% Creating 3D domain 
        comb_R_lh = -flip(comb.R1,2);
        comb_bilat = cat(2,comb_R_lh(:,1:end-1),comb.R1);
        [nrc,ncc] = size(comb_bilat);
        x_comb = zeros(nrc,ncc,2);
        x_comb(:,:,1) = comb_bilat;
        x_comb(:,:,2) = comb_bilat;
        y_comb = zeros(nrc,ncc,2);
        comb_Z = cat(2,comb.Z1(:,1:end-1),comb.Z1);  
        y_comb(:,:,1) = comb_Z;
        y_comb(:,:,2) = comb_Z;
        z_comb = zeros(nrc,ncc,2);
        z_comb(:,:,1) = ones(nrc,ncc)*1e-3;
        z_comb(:,:,2) = -1e-3*ones(nrc,ncc);
        field_comb_3D = zeros(nrc,ncc,2);
        field_comb_3D(:,:,1) = comb_field;
        field_comb_3D(:,:,2) = comb_field;
        %
        noz_R_lh = -flip(noz.R2,2);
        noz_bilat = cat(2,noz_R_lh(:,1:end-1),noz.R2);
        [nrn,ncn] = size(noz_bilat);
        x_noz = zeros(nrn,ncn,2);
        x_noz(:,:,1) = noz_bilat;
        x_noz(:,:,2) = noz_bilat;
        y_noz = zeros(nrn,ncn,2);
        noz_Z = cat(2,noz.Z2(:,1:end-1),noz.Z2);  
        y_noz(:,:,1) = noz_Z;
        y_noz(:,:,2) = noz_Z;
        z_noz = zeros(nrn,ncn,2);
        z_noz(:,:,1) = ones(nrn,ncn)*1e-3;
        z_noz(:,:,2) = -1e-3*ones(nrn,ncn);
        
        field_noz_3D = zeros(nrn,ncn,2);
        field_noz_3D(:,:,1) = noz_field;
        field_noz_3D(:,:,2) = noz_field;
        %% Writing

        if ~isfolder(outdir);mkdir(outdir);end

        if i==1 % writing grids file only once
            x_comb_path = sprintf('/%s/source_blocks/0/x',grid_name1);
            h5create(h5_grids_filename,x_comb_path,size(x_comb));
            y_comb_path = sprintf('/%s/source_blocks/0/y',grid_name1);
            h5create(h5_grids_filename,y_comb_path,size(y_comb));
            z_comb_path = sprintf('/%s/source_blocks/0/z',grid_name1);
            h5create(h5_grids_filename,z_comb_path,size(z_comb));
            h5write(h5_grids_filename,x_comb_path,x_comb);
            h5write(h5_grids_filename,y_comb_path,y_comb);
            h5write(h5_grids_filename,z_comb_path,z_comb);

            x_noz_path = sprintf('/%s//source_blocks/0/x',grid_name2);
            h5create(h5_grids_filename,x_noz_path,size(x_noz));
            y_noz_path = sprintf('/%s/source_blocks/0/y',grid_name2);
            h5create(h5_grids_filename,y_noz_path,size(y_noz));
            z_noz_path = sprintf('/%s/source_blocks/0/z',grid_name2);
            h5create(h5_grids_filename,z_noz_path,size(z_noz));
            h5write(h5_grids_filename,x_noz_path,x_noz);
            h5write(h5_grids_filename,y_noz_path,y_noz);
            h5write(h5_grids_filename,z_noz_path,z_noz);
        end

        field_path_c = sprintf('/%s/%s/fields/0/%s',phase_name,grid_name1,fieldName);
        h5create(h5_filename,field_path_c,size(field_comb_3D));
        h5write(h5_filename,field_path_c,field_comb_3D);
        
        
        
        field_path_n = sprintf('/%s/%s/fields/0/%s',phase_name,grid_name2,fieldName);
        h5create(h5_filename,field_path_n,size(field_noz_3D));
        h5write(h5_filename,field_path_n,field_noz_3D);
    end
        % h5close();
end