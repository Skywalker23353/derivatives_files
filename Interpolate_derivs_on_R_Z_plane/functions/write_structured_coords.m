function write_structured_coords(R,Z,coords_dir,no_of_planes)
    
    dth = 180/no_of_planes;
    theta = 0:dth:180-dth;
    
    for tt = 1:length(theta)
        X_MAT = R*cosd(theta(tt));
        Y_MAT = R*sind(theta(tt));

        fname_dat = sprintf("./%s/coords_pl_%d.dat",coords_dir,tt);
        fid_dat = fopen(fname_dat,'w');
        
        fname_csv = sprintf("./%s/coords_pl_%d.csv",coords_dir,tt);
        fid_csv = fopen(fname_csv,'w');

        for ii = 1:size(R,1)
            for jj = 1:size(R,2)
                fprintf(fid_dat,'%e,%e,%e\n',X_MAT(ii,jj),Y_MAT(ii,jj),Z(ii,jj));
                fprintf(fid_csv,'%e,%e,%e\n',X_MAT(ii,jj),Y_MAT(ii,jj),Z(ii,jj));
            end
        end
        fclose(fid_dat);
        fclose(fid_csv);
        
        fprintf("Coordinates written for plane %d \n",tt);
        clear xvec yvec;
    end
   
 end
