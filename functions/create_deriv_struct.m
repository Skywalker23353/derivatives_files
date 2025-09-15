function create_deriv_struct(case_dir,deriv_dir,fname,dUr_dr,dUth_dr,dUZ_dr,dUr_dz,dUth_dz,dUZ_dz,R,Z)
    save(sprintf("%s/derivatives_files/%s/%s.mat",case_dir,deriv_dir,fname),"dUr_dr","dUth_dr","dUZ_dr","dUr_dz","dUth_dz","dUZ_dz","R","Z");
    fprintf("%s/derivatives_files/%s/%s.mat CREATED \n",case_dir,deriv_dir,fname);
    clear dUr_dr dUth_dr dUZ_dr dUr_dz dUth_dz dUZ_dz R Z ;
end
