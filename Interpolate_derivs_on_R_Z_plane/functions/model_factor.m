function factor = model_factor(global_hrr_mean,volumetric_hrr,R_MAT,Z_MAT)
    r_vec = R_MAT(1,:)';
    z_vec = Z_MAT(:,1)';
    integrand = R_MAT.*volumetric_hrr;
    integrand = trapz(r_vec,integrand,2);
    global_hrr = trapz(z_vec,integrand,1);
    global_hrr = 2*pi*global_hrr;
    factor = global_hrr_mean./global_hrr;
end