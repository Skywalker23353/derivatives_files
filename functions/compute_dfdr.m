function [df_dr] = compute_dfdr(f,r)
    
    [NZ,NR] = size(r);
    df_dr = zeros(NZ,NR);
    r_vec(1,:) = r(1,:);

    for ii = 1:NZ
        temp_vec(1,:) = f(ii,:);
        temp_vec_deriv = compute_deriv(temp_vec,r_vec);
        df_dr(ii,:) = temp_vec_deriv(1,:);
        clear temp_vec_deriv temp_vec;
    end
 
end