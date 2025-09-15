function [df_dx] = compute_deriv(f,x)

    eps_tol_deriv = 1e-8;
    dx_deta = compute_deriv_8th_order_stencil(x);
    
    df_deta = compute_deriv_8th_order_stencil(f);
    df_dx = df_deta./dx_deta;

    for ii=1:length(dx_deta)
        if abs(dx_deta(ii)) < eps_tol_deriv 
            df_dx(ii) = 0.0;
        end
    end

end