function [df_dz] = compute_dfdz(f,z)
    
    [NZ,NR] = size(z);
    df_dz = zeros(NZ,NR);
%     df_dz = zeros(NZ,NR);
    %z_vec(1,:) = z(1,:);
    z_vec(1,:) = z(:,1);
    
    for ii = 1:NR
        temp_vec(1,:) = f(:,ii);
        temp_vec_deriv = compute_deriv(temp_vec,z_vec);
        df_dz(:,ii) = temp_vec_deriv(1,:);
        clear temp_vec_deriv temp_vec;
    end
%     for ii = 1:NR
%         temp_vec(1,:) = f(:,ii);
%         temp_vec_deriv = compute_deriv(temp_vec,z_vec);
%         df_dz(:,ii) = temp_vec_deriv(1,:);
%         clear temp_vec_deriv temp_vec;
%     end
 
end