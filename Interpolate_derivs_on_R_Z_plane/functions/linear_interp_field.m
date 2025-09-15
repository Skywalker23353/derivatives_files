function data_vec_interp = linear_interp_field(C_field_vec, C_vec, data_vec)
    % Linear interpolation function
    
    data_vec_interp = zeros(length(C_field_vec), 1);
    
    for i = 1:length(C_field_vec)
        C_val = C_field_vec(i);
        
        % Find the index where C_vec >= C_val
        upper_idx = find(C_vec >= C_val, 1);
        
        % Handle edge cases
        if isempty(upper_idx)
            % C_val is larger than all values in C_vec
            data_vec_interp(i) = data_vec(end);
            continue;
        elseif upper_idx == 1
            % C_val is smaller than or equal to first value
            data_vec_interp(i) = data_vec(1);
            continue;
        end
        
        % Linear interpolation between two points
        lower_idx = upper_idx - 1;
        
        % Calculate weights (inverse distance weighting)
        total_distance = C_vec(upper_idx) - C_vec(lower_idx);
        m = (C_val - C_vec(lower_idx)) / total_distance;
        
        % Interpolated value
        data_vec_interp(i) = (1-m) * data_vec(lower_idx) + m * data_vec(upper_idx);
    end
end