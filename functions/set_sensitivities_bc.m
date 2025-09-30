function dspecies_dc_modified = set_sensitivities_bc(dspecies_dc, Z_MAT, z_ref)
    % Apply downstream replacement for species derivatives
    % Replaces all values downstream of z_ref with the value at z_ref
    %
    % INPUTS:
    %   dspecies_dc - 2D matrix of species derivatives with respect to c
    %                 (rows = Z locations, columns = radial positions)
    %   Z_MAT       - 2D matrix of Z coordinates (axial direction)
    %   z_ref       - Reference Z value for downstream replacement
    %
    % OUTPUT:
    %   dspecies_dc_modified - Modified derivative matrix with downstream replacement
    
    fprintf('  Applying downstream replacement at z_ref = %.4f...\n', z_ref);
    
    % Initialize output matrix
    dspecies_dc_modified = dspecies_dc;
    
    % Since rows correspond to Z locations, extract Z values from first column of Z_MAT
    z_values = Z_MAT(:, 1);
    
    % Find the index closest to z_ref
    [~, z_ref_idx] = min(abs(z_values - z_ref));
    replacement_val = dspecies_dc(z_ref_idx,end);
    
    dspecies_dc_modified(:,end) = replacement_val;
    
    fprintf('   ✓ C = 1 boundary replacement completed\n');
end