function [species_struct, successful_species] = compute_species_numerators(species_configs, field_configs, data_dir, C_MAT, Z_MAT, omega_dot_k_scaling)
    % Compute species numerators (derivatives with respect to c)
    fprintf('\nProcessing %d species (numerators)...\n', size(species_configs, 1));
    species_struct = struct();
    successful_species = {};

    for i = 1:size(species_configs, 1)
        species_name = species_configs{i, 1};
        smooth_suffix = field_configs{i, 2};  % Note: using field_configs smooth suffix
        file_prefix = species_configs{i, 3};
        species_latex = species_configs{i, 4};

        fprintf('\n--- Processing Species %d/%d: %s ---\n', i, size(species_configs, 1), species_name);

        % Load species reaction rate data
        species_filename = sprintf('%s%s', file_prefix, smooth_suffix);
        species_path = sprintf('%s/%s.mat', data_dir, species_filename);

        if ~exist(species_path, 'file')
            fprintf('  Species file not found: %s\n', species_path);
            continue;
        end

        try
            % Load species data
            fprintf('  Loading species reaction rate data...\n');
            species_data = load(species_path);

            % Compute species derivative with respect to c
            fprintf('  Computing derivative d%s/dc...\n', species_name);
            dspecies_dc = compute_dfdr(species_data.DF, C_MAT);
            dspecies_dc = set_boundary_to_zero(dspecies_dc, 'BoundaryWidth', 1, 'Boundaries', {'left','right'});
            dspecies_dc = dspecies_dc./omega_dot_k_scaling;
            % Initialize species_struct with numerator information
            species_struct.(species_name).actual_data = species_data.DF;
            species_struct.(species_name).derivative_wrt_C = dspecies_dc;
            species_struct.(species_name).latex_name = species_latex;
            species_struct.(species_name).C_MAT = C_MAT;
            species_struct.(species_name).Z_MAT = Z_MAT;
            species_struct.(species_name).derivative_wrt_fields = struct();  % Initialize empty

            successful_species{end+1} = species_name;
            fprintf('   ✓ Successfully loaded and computed numerator for %s\n', species_name);

        catch ME
            fprintf('   ✗ Error processing species %s: %s\n', species_name, ME.message);
            continue;
        end
    end
end