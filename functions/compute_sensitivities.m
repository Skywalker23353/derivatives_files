function species_struct = compute_sensitivities(species_struct, fields_struct, successful_species, successful_fields, threshold_and_smooth_results_flag, sensitivities_dir, smth_window,z_ref)
    % Compute final derivatives/sensitivities: d(species_reaction_rate)/d(field)
    fprintf('\nComputing sensitivities (final derivatives)...\n');
    
    for i = 1:length(successful_species)
        species_name = successful_species{i};
        species_info = species_struct.(species_name);
        Z_MAT = species_info.Z_MAT;

        fprintf('\n--- Computing sensitivities for %s ---\n', species_name);

        % Get species derivative (numerator)
        dspecies_dc = species_info.derivative_wrt_C;

        % Process each field for this species
        successful_combinations = {};

        for j = 1:length(successful_fields)
            field_name = successful_fields{j};
            field_info = fields_struct.(field_name);

            fprintf('  Computing d%s/d%s...\n', species_name, field_name);

            try
                % Get field derivative (denominator)
                df_dc = field_info.derivative_wrt_C;

                % Handle zero values in denominator
                zero_idx = find(abs(df_dc) < 1e-16);
                if ~isempty(zero_idx)
                    fprintf('    Handling %d zero values in denominator...\n', length(zero_idx));
                    df_dc(zero_idx) = 1e-8;
                end

                % Compute final derivative: d(species_reaction_rate)/d(field)
                final_derivative = dspecies_dc ./ df_dc;

                % Apply thresholding and smoothing if requested
                if threshold_and_smooth_results_flag
                    final_derivative = custom_thresholding_and_smoothing(final_derivative, sensitivities_dir, species_name, field_info.short_name, smth_window);
                end
                final_derivative = set_sensitivities_bc(final_derivative, Z_MAT, z_ref);
%                 final_derivative =  apply_smoothing_ignore_boundaries_1(final_derivative,3,{'all'},1);
                % Store results in species_struct
                dfield_name = ['d', species_name, '_d', field_info.short_name];
                species_struct.(species_name).derivative_wrt_fields.(dfield_name) = final_derivative;

                successful_combinations{end+1} = sprintf('%s_%s', species_name, field_name);
                fprintf('     ✓ Successfully computed d%s/d%s\n', species_name, field_name);

            catch ME
                fprintf('     ✗ Error computing d%s/d%s: %s\n', species_name, field_name, ME.message);
                continue;
            end
        end

        fprintf('   Successfully computed sensitivities for %s with %d field combinations\n', species_name, length(successful_combinations));
    end
end