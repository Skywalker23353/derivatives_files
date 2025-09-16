function save_species_results(species_struct, fields_struct, sensitivities_dir)
    % Save results for each species-field combination
    species_names = fieldnames(species_struct);
    for i = 1:length(species_names)
        species_name = species_names{i};
        species_data = species_struct.(species_name);
        field_names = fieldnames(fields_struct);
        for j = 1:length(field_names)
            field_name = field_names{j};
            field_data = fields_struct.(field_name);
            dfield_name = ['d', species_name, '_d', field_data.short_name];
            if isfield(species_data.derivative_wrt_fields, dfield_name)
                sensitivity_data = struct();
                sensitivity_data.sensitivity = species_data.derivative_wrt_fields.(dfield_name);
                sensitivity_data.species_name = species_name;
                sensitivity_data.field_name = field_name;
                sensitivity_data.species_latex = species_data.latex_name;
                sensitivity_data.field_latex = field_data.latex_label;
                output_filename = sprintf('d%s_d%s.mat', species_name, field_data.short_name);
                output_path = fullfile(sensitivities_dir, output_filename);
                try
                    save(output_path, '-struct', 'sensitivity_data');
                    fprintf('   Saved: %s\n', output_filename);
                catch ME
                    fprintf('   Error saving %s: %s\n', output_filename, ME.message);
                end
            end
        end
    end
end