function save_individual_surface_figures(species_struct, fields_struct, species_name, successful_fields, D, surface_figures_dir, plot_surface_figures_flag, save_surface_figures_flag)
    % 4 figures for species derivatives
    species_data = species_struct.(species_name);
    fig_num_base = 3000;
    fig_handle = fig_num_base;
    if plot_surface_figures_flag
        plot_surf_custom(fig_handle, species_data.C_MAT, species_data.Z_MAT/D, species_data.derivative_wrt_C, ...
            '$c$', '$z/D$', 'Species Derivative', sprintf('%s Species Derivative', species_name));
    end
    if save_surface_figures_flag
        saveas(gcf, fullfile(surface_figures_dir, sprintf('surf_species_derivative_%s.png', species_name)));
        saveas(gcf, fullfile(surface_figures_dir, sprintf('surf_species_derivative_%s.fig', species_name)));
    end

    % 6 figures for field derivatives
    field_names = fieldnames(fields_struct);
    for i = 1:length(field_names)
        field_name = field_names{i};
        field_data = fields_struct.(field_name);
        fig_handle = fig_num_base + i;
        if plot_surface_figures_flag
            plot_surf_custom(fig_handle, field_data.C_MAT, field_data.Z_MAT/D, field_data.derivative_wrt_C, ...
                '$c$', '$z/D$', 'Field Derivative', sprintf('%s Field Derivative (%s)', species_name, field_name));
        end
        if save_surface_figures_flag
            saveas(gcf, fullfile(surface_figures_dir, sprintf('surf_field_derivative_%s_%s.png', species_name, field_name)));
            saveas(gcf, fullfile(surface_figures_dir, sprintf('surf_field_derivative_%s_%s.fig', species_name, field_name)));
        end
    end

    % 24 figures for final derivatives
    for i = 1:length(field_names)
        field_name = field_names{i};
        field_data = fields_struct.(field_name);
        dfield_name = ['d', species_name, '_d', field_data.short_name];
        fig_handle = fig_num_base + 100 + i;
        if isfield(species_data.derivative_wrt_fields, dfield_name)
            final_derivative = species_data.derivative_wrt_fields.(dfield_name);
            if plot_surface_figures_flag
                plot_surf_custom(fig_handle, species_data.C_MAT, species_data.Z_MAT/D, final_derivative, ...
                    '$c$', '$z/D$', 'Final Derivative', sprintf('%s Final Derivative (%s)', species_name, field_name));
            end
            if save_surface_figures_flag
                saveas(gcf, fullfile(surface_figures_dir, sprintf('surf_final_derivative_%s_%s.png', species_name, field_name)));
                saveas(gcf, fullfile(surface_figures_dir, sprintf('surf_final_derivative_%s_%s.fig', species_name, field_name)));
            end
        end
    end
end