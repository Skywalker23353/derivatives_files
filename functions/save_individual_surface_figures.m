function save_individual_surface_figures(species_struct, fields_struct, species_name, successful_fields, D, surface_figures_dir, plot_surface_figures_flag, save_surface_figures_flag,idx)
    % 4 figures for species derivatives
    field_names = fieldnames(fields_struct);
    species_data = species_struct.(species_name);
    species_latex_name = species_struct.(species_name).latex_name;
    fig_num_base = 3000 + idx;
    fig_handle = fig_num_base;
    if plot_surface_figures_flag
        figname = sprintf('$\\frac{\\partial %s}{\\partial c}$',species_latex_name);
        myutils.plot_surf_field(fig_handle, species_data.C_MAT, species_data.Z_MAT/D, species_data.derivative_wrt_C, '$c$', '$z/D$',figname);
    end
    if save_surface_figures_flag
        saveas(gcf, fullfile(surface_figures_dir, sprintf('d%s_dc.png', species_name)));
        saveas(gcf, fullfile(surface_figures_dir, sprintf('d%s_dc.fig', species_name)));
    end

    % 24 figures for final derivatives
    for i = 1:length(field_names)
        field_name = field_names{i};
        field_latex_label = fields_struct.(field_name).latex_label;
        field_data = fields_struct.(field_name);
        dfield_name = ['d', species_name, '_d', field_data.short_name];
        species_latex_name = species_struct.(species_name).latex_name;
        fig_handle = fig_num_base + 100 + i + 10*idx;
        if isfield(species_data.derivative_wrt_fields, dfield_name)
            final_derivative = species_data.derivative_wrt_fields.(dfield_name);
            if plot_surface_figures_flag
                figname = sprintf('$\\frac{\\partial %s}{\\partial %s}$',species_latex_name,field_latex_label);
                myutils.plot_surf_field(fig_handle, species_data.C_MAT, species_data.Z_MAT/D, final_derivative, '$c$', '$z/D$',figname);
            end
            if save_surface_figures_flag
                saveas(gcf, fullfile(surface_figures_dir, sprintf('d%s_d%s.png', species_name, field_name)));
                saveas(gcf, fullfile(surface_figures_dir, sprintf('d%s_d%s.fig', species_name, field_name)));
            end
        end
    end
end