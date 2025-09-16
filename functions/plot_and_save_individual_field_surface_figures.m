function plot_and_save_individual_field_surface_figures(fields_struct,D, surface_figures_dir, plot_surface_figures_flag, save_surface_figures_flag)
% 6 figures for field derivatives
    fig_num_base = 3000;
    field_names = fieldnames(fields_struct);
    for i = 1:length(field_names)
        field_name = field_names{i};
        field_data = fields_struct.(field_name);
        field_latex_label = fields_struct.(field_name).latex_label;
        fig_handle = fig_num_base + i*20;
        if plot_surface_figures_flag
            figname = sprintf('$\\frac{\\partial %s}{\\partial c}$',field_latex_label);
            myutils.plot_surf_field(fig_handle, field_data.C_MAT, field_data.Z_MAT/D, field_data.derivative_wrt_C, '$c$', '$z/D$',figname);
        end
        if save_surface_figures_flag
            saveas(gcf, fullfile(surface_figures_dir, sprintf('d%s_dc.png',field_name)));
            saveas(gcf, fullfile(surface_figures_dir, sprintf('d%s_dc.fig',field_name)));
        end
    end
end