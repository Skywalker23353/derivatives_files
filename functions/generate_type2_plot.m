function generate_type2_plot(species_struct, fields_struct, species_name, successful_fields, D, figures_dir)
    % Generate Type 2 plot: All field derivatives for this species (6 subplots)
    % This creates 4 figures total (one per species)

    species_data = species_struct.(species_name);
    field_names = fieldnames(fields_struct);

    % Create figure
    figure_num = 2000 + find(strcmp({'CH4', 'O2', 'CO2', 'H2O'}, species_name));
    figure(figure_num);
    set(gcf, 'WindowState', 'maximized');

    % Calculate subplot arrangement
    n_fields = length(field_names);
    subplot_rows = 2;
    subplot_cols = 3;

    for i = 1:length(field_names)
        field_name = field_names{i};
        field_data = fields_struct.(field_name);
        dfield_name = ['d', species_name, '_d', field_data.short_name];

        subplot(subplot_rows, subplot_cols, i);
        if isfield(species_data.derivative_wrt_fields, dfield_name)
            final_derivative = species_data.derivative_wrt_fields.(dfield_name);
            try
                myutils.plot_contourf(gcf, species_data.C_MAT, species_data.Z_MAT/D, final_derivative, ...
                    '$c$', '$z/D$', sprintf('$\\frac{\\partial %s}{\\partial %s}$', species_data.latex_name, field_data.latex_label));
            catch
                contourf(species_data.C_MAT, species_data.Z_MAT/D, final_derivative, 20);
                colorbar;
                xlabel('$c$', 'Interpreter', 'latex');
                ylabel('$z/D$', 'Interpreter', 'latex');
                title(sprintf('$\\frac{\\partial %s}{\\partial %s}$', species_data.latex_name, field_data.latex_label), 'Interpreter', 'latex');
            end
        end
    end

    % Add overall title
    sgtitle(sprintf('%s Reaction Rate Derivatives w.r.t. All Fields', species_name), 'FontSize', 16, 'FontWeight', 'bold');

    % Save figure
    try
        saveas(gcf, fullfile(figures_dir, sprintf('Type2_%s_all_derivatives.png', species_name)));
        saveas(gcf, fullfile(figures_dir, sprintf('Type2_%s_all_derivatives.fig', species_name)));
        fprintf('   Type 2 figure saved: %s all derivatives\n', species_name);
    catch ME
        fprintf('   Error saving Type 2 figure for %s: %s\n', species_name, ME.message);
    end
end