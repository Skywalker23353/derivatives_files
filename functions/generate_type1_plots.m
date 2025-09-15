function generate_type1_plots(species_struct, fields_struct, species_name, successful_fields, D, figures_dir)
    % Generate Type 1 plots: Individual derivative analysis (3 subplots each)
    % This creates 24 figures total (6 fields × 4 species)

    species_data = species_struct.(species_name);
    field_names = fieldnames(fields_struct);

    for i = 1:length(field_names)
        field_name = field_names{i};
        field_data = fields_struct.(field_name);
        dfield_name = ['d', species_name, '_d', field_data.short_name];

        % Create figure
        figure_num = 1000 + i + length(successful_fields) * find(strcmp({'CH4', 'O2', 'CO2', 'H2O'}, species_name));
        figure(figure_num);
        set(gcf, 'WindowState', 'maximized');

        % Subplot 1: Species derivative w.r.t. c
        subplot(1, 3, 1);
        try
            myutils.plot_contourf(gcf, species_data.C_MAT, species_data.Z_MAT/D, species_data.derivative_wrt_C, ...
                '$c$', '$z/D$', sprintf('$\\frac{\\partial \\langle %s|c\\rangle}{\\partial c}$', species_data.latex_name));
        catch
            contourf(species_data.C_MAT, species_data.Z_MAT/D, species_data.derivative_wrt_C, 20);
            colorbar;
            xlabel('$c$', 'Interpreter', 'latex');
            ylabel('$z/D$', 'Interpreter', 'latex');
            title(sprintf('$\\frac{\\partial \\langle %s|c\\rangle}{\\partial c}$', species_data.latex_name), 'Interpreter', 'latex');
        end

        % Subplot 2: Field derivative w.r.t. c
        subplot(1, 3, 2);
        try
            myutils.plot_contourf(gcf, field_data.C_MAT, field_data.Z_MAT/D, field_data.derivative_wrt_C, ...
                '$c$', '$z/D$', sprintf('$\\frac{\\partial \\langle %s|c\\rangle}{\\partial c}$', field_data.latex_label));
        catch
            contourf(field_data.C_MAT, field_data.Z_MAT/D, field_data.derivative_wrt_C, 20);
            colorbar;
            xlabel('$c$', 'Interpreter', 'latex');
            ylabel('$z/D$', 'Interpreter', 'latex');
            title(sprintf('$\\frac{\\partial \\langle %s|c\\rangle}{\\partial c}$', field_data.latex_label), 'Interpreter', 'latex');
        end

        % Subplot 3: Final derivative (sensitivity)
        subplot(1, 3, 3);
        if isfield(species_data.derivative_wrt_fields, dfield_name)
            final_derivative = species_data.derivative_wrt_fields.(dfield_name);
            try
                myutils.plot_contourf(gcf, species_data.C_MAT, species_data.Z_MAT/D, final_derivative, ...
                    '$c$', '$z/D$', sprintf('$\\frac{\\partial \\langle %s|c\\rangle}{\\partial %s}$', species_data.latex_name, field_data.latex_label));
            catch
                contourf(species_data.C_MAT, species_data.Z_MAT/D, final_derivative, 20);
                colorbar;
                xlabel('$c$', 'Interpreter', 'latex');
                ylabel('$z/D$', 'Interpreter', 'latex');
                title(sprintf('$\\frac{\\partial \\langle %s|c\\rangle}{\\partial %s}$', species_data.latex_name, field_data.latex_label), 'Interpreter', 'latex');
            end
        end

        % Add overall title
        sgtitle(sprintf('%s Derivative Analysis w.r.t. %s', species_name, field_name), 'FontSize', 16, 'FontWeight', 'bold');

        % Save figure
        try
            saveas(gcf, fullfile(figures_dir, sprintf('Type1_%s_vs_%s_analysis.png', species_name, field_name)));
            saveas(gcf, fullfile(figures_dir, sprintf('Type1_%s_vs_%s_analysis.fig', species_name, field_name)));
            fprintf('   Type 1 figure saved: %s vs %s\n', species_name, field_name);
        catch ME
            fprintf('   Error saving Type 1 figure for %s vs %s: %s\n', species_name, field_name, ME.message);
        end
    end
end