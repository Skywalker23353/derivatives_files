function plot_surface_field_results(field_result, heatRelease, field_name, D, figure_offset, figures_dir, C_MAT, Z_MAT)
    % Plot original field, its derivative, and final result in a single
    % maximized figure
    
    base_figure = 200 + figure_offset * 10;
    
    % Create maximized figure with 3 subplots
    figure(base_figure);
    set(gcf, 'WindowState', 'maximized');  % Maximize the figure
    set(gcf, 'Position', get(0, 'Screensize'));  % Alternative maximization method

    myutils.plot_surf_field(gcf, C_MAT, Z_MAT/D, field_result.field_derivative, ...
            '$c$', '$z/D$', sprintf('$\\frac{\\partial \\langle %s|c\\rangle}{\\partial c}$', field_result.latex_label));
    
    % Adjust subplot spacing
    set(gcf, 'Units', 'normalized');
   
end