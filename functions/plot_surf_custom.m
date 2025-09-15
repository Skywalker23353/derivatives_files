function plot_surf_custom(fig_handle, X, Y, Z, x_label, y_label, z_label, plot_title)
% plot_surf_custom - Plots a surface using surf() with custom labels and title
% Usage:
%   plot_surf_custom(fig_handle, X, Y, Z, x_label, y_label, z_label, plot_title)

figure(fig_handle);
surf(X, Y, Z, 'EdgeColor', 'none');
colorbar;
xlabel(x_label, 'Interpreter', 'latex');
ylabel(y_label, 'Interpreter', 'latex');
zlabel(z_label, 'Interpreter', 'latex');
title(plot_title, 'Interpreter', 'latex');
view(2); % Top-down view for consistency with contourf
set(gcf, 'WindowState', 'maximized');
end
