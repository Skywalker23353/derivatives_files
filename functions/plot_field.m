function [h] = plot_field(indx,X_MAT,Y_MAT,F_MAT,C_MAT,c_indx)
    
    f_size = 24;
    figure(indx);
    [pc,h] = contourf(X_MAT,Y_MAT,F_MAT,200);
    set(h,'LineColor','None');
    shading flat
    colormap jet
    colorbar;
%   caxis([0,0.3]);
    if c_indx == 1
        hold on;
        contour(X_MAT,Y_MAT,C_MAT,[0,0],'k');
    end
    xlabel('$x/D$','Interpreter','latex','FontSize',f_size);
    ylabel('$\frac{y}{D}$','Interpreter','latex','Rotation',0,'FontSize',f_size);
    set(gca,'FontSize',f_size);

end

