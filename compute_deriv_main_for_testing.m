clear;clc;
%%
case_dir = '/work/home/satyam/satyam_files/Burner_6x_LSA';
grid_mesh_path = case_dir;
load(sprintf('%s/structured_grid_generate_and_interp_scripts/structured_grid_from_LES_grid_with_zero.mat',grid_mesh_path));
%%

[~,col] = size(R1);
mid = (col+1)/2;
r = R1(:,:);
r_h = R1(:,mid:end);
z_h = Z1(:,mid:end);
f = r_h.^3.*z_h;
f_flip = flip(f,2);
fn = [f_flip(:,1:end-1),f];
deriv = 3*(r_h.^2).*z_h;
deriv_flip = flip(deriv,2);
actual_deriv = [deriv_flip(:,1:end-1),deriv];
%test_deriv = compute_dfdr(fn,R1);
test_deriv = compute_dfdr(fn,r);
test_deriv(:,1:mid) = -test_deriv(:,1:mid); 
fname = 'test_deriv_files';
mkdir(sprintf("./%s",fname));
%save(sprintf("./%s/testfn2.mat",fname),'fn','actual_deriv','test_deriv');
%%

% [~,col] = size(R1);
% mid = (col+1)/2;
% r = R1(1,:);
% r_h = R1(1,mid:end);
% z_h = Z1(1,mid:end);
% f = r_h.^3;
% f_flip = flip(f);
% fn = [f_flip(1,1:end-1),f];
% deriv = 3*(r_h.^2);
% deriv_flip = flip(deriv);
% actual_deriv = [deriv_flip(1,1:end-1),deriv];
% %test_deriv = compute_dfdr(fn,R1);
% test_deriv = compute_deriv(fn,r);
% test_deriv(1,1:mid) = -test_deriv(1,1:mid); 
% fname = 'test_deriv_files';
% mkdir(sprintf("./%s",fname));
% save(sprintf("./%s/testfn1.mat",fname),'fn','actual_deriv','test_deriv');
% %%
% fig1 = figure(1);clf;
% set(fig1,'Color',[1 1 1]);
% ax1 = axes('Position',[0.1 0.15 0.35 0.75]);
% ph1 = plot(r,fn,'s-b',"DisplayName",'$fn = r^3$',"MarkerSize",20);
% lh1 = legend(ax1,"Box","off","Location","north");
% title(ax1,"$Function$");
% xlabel(ax1,"$r$");
% ylabel(ax1,"$z$");
% ax2 = axes('Position',[0.55 0.15 0.35 0.75]);
% ph2 = plot(r,actual_deriv,'s-k',"DisplayName","$ Actual$","MarkerSize",20);hold on;
% ph3 = plot(r,test_deriv,'-*r',"DisplayName","$ Test$","MarkerSize",20);
% lh2 = legend("FontSize",30,"Box","off","Location","north");
% title(ax2,"$Derivative$");
% xlabel(ax2,"$r$");
% ylabel(ax2,"$z$");
% hold off;
% axis(ax1,[-0.001 0.001 0 1.2e-9]);
% axis(ax2,[-0.001 0.001 0 3e-6]);
% xticks(ax1,[-0.001 0 0.001]);
% xticklabels(ax1,["0.001","0","0.001"]);
% xticks(ax2,[-0.001 0 0.001]);
% xticklabels(ax2,["0.001","0","0.001"]);
% set(findall(fig1,'-property','Fontsize'),'Fontsize',27);
% set(findall(fig1,'-property','Interpreter'),'Interpreter','latex');
% set(findall(fig1,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
%%
%load("./test_deriv_files/testfn1.mat");
plot_field(1,R1,Z1,fn,0,0);
title(1,"Function");
plot_field(2,R1,Z1,actual_deriv,0,0);

plot_field(3,R1,Z1,test_deriv,0,0);

error = test_deriv - actual_deriv;

plot_field(4,R1,Z1,error,0,0);

%%
close all;