clear;clc;
% close all;
addpath('~/MATLAB');
%%
% Reading the data
data_dir = '/work/home/satyam/satyam_files/CH4_jet_PF/2025_Runs/c_cond_stats/C_cond_fields_800_old';
field = 'Heatrelease_smooth';
heatRelease = load(sprintf('%s/%s.mat',data_dir,field));
field = 'Temperature_smooth';
Temperature = load(sprintf('%s/%s.mat',data_dir,field));
D = 2e-3;
z_axis = heatRelease.Z_MAT(:,1);
z_idx = find((heatRelease.Z_MAT)/D >= 8.5,1);
%% Test function
% fn = Temperature.C_MAT.^2;
% fn1 = sin(Temperature.C_MAT*pi*2);
% dfn_dc = compute_dfdr(fn,Temperature.C_MAT);
% dfn1_dc = compute_dfdr(fn1,Temperature.C_MAT);
%% Compute dq_dc  

dq_dc = compute_dfdr((heatRelease.DF(1:z_idx,:)),heatRelease.C_MAT(1:z_idx,:));
dT_dc = compute_dfdr((Temperature.DF(1:z_idx,:)),Temperature.C_MAT(1:z_idx,:));

zeroIdx = find(abs(dT_dc) < 1e-16);
dT_dc(zeroIdx) = 1e-8;
%% Test fns
% figure(100)
% z_idx = find((Temperature.Z_MAT)/D >= 1,1);
% plot(Temperature.C_MAT(z_idx,:),fn(z_idx,:),'-kx','DisplayName','$f = c^2$','MarkerSize',12);hold on ;
% plot(Temperature.C_MAT(z_idx,:),dfn_dc(z_idx,:),'-bs','DisplayName','$\frac{\partial fn}{\partial c}$','MarkerSize',12);
% plot(Temperature.C_MAT(z_idx,:),fn1(z_idx,:),'-rx','DisplayName','$f = sin(2\pi c)$','MarkerSize',12);
% plot(Temperature.C_MAT(z_idx,:),dfn1_dc(z_idx,:),'-cs','DisplayName','$\frac{\partial fn}{\partial c}$','MarkerSize',12);
% 
% xlabel('$c$','Interpreter','latex');
% legend('Interpreter','latex','FontSize',33);
% set(findall(gcf,'-property','Fontsize'),'Fontsize',33);
% title('$\frac{z}{D} = 1$','Interpreter','latex','FontSize',33);
% hold off;
%%
figure(11)
myutils.plot_field(11,heatRelease.C_MAT(1:z_idx,:),heatRelease.Z_MAT(1:z_idx,:)/D,heatRelease.DF(1:z_idx,:),'$\langle \dot{\omega}_{T}|c\rangle$');
pbaspect([9 16 1]);
% ylim([0 10]);

figure(12)
myutils.plot_field(12,heatRelease.C_MAT(1:z_idx,:),heatRelease.Z_MAT(1:z_idx,:)/D,dq_dc(1:z_idx,:),'$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial c}$')
pbaspect([9 16 1]);
% ylim([0 10]);
%%
figure(3)
myutils.plot_field(3,heatRelease.C_MAT(1:z_idx,:),heatRelease.Z_MAT(1:z_idx,:)/D,Temperature.DF(1:z_idx,:),'$\langle T|c\rangle$')
pbaspect([9 16 1]);


figure(4)
myutils.plot_field(4,heatRelease.C_MAT(1:z_idx,:),heatRelease.Z_MAT(1:z_idx,:)/D,dT_dc(1:z_idx,:),'$\frac{\partial \langle T|c\rangle}{\partial c}$')
pbaspect([9 16 1]);

%%
% figure(5)
% z_idx = find((Temperature.Z_MAT)/D >= 0,1);
% plot(Temperature.C_MAT(z_idx,:),Temperature.DF(z_idx,:),'-kx','DisplayName','$\langle T|c\rangle$','MarkerSize',12);hold on ;
% plot(Temperature.C_MAT(z_idx,:),dT_dc(z_idx,:),'-bs','DisplayName','$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial c}$','MarkerSize',12);
% xlabel('$c$','Interpreter','latex');
% legend('Interpreter','latex','FontSize',33);
% set(findall(gcf,'-property','Fontsize'),'Fontsize',33);
% title('$\frac{z}{D} = 1$','Interpreter','latex','FontSize',33);
% hold off;
%%

% dwt_dT_log = log(abs(dq_dc(1:z_idx,:)./dT_dc(1:z_idx,:)) + 1e-8);
dwt_dT = (dq_dc(1:z_idx,:)./dT_dc(1:z_idx,:));

% figure(6)
% myutils.plot_field(6,Temperature.C_MAT,Temperature.Z_MAT/D,dwt_dT,'$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial T}$')
% pbaspect([9 16 1]);
% % xlim([0.2 1]);
% ylim([0 10]);

%%
% window = 3;
% smooth_dwT_dT1 = myutils.f_return_smooth_field(dwt_dT,window,'row');
% smooth_dwT_dT = myutils.f_return_smooth_field(smooth_dwT_dT1,window,'col');
% figure()
% surf(heatRelease.C_MAT(1:z_idx,:),heatRelease.Z_MAT(1:z_idx,:)/D,dq_dc);
% figure()
% surf(heatRelease.C_MAT(1:z_idx,:),heatRelease.Z_MAT(1:z_idx,:)/D,dT_dc);
% figure()
% surf(heatRelease.C_MAT(1:z_idx,:),heatRelease.Z_MAT(1:z_idx,:)/D,dwt_dT);
% figure()
% surf(Temperature.C_MAT,Temperature.Z_MAT/D,smooth_dwT_dT);
%%
figure(8)
myutils.plot_field(8,heatRelease.C_MAT(1:z_idx,:),heatRelease.Z_MAT(1:z_idx,:)/D,dwt_dT(1:z_idx,:),'$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial T}$')
pbaspect([9 16 1]);
% xlim([0.2 1]);
% ylim([0 10]);


