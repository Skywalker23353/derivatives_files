clear;clc;
% close all;
addpath('~/MATLAB');
%%
% Reading the data
data_dir = '/work/home/satyam/satyam_files/CH4_jet_PF/2025_Runs/c_cond_stats/C_cond_fields_800';
field = 'Heatrelease_smooth';
heatRelease = load(sprintf('%s/%s.mat',data_dir,field));
field = 'Temperature_smooth';
Temperature = load(sprintf('%s/%s.mat',data_dir,field));
D = 2e-3;
%% Compute dq_dc  

dq_dc = compute_dfdr((heatRelease.DF),heatRelease.C_MAT);
dT_dc = compute_dfdr((Temperature.DF),Temperature.C_MAT);

zeroIdx = find(abs(dT_dc) < 1e-16);
dT_dc(zeroIdx) = 1e-8;

%%
figure(110)
myutils.plot_field(110,heatRelease.C_MAT,heatRelease.Z_MAT/D,heatRelease.DF,'$\langle \dot{\omega}_{T}|c\rangle$');
pbaspect([9 16 1]);
% ylim([0 10]);

figure(120)
myutils.plot_field(120,heatRelease.C_MAT,heatRelease.Z_MAT/D,dq_dc,'$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial c}$')
pbaspect([9 16 1]);
% ylim([0 10]);
%%
figure(30)
myutils.plot_field(30,heatRelease.C_MAT,heatRelease.Z_MAT/D,Temperature.DF,'$\langle T|c\rangle$')
pbaspect([9 16 1]);


figure(40)
myutils.plot_field(40,heatRelease.C_MAT,heatRelease.Z_MAT/D,dT_dc,'$\frac{\partial \langle T|c\rangle}{\partial c}$')
pbaspect([9 16 1]);

%%

dwt_dT = (dq_dc./dT_dc);

figure(80)
myutils.plot_field(80,heatRelease.C_MAT,heatRelease.Z_MAT/D,dwt_dT,'$\frac{\partial \langle \dot{\omega}_{T}|c\rangle}{\partial T}$')
pbaspect([9 16 1]);



