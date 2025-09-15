clear;clc;
% close all;
addpath('~/MATLAB')
addpath('~/satyam_files/CH4_jet_PF/2025_Runs/derivatives_files/functions');
%% Input paths
work_dir = '/work/home/satyam/satyam_files/CH4_jet_PF/2025_Runs/derivatives_files';
data_dir = '/work/home/satyam/satyam_files/CH4_jet_PF/2025_Runs/c_cond_stats/C_cond_fields_800';
%% Load coord data
fprintf('Loading coordinate data from CZ_data.mat...\n');
try
    coord_data = load(fullfile(data_dir, 'CZ_data.mat'));
    C_MAT = coord_data.C_MAT;
    Z_MAT = coord_data.Z_MAT;
    fprintf('? Successfully loaded coordinate data\n');
catch ME
    fprintf('? Error loading coordinate data: %s\n', ME.message);
    fprintf('Please ensure CZ_data.mat exists in: %s\n', data_dir);
    return;
end
%% fieldname
fieldname1 = 'dCH4_dCH4';
dCH4_dCH4 = load(sprintf('%s/species_sensitivities/%s.mat',work_dir,fieldname1));

fieldname2 = 'dCO2_dCH4';
dCO2_dCH4 = load(sprintf('%s/species_sensitivities/%s.mat',work_dir,fieldname2));

fieldname3 = 'dH2O_dCH4';
dH2O_dCH4 = load(sprintf('%s/species_sensitivities/%s.mat',work_dir,fieldname3));

fieldname4 = 'dH2O_dH2O';
dH2O_dH2O = load(sprintf('%s/species_sensitivities/%s.mat',work_dir,fieldname4));

fieldname5 = 'dH2O_dO2';
dH2O_dO2 = load(sprintf('%s/species_sensitivities/%s.mat',work_dir,fieldname5));

fieldname6 = 'dH2O_dT';
dH2O_dT = load(sprintf('%s/species_sensitivities/%s.mat',work_dir,fieldname6));

fieldname7 = 'dO2_dCH4';
dO2_dCH4 = load(sprintf('%s/species_sensitivities/%s.mat',work_dir,fieldname7));

fieldname8 = 'dO2_dCO2';
dO2_dCO2 = load(sprintf('%s/species_sensitivities/%s.mat',work_dir,fieldname8));

fieldname9 = 'dO2_dH2O';
dO2_dH2O= load(sprintf('%s/species_sensitivities/%s.mat',work_dir,fieldname9));

fieldname10 = 'dO2_dO2';
dO2_dO2 = load(sprintf('%s/species_sensitivities/%s.mat',work_dir,fieldname10));

fieldname11 = 'dO2_dT';
dO2_dT = load(sprintf('%s/species_sensitivities/%s.mat',work_dir,fieldname11));

%% Plot
% 1. dCH4_dCH4';
figure()
surf(C_MAT,Z_MAT,dCH4_dCH4.sensitivity);
mx1 = 0.05;min1= 0.5;
smoothing_cycle1 = 2;
smoothed_field = thresholding_and_smoothing(dCH4_dCH4.sensitivity,3,mx1,min1,smoothing_cycle1);
figure()
surf(C_MAT,Z_MAT,smoothed_field);

% 2. dCO2_dCH4';
figure()
surf(C_MAT,Z_MAT,dCO2_dCH4.sensitivity);
mx2 = 0.5;min2 = 0.02;
smoothing_cycle2 = 2;
smoothed_field = thresholding_and_smoothing(dCO2_dCH4.sensitivity,3,mx2,min2,smoothing_cycle2);
figure()
surf(C_MAT,Z_MAT,smoothed_field);

%3. dH2O_dCH4
figure()
surf(C_MAT,Z_MAT,dH2O_dCH4.sensitivity);
mx3 = 0.5;min3= 0.02;
smoothing_cycle3 = 5;
smoothed_field = thresholding_and_smoothing(dH2O_dCH4.sensitivity,3,mx3,min3,smoothing_cycle3);
figure()
surf(C_MAT,Z_MAT,smoothed_field);

%4. dH2O_dH2O
figure()
surf(C_MAT,Z_MAT,dH2O_dH2O.sensitivity);
mx4 = 0.05;min4 = 1;
smoothing_cycle4 = 3;
smoothed_field = thresholding_and_smoothing(dH2O_dH2O.sensitivity,3,mx4,min4,smoothing_cycle4);
figure()
surf(C_MAT,Z_MAT,smoothed_field);

%5. dH2O_dH2O
figure()
surf(C_MAT,Z_MAT,dH2O_dO2.sensitivity);
mx5 = 0.2;min5 = 0.9;
smoothing_cycle5 = 4;
smoothed_field = thresholding_and_smoothing(dH2O_dO2.sensitivity,3,mx5,min5,smoothing_cycle5);
figure()
surf(C_MAT,Z_MAT,smoothed_field);
%%
% %6. dH2O_dT
figure()
surf(C_MAT,Z_MAT,dH2O_dT.sensitivity);
mx6 = 0.1;min6 = 0.9;
smoothing_cycle6 = 4;
smoothed_field = thresholding_and_smoothing(dH2O_dT.sensitivity,3,mx6,min6,smoothing_cycle6);
figure()
surf(C_MAT,Z_MAT,smoothed_field);
%%
% %7. dO2_dCH4
figure()
surf(C_MAT,Z_MAT,dO2_dCH4.sensitivity);
mx7 = 0.1;min7 = 0.9;
smoothing_cycle7 = 3;
smoothed_field = thresholding_and_smoothing(dO2_dCH4.sensitivity,3,mx7,min7,smoothing_cycle7);
figure()
surf(C_MAT,Z_MAT,smoothed_field);
%%
% %8. dO2_dCO2
figure()
surf(C_MAT,Z_MAT,dO2_dCO2.sensitivity);
mx8 = 0.9;min8 = 0.9;
smoothing_cycle8 = 4;
smoothed_field = thresholding_and_smoothing(dO2_dCO2.sensitivity,3,mx8,min8,smoothing_cycle8);
figure()
surf(C_MAT,Z_MAT,smoothed_field);
%%
% %9. dO2_dH2O
figure()
surf(C_MAT,Z_MAT,dO2_dH2O.sensitivity);
mx9 = 0.9;min9 = 0.5;
smoothing_cycle9 = 2;
smoothed_field = thresholding_and_smoothing(dO2_dH2O.sensitivity,3,mx9,min9,smoothing_cycle9);
figure()
surf(C_MAT,Z_MAT,smoothed_field);
%%
% %10. dO2_dO2
figure()
surf(C_MAT,Z_MAT,dO2_dO2.sensitivity);
mx10 = 0.8;min10 = 0.7;
smoothing_cycle10 = 3;
smoothed_field = thresholding_and_smoothing(dO2_dO2.sensitivity,3,mx10,min10,smoothing_cycle10);
figure()
surf(C_MAT,Z_MAT,smoothed_field);
%%
%11. dO2_dT
figure()
surf(C_MAT,Z_MAT,dO2_dT.sensitivity);
mx11 = 0.9;min11 = 0.5;
smoothing_cycle11 = 2;
smoothed_field = thresholding_and_smoothing(dO2_dT.sensitivity,3,mx11,min11,smoothing_cycle11);
figure()
surf(C_MAT,Z_MAT,smoothed_field);
close all;
%%
custom_thld_and_smth_param = cell(11,4);
custom_thld_and_smth_param(:,1) = {fieldname1 fieldname2 fieldname3 fieldname4 fieldname5 fieldname6 fieldname7 fieldname8 fieldname9 fieldname10 fieldname11};
custom_thld_and_smth_param(:,2) = {mx1 mx2 mx3 mx4 mx5 mx6 mx7 mx8 mx9 mx10 mx11 };
custom_thld_and_smth_param(:,3) = {min1 min2 min3 min4 min5 min6 min7 min8 min9 min10 min11 };
custom_thld_and_smth_param(:,4) = {smoothing_cycle1 smoothing_cycle2 smoothing_cycle3 smoothing_cycle4 smoothing_cycle5 smoothing_cycle6 smoothing_cycle7 smoothing_cycle8 smoothing_cycle9 smoothing_cycle10 smoothing_cycle11 };
%%
data = custom_thld_and_smth_param;
save("./species_sensitivities/custom_thld_and_smth_param.mat","data");

