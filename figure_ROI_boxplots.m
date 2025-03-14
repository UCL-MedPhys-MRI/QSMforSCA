% figure_ROI_boxplots.m
%
% Script for generating box plots showing regional susceptibility values in
% deep-brain grey matter ROIs, comparing Healthy Controls with SCA
% patients, based on data from the Sickle-UK dataset. Data should be
% generated by SICKLE_EXTRACT_AVERAGES.m
% 
%
%       Copyright (C) University College London, 2025
%
%
% Created by MT Cherukara, January 2024
%
% CHANGELOG:
%
% 2025-03-11 (MTC). Refactored and brought into alignment with the other
%       scripts in this folder, ready for the paper.


clearvars;

% Define UCL colors
colour.g = [143, 153,  62]./255;  % UCL Mid Green 100%
colour.p = [ 98,  32, 113]./255;  % UCL Mid Purple 90%
colour.r = [147,  39,  44]./255;  % UCL Mid Red 100%
colour.b = [  0,  43,  85]./255;  % UCL Mid Blue 100%
colour.r60 = [179, 104, 107]./255;  % UCL Mid Red 60%
colour.b60 = [102, 126, 153]./255;  % UCL Mid Blue 60%

% Load the data
load('SickleUK_QSMData.mat');
load('ROI_names.mat');

% Data Quantities
n_subs = height(tbl_all);
n_rois = length(roi_names);

% Extract a logical vector for healthy controls
is_hc = strcmp(tbl_all.Group,'HC');

% Pre-allocate array for storing data
mat_QSM = zeros(n_subs,n_rois);
mat_R2s = zeros(n_subs,n_rois);

% Loop through ROIs
for rr = 1:n_rois

    % ROI name
    rname = roi_names{rr};

    % Extract susceptibility for this ROI
    vec_QSM = tbl_all.(strcat('QSM_',rname));
    vec_R2s = tbl_all.(strcat('R2s_',rname));

    % Store susceptibility
    mat_QSM(:,rr) = vec_QSM;
    mat_R2s(:,rr) = vec_R2s;

end % for rr = 1:n_rois



%% Box Plot Labels

% Create Label arrays
label_hc = repmat(is_hc,1,n_rois);
label_rois = repmat(1:n_rois,n_subs,1);

% ROI titles
roi_titles = {'L. Caudate';     'R. Caudate';...
              'L. G.P.';        'R. G.P.';...
              'L. Putamen';     'R. Putamen';...
              'L. Thalamus';    'R. Thalamus';...
              'L. Pulvinar';    'R. Pulvinar';...
              'L. Subthal. N.'; 'R. Subthal. N.';...
              'L. S. Nigra';    'R. S. Nigra';...
              'L. Red N.';      'R. Red N.';...
              'L. Dentate';     'R. Dentate'};

%% Box Plot of QSM Values

gsp = 0.01; % spacing variable for the significance stars

% Max and min values for the y-axis
chi_min = min(mat_QSM(:));
chi_max = max(mat_QSM(:));

% Create figure
figure('WindowStyle','normal','Position',[100,200,600+(50*n_rois),600]);

% Beautiful box chart
b1 = boxchart(label_rois(:),mat_QSM(:),'GroupByColor',label_hc(:));
box on; hold on;

% Line at chi=0
plot([0,n_rois+1],[0,0],'k--');

% Labels and axes
axis([0.5,n_rois+0.5,chi_min - (1*gsp),chi_max + (1*gsp)]);
ylabel('Susceptibility (ppm)');
xticks(1:n_rois);
xticklabels(roi_titles);
set(gca,'FontSize',14);
b1(1).BoxFaceColor = colour.g;
b1(2).BoxFaceColor = colour.p;
b1(1).MarkerColor = colour.g;
b1(2).MarkerColor = colour.p;

% Create the legend at the end
legend('Healthy Controls','SCA Patients','','Location','SouthWest');


%% Box Plot of R2-star Values

gsp = 0.1; % spacing variable for the significance stars

% Max and min values for the y-axis
chi_min = min(mat_R2s(:));
chi_max = max(mat_R2s(:));

% Create figure
figure('WindowStyle','normal','Position',[100,300,600+(50*n_rois),600]);

% Beautiful box chart
b2 = boxchart(label_rois(:),mat_R2s(:),'GroupByColor',label_hc(:));
box on; hold on;

% Line at chi=0
plot([0,n_rois+1],[0,0],'k--');

% Labels and axes
axis([0.5,n_rois+0.5,chi_min - (1*gsp),chi_max + (1*gsp)]);
ylabel('R_2* (s^-^1)');
xticks(1:n_rois);
xticklabels(roi_titles);
set(gca,'FontSize',14);
b2(1).BoxFaceColor = colour.g;
b2(2).BoxFaceColor = colour.p;
b2(1).MarkerColor = colour.g;
b2(2).MarkerColor = colour.p;

% Create the legend at the end
legend('Healthy Controls','SCA Patients','','Location','NorthEast');
