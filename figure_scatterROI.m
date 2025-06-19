% figure_scatterROI.m
%
% Script for generating the scatter graph (with regression line and
% confidence intervals) for susceptibility (or R2*) values in an ROI, with
% data from the Sickle-UK dataset. Stats should be created by
% SICKLE_EXTRACT_AVERAGES.m
% 
%
%       Copyright (C) University College London, 2025
%
%
% Created by MT Cherukara, March 2024
%
% CHANGELOG:
%
% 2025-03-11 (MTC). Refactored and brought into alignment with the other
%       scripts in this folder, ready for the paper.
%
% 2025-06-19 (MTC). Forked this script into two, one for the
%       single-variable scatter plot (this script) and a separate one
%       (FIGURE_SCATTERBOX.m) for split-variable scatter graph and box
%       plot.


clearvars;
close all;

%% Load the Data

% Define UCL colors
colour.a = [143, 153,  62]./255;  % UCL Mid Green 100%
colour.b = [ 98,  32, 113]./255;  % UCL Mid Purple 90%
% colour.a = [  0,  43,  85]./255;  % UCL Mid Blue 100%
% colour.b = [147,  39,  44]./255;  % UCL Mid Red 100%

% Load the data
load('SickleUK_QSMData.mat');
load('ROI_names.mat');

% Data Quantities
n_subs = height(tbl_all);
n_rois = length(roi_names);


%% Choose your ROI and variable to split on (Sex or Group)

% Choose ROI
pick_roi = 'Red_nucleus_L';
ylim_roi = [-0.02,0.12];

% % Choose ROI
% pick_roi = 'Caudate_R';
% ylim_roi = [0.0,0.048];

% Define split
vec_split = strcmp(tbl_all.Group,'SS');

% Choose independent variable
pick_ind = 'Design_fluency';


%% Calculate linear fit parameters

% Pull out the data
vec_susc = tbl_all.(strcat('QSM_',pick_roi));
vec_lage = tbl_all.(pick_ind);


% x-values
x1 = 0; x2 = 14;    % Design fluency
% x1 = 15; x2 = 36;   % Pegboard score
% x1 = 2; x2 = 3.6;   % Log(Age)

vec_x = linspace(x1,x2,100); 

% Make vector of x-values
vec_xconf = [vec_x, vec_x(end:-1:1)];

% Clean out NaNs prior to polynomial fitting
vec_fx = vec_lage(~isnan(vec_lage));
vec_fy = vec_susc(~isnan(vec_lage));

% Figure out line of best fit
[fit_p,fit_s] = polyfit(vec_fx,vec_fy,1);

% Apply the polyfit to vec_x
[vec_y, vec_d] = polyval(fit_p,vec_x,fit_s);

% Reversed y-values
vec_yrev = vec_y(end:-1:1);

% y-values of the confidence intervals
vec_yconf = [vec_y - vec_d, vec_yrev + vec_d];

%% Plot the scatter graph

% Create figure window
f1 = figure('WindowStyle','normal','Position',[100,200,550,450]);

% Plot confidence interval shape
p1 = fill(vec_xconf,vec_yconf,[0.5,0.5,0.5,],'FaceAlpha',0.3);
p1.EdgeColor = 'none';
hold on; box off;

% Plot regression line
plot(vec_x,vec_y,'k-','LineWidth',2);

% Plot scatter points
scatter(vec_lage(~vec_split), vec_susc(~vec_split), 70, 'o','filled','MarkerFaceColor',colour.a);
scatter(vec_lage(vec_split),  vec_susc(vec_split),  70, 'd','filled','MarkerFaceColor',colour.b);

% Axes and labels
xlim([x1,x2]);
ylim(ylim_roi);
xlabel('Pegboard Time (s)');
% xlabel('Design Fluency Score');
ylabel('Susceptibility (ppm)');
legend('','','Healthy Controls','SCA Patients','Location','SouthWest');
legend boxoff;
set(gca,'FontSize',16);

