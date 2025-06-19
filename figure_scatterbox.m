% figure_scatterbox.m
%
% Script for generating the scatter graph (with regression line and
% confidence intervals) for susceptibility (or R2*) values in an ROI, with
% data split into two groups (e.g. SCA / HC) and also a comparative bar
% chart next to each one. Data is taken from the Sickle-UK dataset. Stats
% should be created by SICKLE_EXTRACT_AVERAGES.m
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
% colour.a = [143, 153,  62]./255;  % UCL Mid Green 100%
% colour.b = [ 98,  32, 113]./255;  % UCL Mid Purple 90%
colour.a = [  0,  43,  85]./255;  % UCL Mid Blue 100%
colour.b = [147,  39,  44]./255;  % UCL Mid Red 100%

% Load the data
load('SickleUK_QSMData.mat');
load('ROI_names.mat');

% Data Quantities
n_subs = height(tbl_all);
n_rois = length(roi_names);


%% Choose your ROI and variable to split on (Sex or Group)

% % Choose ROI
% pick_roi = 'Red_nucleus_L';
% ylim_roi = [-0.02,0.12];

% % Choose ROI
% pick_roi = 'Caudate_R';
% ylim_roi = [0.0,0.048];

% Choose ROI
pick_roi = 'Caudate_L';
ylim_roi = [0.005,0.055];

% % Choose ROI
% pick_roi = 'Substantia_nigra_L';
% ylim_roi = [0.0,0.11];

% Define split
% vec_split = strcmp(tbl_all.Group,'SS');
vec_split = strcmp(tbl_all.Sex,'F');

% Choose independent variable
pick_ind = 'Log_Age';


%% Calculate linear fit parameters

% Pull out the data
vec_susc = tbl_all.(strcat('QSM_',pick_roi));
vec_lage = tbl_all.(pick_ind);

% x-values
% x1 = 0; x2 = 14;    % Design fluency
% x1 = 15; x2 = 36;   % Pegboard score
x1 = 2; x2 = 3.6;   % Log(Age)

vec_x = linspace(x1,x2,100); 

% Make vector of x-values
vec_xconf = [vec_x, vec_x(end:-1:1)];

% Clean out NaNs prior to polynomial fitting
vec_fx = vec_lage(~isnan(vec_lage));
vec_fy = vec_susc(~isnan(vec_lage));
vec_split = vec_split(~isnan(vec_lage));

% Split the data
vec_x1 = vec_fx(vec_split == 0);
vec_y1 = vec_fy(vec_split == 0);
vec_x2 = vec_fx(vec_split == 1);
vec_y2 = vec_fy(vec_split == 1);

% Calculate linear fit
[fit_p1,fit_s1] = polyfit(vec_x1,vec_y1,1);
[fit_p2,fit_s2] = polyfit(vec_x2,vec_y2,1);

% Apply both linear fits to our ideal vector x to generate line of best fit
[vec_f1, vec_d1] = polyval(fit_p1,vec_x,fit_s1);
[vec_f2, vec_d2] = polyval(fit_p2,vec_x,fit_s2);

% Reversed y-values
vec_yrev1 = vec_f1(end:-1:1);
vec_yrev2 = vec_f2(end:-1:1);

% y-values of the confidence intervals
vec_yconf1 = [vec_f1 - vec_d1, vec_yrev1 + vec_d1];
vec_yconf2 = [vec_f2 - vec_d2, vec_yrev2 + vec_d2];

%% Plot the scatter graph

% Create figure window
f1 = figure('WindowStyle','normal','Position',[100,200,550,450]);

% Plot confidence interval shapes
p1 = fill(vec_xconf,vec_yconf1,colour.a,'FaceAlpha',0.3);
p1.EdgeColor = 'none';
hold on; box off;
p2 = fill(vec_xconf,vec_yconf2,colour.b,'FaceAlpha',0.3);
p2.EdgeColor = 'none';

% Plot regression lines
plot(vec_x,vec_f1,'Color',colour.a,'LineWidth',2);
plot(vec_x,vec_f2,'Color',colour.b,'LineWidth',2);

% Plot scatter points
scatter(vec_x1, vec_y1, 70, 'o','filled','MarkerFaceColor',colour.a);
scatter(vec_x2, vec_y2, 70, 'd','filled','MarkerFaceColor',colour.b);

% Axes and labels
xlim([x1,x2]);
ylim(ylim_roi);
% xlabel('Pegboard Time (s)');
% xlabel('Design Fluency Score');
xlabel('Log(Age)');
ylabel('Susceptibility (ppm)');
% legend('','','','','Healthy Controls','SCA Patients','Location','SouthEast');
legend('','','','','Male','Female','Location','SouthEast');
legend boxoff;
set(gca,'FontSize',16);


%% Make the Box Plot

% Make a figure
f2 = figure('WindowStyle','normal','Position',[700,200,350,450]);

% Box plot
b1 = boxchart(vec_susc,'GroupByColor',vec_split);

% Colours
b1(1).BoxFaceColor = colour.a;
b1(2).BoxFaceColor = colour.b;
b1(1).MarkerColor = colour.a;
b1(2).MarkerColor = colour.b;

% Axes and labels
ylim(ylim_roi);
xlabel('   ');
% xticks(1);
xticklabels(' ');
ylabel('Susceptibility (ppm)');
set(gca,'FontSize',16);
set(gca,'yaxislocation','right');

