% figure_scatterbox.m
%
% Script for generating the scatter graph (with regression line) and box
% plots of the QSM (or R2-star) data from the Sickle-UK dataset. Stats
% should be created by SICKLE_EXTRACT_AVERAGES.m and SICKLE_GLM_ANALYSIS.m 
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
load('SickleUK_GLMResults_QSM.mat'); 
load('ROI_names.mat');

% Data Quantities
n_subs = height(tbl_all);
n_rois = length(roi_names);


%% Choose your ROI and variable to split on (Sex or Group)

% Choose ROI
pick_roi = 'Red_nucleus_L';
ylim_roi = [0.005,0.055]; 

% Define split
vec_split = strcmp(tbl_all.Sex,'F');
name_split = 'Sex_F';
leg_labels = {'','','','','Male','Female'};

% % Choose ROI
% pick_roi = 'Substantia_nigra_L';
% ylim_roi = [0.0,0.12];
% 
% % Define split
% vec_split = strcmp(tbl_all.Group,'SS');
% name_split = 'Group_SS';
% leg_labels = {'','','','','Healthy Controls','SCA Patients'};


%% Extract Data

% Pull out the data
vec_susc = tbl_all.(strcat('QSM_',pick_roi));
vec_lage = tbl_all.Log_Age;
% vec_lage = tbl_all.Design_fluency;

% Pull out the current model
pick_mdl = mdl.(pick_roi);

% Find intercept
if pick_mdl.Coefficients{'(Intercept)','pValue'} < 0.05
    a0 = pick_mdl.Coefficients{'(Intercept)','Estimate'};
    a0_p = a0 + pick_mdl.Coefficients{'(Intercept)','SE'};
    a0_m = a0 - pick_mdl.Coefficients{'(Intercept)','SE'};
else
    a0 = 0;
    a0_p = 0;
    a0_m = 0;
end

% Find slope (based on log_age)
a1 = pick_mdl.Coefficients{'Log_Age','Estimate'};
a1_p = a1 + pick_mdl.Coefficients{'Log_Age','SE'};
a1_m = a1 - pick_mdl.Coefficients{'Log_Age','SE'};

% Find extra term to split on
a2 =  pick_mdl.Coefficients{name_split,'Estimate'};
a2_p = a2 + pick_mdl.Coefficients{name_split,'SE'};
a2_m = a2 - pick_mdl.Coefficients{name_split,'SE'};

% x-values (age)
x1 = 2;
x2 = 3.6;
vec_x = linspace(x1,x2,100);
vec_xconf = [vec_x, vec_x(end:-1:1)];

% Construct best-fit line (for the ones without the split)
vec_y1 = a0 + (a1.*vec_x) - a2;
vec_y1_lo = a0 + (a1_m.*vec_x) - a2;
vec_y1_hi = a0 + (a1_p.*vec_x) - a2;
vec_y1conf = [vec_y1_lo, vec_y1_hi(end:-1:1)];

% Construct best-fit line (for the ones with the split)
vec_y2 = a0 + (a1.*vec_x);
vec_y2_lo = a0 + (a1_m.*vec_x);
vec_y2_hi = a0 + (a1_p.*vec_x);
vec_y2conf = [vec_y2_lo, vec_y2_hi(end:-1:1)];


%% Make the Scatter Plot

% Make a figure
figure('WindowStyle','normal','Position',[100,200,550,450]);

% Plot the error-bound line
p1 = fill(vec_xconf,vec_y1conf,colour.a,'FaceAlpha',0.3);
p1.EdgeColor = 'none';
hold on; box off;
p2 = fill(vec_xconf,vec_y2conf,colour.b,'FaceAlpha',0.3);
p2.EdgeColor = 'none';

% Plot scatter points
scatter(vec_lage(~vec_split), vec_susc(~vec_split), 70, 'o','filled','MarkerFaceColor',colour.a);
scatter(vec_lage(vec_split),  vec_susc(vec_split),  70, 'd','filled','MarkerFaceColor',colour.b);

% Plot regression lines
plot(vec_x,vec_y1,'-','LineWidth',2,'Color',colour.a);
plot(vec_x,vec_y2,'-','LineWidth',2,'Color',colour.b);

% Axes and labels
xlim([2,3.6]);
ylim(ylim_roi);
xlabel('Log(Age)');
ylabel('Susceptibility (ppm)');
legend(leg_labels,'Location','SouthEast');
set(gca,'FontSize',18);

% % Background Colour
% set(gcf,'color',[201, 147, 150]./255);
% set(gca,'color',[201, 147, 150]./255);


% Make the Box Plot

% Make a figure
figure('WindowStyle','normal','Position',[700,200,350,450]);

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
ylabel('Susceptibility (ppm)');
set(gca,'FontSize',18);
set(gca,'yaxislocation','right');

% % Background Colour
% set(gcf,'color',[201, 147, 150]./255);
% set(gca,'color',[201, 147, 150]./255);

