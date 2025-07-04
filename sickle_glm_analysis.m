% sickle_glm_analysis.m
%
% Group level analysis of QSM data from the Sickle-UK cohort. Data should
% be processed using SICKLE_QSM_PIPELINE.m and SICKLE_EXTRACT_AVERAGES.m
% and stored as a MATLAB Table in a .mat file, which can be read in to this
% script. 
%
% We are trying to fit a general linear model that takes ROI mean
% susceptibility as the response variable, and predictors 'age', 'group'
% (i.e. SCA status) and the subject's average score from the 'pegboard'
% task and design fluency tests.
% 
%
%       Copyright (C) University College London, 2025
%
%
% Created by MT Cherukara, October 2023
%
% CHANGELOG:
%
% 2025-03-11 (MTC). Refactored and brought into alignment with the other
%       scripts in this folder, ready for the paper.


clearvars

%% Set-up Script Options

% ROIs
load('ROI_names.mat');
n_rois = length(roi_names);

% Load in data table
load('SickleUK_QSMData.mat');
n_subs = height(tbl_all);

% Choose Modality name 'QSM' or 'R2s'
mname = 'QSM';

% Pre-allocate ROI name arrays
roi_names_qsm = cell(size(roi_names));
roi_names_r2s = cell(size(roi_names));

% Generate table heading names
for rr = 1:n_rois

    roi_names_qsm{rr} = strcat('QSM_',roi_names{rr});
    roi_names_r2s{rr} = strcat('R2s_',roi_names{rr});

end % for rr = 1:n_rois


%% Loop Through ROIs and Fit the GLM

% Pre-allocate structure for storing the results
mdl = struct();

% Create an empty array for holding the following p-values of each fitting
% parameter (currently 6) 
arr_results = ones(6,n_rois);


for rr = 1:n_rois

    rname = roi_names{rr};

    % Specify model
    modelspec = [strcat(mname,'_',rname), ' ~ Sex + Group + Log_Age + Pegboard_R + Design_fluency '];

    % Fit GLM
    mdl.(rname) = fitglm(tbl_all,modelspec);

    % Store p-values
    arr_results(:,rr) = mdl.(rname).Coefficients.pValue;

end % rr = 1:n_rois

% Convert results into a table
tbl_results = array2table(arr_results,'VariableNames',roi_names);

% Add row-names to the table
tbl_results.Properties.RowNames = mdl.(roi_names{1}).CoefficientNames;



%% Loop through Results and print

res_RS = zeros(n_rois,1);
res_pv = zeros(n_rois,1);

clc;

for rr = 1:n_rois

    rname = roi_names{rr};

    % Identify which coefficients are significant
    sigCoefs = mdl.(rname).Coefficients.pValue < 0.05;

    % Store adjusted Rsquared values and p-value
    res_RS(rr) = mdl.(rname).Rsquared.Adjusted;
    res_pv(rr) = coefTest(mdl.(rname));

    % if res_pv(rr) < 0.05
    
        % Print Model statistics
        fprintf('\n%20s ',rname);
        fprintf('adjusted R^2 = %f\t',res_RS(rr));
        fprintf('p-value = %.2e \n',res_pv(rr));
    
        if any(sigCoefs(2:end))
            sigCoefs(1) = false;
            nameCoefs = mdl.(rname).CoefficientNames(sigCoefs)';
            stats = [nameCoefs, num2cell(mdl.(rname).Coefficients.pValue(sigCoefs))]';
            fprintf('%20s has significant covariance in',rname);
            fprintf('\n%24s (%8.3g)',stats{:});
            fprintf('.\n');
        end
    % end

end % rr = 1:n_rois


%% Save the Results
save(strcat('SickleUK_GLMResults_',mname,'.mat'),'tbl_results','res_RS','res_pv','mdl');


%% Leave-One-Out Analysis for Variance Explained

% Variables
var_names = {'Log_Age'; 'Sex'; 'Group'; 'Pegboard_R'; 'Design_fluency' };
n_var = length(var_names);

% Pre-allocate table to store variance explained results
tbl_varexp = table(roi_names,res_pv,res_RS);

% Loop over variables to leave out
for vv = 1:n_var

    % Pull out variable name
    vname = var_names{vv};

    % Make a list of all the others
    list_var = var_names;
    idx = strcmp(var_names,vname);
    list_var(idx) = [];

    % Construct model spec
    modelspec_1 = [list_var{1},' + ',list_var{2},' + ',list_var{3},' + ',list_var{4}];

    % Pre-allocate vector for R-squared
    vec_RS = zeros(n_rois,1);

    % Loop over ROIs and calculate the GLM
    for rr = 1:n_rois

        rname = roi_names{rr};

        % Construct model spec
        modelspec = [mname,'_',rname,' ~ ',modelspec_1];

        % Fit GLM
        varmdl = fitglm(tbl_all,modelspec);

        % Store R-squraed vector
        vec_RS(rr) = varmdl.Rsquared.Adjusted;

    end % for rr = 1:n_rois 

    % Set negative R^2 values to 0
    vec_RS(vec_RS<0) = 0;

    % Calculate variance explained by this variable
    vec_VE = res_RS - vec_RS;
    vec_VE(vec_VE<0) = 0;

    % Put it in the table
    tbl_varexp = addvars(tbl_varexp,vec_RS,vec_VE,...
                         'NewVariableNames',{strcat('RS_',vname),strcat('VE_',vname)});

end % for vv = 1:n_var 

% Save the variance explained data
save('SickleUK_GLMResults_VarExp.mat','tbl_varexp');



%% Calculate Regional Averages (for Table)

% Extract a logical vector for healthy controls
is_hc = strcmp(tbl_all.Group,'HC');

vec_age_hc = tbl_all.Age(is_hc);
vec_age_ss = tbl_all.Age(~is_hc);

av_susc_hc = zeros(n_rois,2);
av_susc_ss = zeros(n_rois,2);
mat_susc = zeros(n_subs,n_rois);
vec_pvals = zeros(n_rois,1);

% Loop through ROIs
for rr = 1:n_rois

    % ROI name
    rname = roi_names{rr};

    % Extract susceptibility for this ROI
    vec_susc = tbl_all.(strcat(mname,'_',rname));

    % Separate HC and SCD values
    vec_hc = vec_susc(is_hc);
    vec_ss = vec_susc(~is_hc);

    % Calculate averages
    av_susc_hc(rr,1) = mean(vec_hc,'omitnan');
    av_susc_hc(rr,2) = std(vec_hc,'omitnan');
    av_susc_ss(rr,1) = mean(vec_ss,'omitnan');
    av_susc_ss(rr,2) = std(vec_ss,'omitnan');

    % Store susceptibility
    mat_susc(:,rr) = vec_susc;

    % Do a t-test
    [~,vec_pvals(rr)] = ttest2(vec_hc,vec_ss);

end % for rr = 1:n_rois

tbl_av = table(roi_names,av_susc_hc,av_susc_ss,vec_pvals);


%% Zero-Order Correlations of Cognitive Data

% Extract a logical vector for healthy controls
is_hc = strcmp(tbl_all.Group,'HC');

% Specify which variables we want to examine
var_names = {'Age';'Pegboard_L';'Pegboard_R';'Design_fluency' };

for vv = 1:length(var_names)
    
    % Pull out the data
    vec_cogn = tbl_all.(var_names{vv});

    % Do the t-test
    [h1, p1] = ttest2(vec_cogn(is_hc),vec_cogn(~is_hc));

    % Print the result
    fprintf('%s p-value: %.4f \n',var_names{vv},p1);


end


%% Zero-Order Correlations of Susceptibility Data

% Extract a logical vector for healthy controls
is_hc = strcmp(tbl_all.Group,'HC');
% is_hc = strcmp(tbl_all.Sex,'F');

for vv = 1:length(roi_names)
    
    % Pull out the data
    vec_cogn = tbl_all.(strcat(mname,'_',roi_names{vv}));

    % Do the t-test
    [h1, p1] = ttest2(vec_cogn(is_hc),vec_cogn(~is_hc));

    % Print the result
    if p1 < 0.05
        fprintf('%s p-value: %.4f \n',roi_names{vv},p1);
    end


end


