% sickle_sca_prediction.m
%
% Group level analysis of QSM (and R2*) data from the Sickle-UK cohort.
% This script is based significantly on SICKLE_GLM_ANALYSIS.m, but performs
% the inference the other way around.
%
% Data should be processed using SICKLE_QSM_PIPELINE.m and
% SICKLE_EXTRACT_AVERAGES.m and stored as a MATLAB Table in a .mat file,
% which can be read into this script.
% 
% In this version, we are trying to fit a general linear model that takes
% basal ganglia (average) susceptibility and age as predictor variables,
% to infer on either 'SCA status' or pegboard score as the response
% variable.
% 
%
%       Copyright (C) University College London, 2025
%
%
% Created by MT Cherukara, June 2025
%
% CHANGELOG:


clearvars;

%% Set-up script options

% ROIs
load('ROI_names_BG.mat');
n_rois = length(roi_names);

% Load in data table
load('SickleUK_QSMdata_BGanglia.mat');
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

% Make an 'SCA_status' logical array and add it to the table
vec_sca = strcmp(tbl_all.Group,'SS');
tbl_all = addvars(tbl_all,vec_sca,'NewVariableNames','SCA_status');


%% Loop Through ROIs and Fit the GLM

% Pre-allocate structure for storing the results
mdl = struct();

% Create an empty array for holding the following p-values of each fitting
% parameter (currently 3) 
arr_results = ones(4,n_rois);

for rr = 1:n_rois

    rname = roi_names{rr};
    vname = strcat(mname,'_',rname);

    % % Specify model and fit GLM (SCA_status)
    % modelspec = [' SCA_status ~ Log_Age:',vname,' + Sex:',vname,' + ',vname];
    % mdl.(rname) = fitglm(tbl_all,modelspec,Distribution="binomial");

    % Specify model and fit GLM (Pegboard)
    modelspec = [' Pegboard_R ~ Log_Age:',vname,' + Sex:',vname,' + ',vname];
    mdl.(rname) = fitglm(tbl_all,modelspec);

    % Store p-values
    arr_results(:,rr) = mdl.(rname).Coefficients.pValue;

end % for rr = 1:n_rois

% Convert results into a table
tbl_results = array2table(arr_results,'VariableNames',roi_names);

% Add row-names to the table
tbl_results.Properties.RowNames = mdl.(roi_names{1}).CoefficientNames;


%% Loop Through Results and Print

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
        fprintf('\n%24s ',rname);
        fprintf('adjusted R^2 = %f\t',res_RS(rr));
        fprintf('p-value = %.2e \n',res_pv(rr));
    
        if any(sigCoefs(2:end))
            sigCoefs(1) = false;
            nameCoefs = mdl.(rname).CoefficientNames(sigCoefs)';
            stats = [nameCoefs, num2cell(mdl.(rname).Coefficients.pValue(sigCoefs))]';
            fprintf('%24s has significant covariance in',rname);
            fprintf('\n%40s (%8.3g)',stats{:});
            fprintf('.\n');
        end
    % end

end % rr = 1:n_rois
