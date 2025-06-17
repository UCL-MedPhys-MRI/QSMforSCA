% sickle_extract_averages.m
%
% Uses MRIcloud derived segmentations of deep-brain grey matter regions,
% and QSMs created by SICKLE_QSM_PIPELINE.m to extract regional average
% values of susceptibility in the Sickle-UK data.
%
% The QSMs and the MRIcloud dseg maps must be stored in BIDS format. You
% also need the file "SickleUK_SubjectData.csv"
% 
%
%       Copyright (C) University College London, 2025
%
%
% Created by MT Cherukara, November 2023
%
% CHANGELOG:
%
% 2024-01-19 (MTC). Added a section that averages out the two globus
%       pallidus regions into a single ROI
%
% 2025-06-17 (MTC). Refactored and updated to be ready for release with the
%       paper

clearvars;


%% Set Up Script Options

% Data directory
dir_data = '/media/cherukara/DATA/Sickle_UK/SickleUK_Data_BIDS/';

% Load a preliminary table which has the ROI names
tbl_names = readtable(strcat(dir_data,'sub-02/swi/sub-02_desc-MRIcloud_stats.txt'));

% extract ROI names (omit the last 4) 
roi_names = tbl_names.labelname(1:21);

% Load the subject demographic data
%       Note that subjects 17 and 34 are missing
tbl_subs = readtable('./SickleUK_SubjectData.csv');
tbl_subs.Properties.RowNames = tbl_subs.Number;

% Numbers
n_subs = numel(tbl_subs.Number);
n_rois = numel(roi_names);


%% Calculate Averages

% Pre-allocate arrays of ROI average data
arr_qsm_mean = zeros(n_subs,n_rois);
arr_r2s_mean = zeros(n_subs,n_rois);
arr_roi_size = zeros(n_subs,n_rois);

% Loop through subjects
for ss = 1:n_subs

    sname = tbl_subs.Number{ss};
    
    % Load the data for this subject
    mat_qsm = niftiread(strcat(dir_data,sname,'/swi/',sname,'_reg-T1_susc.nii.gz'));
    mat_r2s  = niftiread(strcat(dir_data,sname,'/anat/',sname,'_reg-T1_R2smap.nii.gz'));
    mat_seg = single(niftiread(strcat(dir_data,sname,'/anat/',sname,'_desc-MRIcloud_dseg.nii')));

    % Loop through ROIs and calculate mean values
    for rr = 1:n_rois

        % Linearize the data
        vec_mask = mat_seg(:) == rr;
        vec_qsm = mat_qsm(:);
        vec_r2s = mat_r2s(:);

        % Apply the mask
        vec_qsm(vec_mask == 0) = [];
        vec_r2s(vec_mask == 0) = [];

        % Calculate and store averages
        arr_qsm_mean(ss,rr) = mean(vec_qsm,'omitnan');
        arr_r2s_mean(ss,rr) = mean(vec_r2s,'omitnan');
        arr_roi_size(ss,rr) = nnz(vec_mask);
        % arr_features(ss,rr) = mtc_extractfeatures(mat_qsm(:),vec_mask);


    end % for rr = 1:n_rois

end % for ss = subs


%% Create new ROI names

% Pre-allocate ROI name arrays
roi_names_qsm = cell(size(roi_names));
roi_names_r2s = cell(size(roi_names));

% Generate table heading names
for rr = 1:n_rois

    roi_names_qsm{rr} = strcat('QSM_',roi_names{rr});
    roi_names_r2s{rr} = strcat('R2s_',roi_names{rr});

end % for rr = 1:n_rois



%% Store the data in tables

% Convert average arrays into tables
tbl_qsmdata = array2table(arr_qsm_mean,'VariableNames',roi_names_qsm,'RowNames',tbl_subs.Number);
tbl_r2sdata = array2table(arr_r2s_mean,'VariableNames',roi_names_r2s,'RowNames',tbl_subs.Number);


%% Combine Globus Pallidus Regions

% Pull out relevant data
gpl_sz = arr_roi_size(:,[3,5]);
gpr_sz = arr_roi_size(:,[4,6]);
gpl_qsm = arr_qsm_mean(:,[3,5]);
gpr_qsm = arr_qsm_mean(:,[4,6]);
gpl_r2s = arr_r2s_mean(:,[3,5]);
gpr_r2s = arr_r2s_mean(:,[4,6]);

% Calculate new mean values
vec_gp_l_qsm = sum(gpl_sz.*gpl_qsm,2)./sum(gpl_sz,2);
vec_gp_r_qsm = sum(gpr_sz.*gpr_qsm,2)./sum(gpr_sz,2);
vec_gp_l_r2s = sum(gpl_sz.*gpl_r2s,2)./sum(gpl_sz,2);
vec_gp_r_r2s = sum(gpr_sz.*gpr_r2s,2)./sum(gpr_sz,2);

% Add data to the tables
tbl_qsmdata = addvars(tbl_qsmdata,vec_gp_l_qsm,vec_gp_r_qsm,'NewVariableNames',{'QSM_Globus_pallidus_L';'QSM_Globus_pallidus_R'});
tbl_r2sdata = addvars(tbl_r2sdata,vec_gp_l_r2s,vec_gp_r_r2s,'NewVariableNames',{'R2s_Globus_pallidus_L';'R2s_Globus_pallidus_R'});


%% Average Across the Basal Ganglia
% Following Carpenter, 2016, we define basal ganglia as Caudate, Globus
% Pallidus (int. + ext.), Putamen, Pulvinar, Substantia Nigra

% List of ROI indices
ind_bgl = [1,3,5,7,11,15];
ind_bgr = [2,4,6,8,12,16];
ind_bgt = [ind_bgl, ind_bgr];

% Pull out sizes
bgl_sz = arr_roi_size(:,ind_bgl);
bgr_sz = arr_roi_size(:,ind_bgr);
bgt_sz = arr_roi_size(:,ind_bgt);

% Pull out QSM and R2s data
bgl_qsm = arr_qsm_mean(:,ind_bgl);
bgr_qsm = arr_qsm_mean(:,ind_bgr);
bgt_qsm = arr_qsm_mean(:,ind_bgt);
bgl_r2s = arr_r2s_mean(:,ind_bgl);
bgr_r2s = arr_r2s_mean(:,ind_bgr);
bgt_r2s = arr_r2s_mean(:,ind_bgt);

% Calculate new mean values
vec_bg_l_qsm = sum(bgl_sz.*bgl_qsm,2)./sum(bgl_sz,2);
vec_bg_r_qsm = sum(bgr_sz.*bgr_qsm,2)./sum(bgr_sz,2);
vec_bg_t_qsm = sum(bgt_sz.*bgt_qsm,2)./sum(bgt_sz,2);
vec_bg_l_r2s = sum(bgl_sz.*bgl_r2s,2)./sum(bgl_sz,2);
vec_bg_r_r2s = sum(bgr_sz.*bgr_r2s,2)./sum(bgr_sz,2);
vec_bg_t_r2s = sum(bgt_sz.*bgt_r2s,2)./sum(bgt_sz,2);

% Add data to the tables
tbl_qsmdata = addvars(tbl_qsmdata,vec_bg_l_qsm,vec_bg_r_qsm,vec_bg_t_qsm,...
                      'NewVariableNames',{'QSM_Basal_ganglia_L';'QSM_Basal_ganglia_R';'QSM_Basal_ganglia_T'});
tbl_r2sdata = addvars(tbl_r2sdata,vec_bg_l_r2s,vec_bg_r_r2s,vec_bg_t_r2s,...
                      'NewVariableNames',{'R2s_Basal_ganglia_L';'R2s_Basal_ganglia_R';'R2s_Basal_ganglia_T'});


%% Finally Assembly of the Data

% Merge the tables
tbl_all = join(tbl_subs,tbl_qsmdata,'Keys','Row');
tbl_all = join(tbl_all,tbl_r2sdata,'Keys','Row');

% Save the table
writetable(tbl_all,strcat(dir_data,'SickleUK_QSMData_BGanglia.csv'),'WriteRowNames',true);
save('SickleUK_QSMdata_BGanglia.mat','tbl_all');