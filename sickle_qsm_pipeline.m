% sickle_qsm_pipeline.m
%
% QSM processing pipeline for the Sickle-UK dataset, with data are stored
% in BIDS format. 
%
% For the paper, Cherukara, M.T. et al., 'Changes in Deep-Brain Magnetic
% Susceptibility and Motor Function in Children and Young People with
% Sickle Cell Anaemia', 2025 (UNDER REVIEW)
% 
%
%       Copyright (C) University College London, 2025
%
%
% Created by MT Cherukara, September 2023
%
% CHANGELOG:
%
% 2025-03-11 (MTC). Updated this version to be ready for GitHub and release
%       of the code alongside the paper.

clearvars;

% Add toolboxes to path - replace these with your own path
addpath(genpath('/home/cherukara/Documents/Coding/Toolboxes/MEDI_toolbox/'));
addpath(genpath('/home/cherukara/Documents/Coding/Toolboxes/STISuite_V3.0/'));
addpath(genpath('/home/cherukara/Documents/Coding/Toolboxes/Segue/'));
addpath(genpath('/home/cherukara/Documents/Coding/Toolboxes/AK_QSM/'));


%% SELECT OPTIONS

% Choose subject ID numbers
subs = 4:47;

% Optional saving of data
save_noise = 0;
save_multiecho = 0;

% Start the timer
tStart = tic;

% Data directory
dir_data = '/media/cherukara/DATA/Sickle_UK/SickleUK_Data_BIDS/';

% Pre-define some parameters
Params.TEs = [4.94, 12.29, 19.64, 26.99, 34.34].*1e-3;
Params.Resolution = [1, 1, 1];
Params.Orientation = [0, 0, 1];
Params.Threshold = 6e-2;
Params.Vsize = 26;          % V-SHARP maximum kernel size
Params.ErodeRad = 4;        % Brain mask erosion kernel radius
Params.PCAwin = 2;          % PCA denoising window size
Params.Alpha = 0.05;        % Iterative Tikhonov regularization parameter
Params.B0 = 1.5;

% Phase-Susceptibility Scaling Parameter
PhaScale = 2*pi*Params.B0*42.58.*(Params.TEs(2)-Params.TEs(1));




%% LOOP OVER SUBJECTS
for ss = subs

    tSub = tic;

    % Create subject name
    if ss < 10
        sname = strcat('sub-0',num2str(ss));
    else
        sname = strcat('sub-',num2str(ss));
    end

    % Path to subject directory
    dir_sub = strcat(dir_data,sname,'/');

    % Derivatives directory
    dir_deriv = strcat(dir_data,'derivatives/qsmproc/',sname,'/');

    % Create a subject specific derivatives directory
    [~, ~] = mkdir(dir_deriv);


    %% 1. Load Magnitude and Phase Data

    % Load the magnitude data
    try
        arr_mag = niftiread(strcat(dir_sub,'swi/',sname,'_part-mag_GRE.nii.gz'));
        inf_mag = niftiinfo(strcat(dir_sub,'swi/',sname,'_part-mag_GRE.nii.gz'));
    catch
        arr_mag = niftiread(strcat(dir_sub,'swi/',sname,'_part-mag_GRE.nii'));
        inf_mag = niftiinfo(strcat(dir_sub,'swi/',sname,'_part-mag_GRE.nii'));
    end

    % Load the phase data
    try
        arr_pha = niftiread(strcat(dir_sub,'swi/',sname,'_part-phase_GRE.nii.gz'));
        inf_pha = niftiinfo(strcat(dir_sub,'swi/',sname,'_part-phase_GRE.nii.gz'));
    catch
        arr_pha = niftiread(strcat(dir_sub,'swi/',sname,'_part-phase_GRE.nii'));
        inf_pha = niftiinfo(strcat(dir_sub,'swi/',sname,'_part-phase_GRE.nii'));
    end

    % Convert to single precision
    arr_mag = single(arr_mag);
    arr_pha = single(arr_pha);

    % Remove NaNs
    arr_mag(isnan(arr_mag)) = 0;
    arr_pha(isnan(arr_pha)) = 0;
    
    % Scale the phase from 0 to 2*pi
    arr_pha = arr_pha - min(arr_pha(:));
    arr_pha = 2*pi*arr_pha./max(arr_pha(:));

    % Load the BET-calculated mask
    arr_mask = niftiread(strcat(dir_sub,'anat/',sname,'_desc-brain_mask.nii'));
    inf_mask = niftiinfo(strcat(dir_sub,'anat/',sname,'_desc-brain_mask.nii'));

    % Convert to LOGICAL, for now
    arr_mask = logical(arr_mask);

    
    %% 2. Non-Linear Fitting of Complex Data


    % Make complex data
    arr_comp = arr_mag .* exp(-1i*arr_pha);

    % iField correction using MEDI function
    arr_comp = iField_correction_new(arr_comp,Params.Resolution,arr_mask);
    fprintf('Completed iField Correction \n');


    % Apply MP-PCA denoising using Vertaart and Does script
    [arr_comp,arr_var,arr_P] = denoiseCV(arr_comp,Params.PCAwin.*[1,1,1]);
    fprintf('Completed MP-PCA Denoising \n');

    % Save out denoised multi-echo data
    if save_multiecho == 1
        niftiwrite(angle(arr_comp), strcat(dir_sub,'swi/',sname,'_part-phase_denoise-MPPCA_GRE.nii'), inf_pha);
        niftiwrite(abs(arr_comp), strcat(dir_sub,'swi/',sname,'_part-mag_denoise-MPPCA_GRE.nii'), inf_mag);
    end

    % Calculate R2* using ARLO
    arr_r2s = arlo(Params.TEs,abs(arr_comp));

    % Create and update an INFO struct for R2* data
    inf_r2s = inf_mag;
    inf_r2s.ImageSize = size(arr_r2s);
    inf_r2s.PixelDimensions = Params.Resolution;

    % Save the R2* map
    niftiwrite(single(arr_r2s), strcat(dir_sub,'anat/',sname,'_acq-ARLO_R2smap.nii'), inf_r2s);
    fprintf('Completed R2* Mapping \n');

    % Non-linear fitting using MEDI function
    [arr_linfit, arr_noise] = Fit_ppm_complex_bipolar(arr_comp);

    % Make noise map absolute
    arr_noise = abs(arr_noise);

    % Create and update an INFO struct for the fieldmap data
    inf_linfit = inf_mag;
    inf_linfit.ImageSize = size(arr_linfit);
    inf_linfit.PixelDimensions = Params.Resolution;

    % Save out the non-linear fit data
    niftiwrite(single(arr_linfit), strcat(dir_deriv,sname,'_acq-MEDI_linfit.nii'), inf_linfit);

    % Save the noise map
    if save_noise == 1

        % Create and update an INFO struct for the noise data
        inf_noise = inf_mag;
        inf_noise.ImageSize = size(arr_noise);
        inf_noise.PixelDimensions = Params.Resolution;

        % Save out the noise data
        niftiwrite(single(arr_noise), strcat(dir_deriv,sname,'_acq-MEDI_noise.nii'), inf_noise);

    end

    % Update stage
    fprintf('Completed Non-Linear Fitting \n');

      

    %% 3. Phase Unwrapping


    % SEGUE unwrapping as implemented in AK's MATLAB toolbox
    SParams.Phase = double(arr_linfit);
    SParams.Mask = double(arr_mask);
    arr_unwrapped = Segue(SParams);

    % Save the unwrapped phase
    inf_unwrapped = inf_pha;
    inf_unwrapped.ImageSize = size(arr_unwrapped);
    inf_unwrapped.PixelDimensions = Params.Resolution;
    niftiwrite(single(arr_unwrapped), ...
        strcat(dir_sub,'swi/',sname,'_part-phase_unwrapped-SEGUE_GRE.nii'),...
        inf_unwrapped);

    % Update
    fprintf('Completed Phase Unwrapping \n');
    

    %% 4. Background Field Removal

    % First, erode the mask - we are always going to do this
    arr_mask = imerode(arr_mask,strel('disk',Params.ErodeRad));

    % Save out the mask
    if save_mask == 1
        niftiwrite( int16(arr_mask), strcat(dir_sub,'anat/',sname,'_method-',method_ms,'_mask.nii'), inf_mask);
    end

    % Perform V-SHARP background field removal using Liu's script
    arr_fieldmap = V_SHARP(arr_unwrapped,arr_mask,'voxelsize',Params.Resolution,'smvsize',Params.Vsize);

    % Create and update INFO struct for the local fieldmap
    inf_fieldmap = inf_linfit;

    % Save the local fieldmap
    niftiwrite( single(arr_fieldmap),...
                strcat(dir_sub,'swi/',sname,'_unwrapped-SEGUE',...
                '_bfr-VSHARP_fmap.nii'), inf_fieldmap);

    % Update
    fprintf('Completed Background Field Removal \n');
    

    %% 5. Susceptibility Calculation

    % Set some parameters
    Params.MatrixSize = size(arr_fieldmap);

    % Call AK's IterTik function
    arr_susc = ak_Tikhonov_iter(arr_fieldmap.*arr_mask, arr_noise.*arr_mask, Params);
    arr_susc = arr_susc./PhaScale;


    % Create and update INFO struct for the susceptibility map
    inf_susc = inf_mag;
    inf_susc.ImageSize = size(arr_susc);
    inf_susc.PixelDimensions = Params.Resolution;

    % Save the QSM using NIFTIWRITE
    niftiwrite( arr_susc, strcat(dir_sub,'swi/',sname,'_unwrapped-SEGUE',...
                '_bfr-VSHARP','_susc-iterTik_susc.nii'), inf_susc);
    

    % Timings
    fprintf('Completed Susceptibility Calculation \n');
   
    tSubend = toc(tSub);
    fprintf('Completed subject %d \n\t Total time elapsed: %.4f \n', ss,tSubend);


end % for ss = subs

tEnd = toc(tStart);
fprintf('Total time elapsed: %.4f \n', tEnd);
