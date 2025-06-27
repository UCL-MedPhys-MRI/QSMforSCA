# QSM for the Sickle UK Dataset

Copyright Matthew Cherukara, 11 March 2025. See LICENSE.

MATLAB code accompanying the paper
- Cherukara MT, Hamdule S, Kawadler JM, Kirkham FJ, Shmueli K. Changes in Deep-Brain Magnetic Susceptibility and Cognitive Function in Children and Young People with Sickle Cell Anaemia. *SUBMITTED* (2025).

## Dependencies

- STI-Suite 3.0 [Berkeley](https://people.eecs.berkeley.edu/~chunlei.liu/software.html)  
- MEDI-Toolbox [Cornell](http://pre.weill.cornell.edu/mri/pages/qsm.html)  
- SEGUE [UCL](https://xip.uclb.com/product/SEGUE)
- MRI susceptibility calculation methods [UCL](https://xip.uclb.com/product/mri_qsm_tkd)

You will need to edit **sickle_qsm_pipeline.m** with paths to your local installations of the dependencies. Please properly credit the authors of any toolboxes used.

## How To Use

Requires multi-echo complex data, stored as NIFTIs, in BIDS format, with separate files for the magntiude and phase data. Such as:
```
dataset
└── rawdata 
    ├── sub-01
        │   ├── anat
        │   |   ├── sub-01_T1w.nii
        │   |   └── etc.
        |   └── swi
        |       ├── sub-01_part-mag_GRE.nii
        │       └── sub-01_part-phase_GRE.nii
        └── sub-02
            └── etc.
```

First, use **sickle_qsm_pipeline.m** to reconstruct QSMs using the pipeline specified in the paper. The outputs will be saved in a BIDS-formatted **derivatives/** folder within your BIDS directory. You will need to generate a brain mask for each subject (e.g. using FSL BET).

Then, use [MRIcloud](https://braingps.mricloud.org/) to generate ROI masks. You will need to upload a T1 and QSM for each subject. The file **ROI_names.mat** contains a MATLAB cell array of the names of the ROIs (with the internal and external globus pallidus combined into a single ROI).

The script **sickle_extract_averages.m** can then loop through all subjects and ROIs and extract the average susceptibility (and $R_2^*$) in each region, storing it all in a MATLAB table and CSV file (examples are provided as **SickleUK_QSMData.mat** and **SickleUK_QSMData.csv**). 

Finally, use **sickle_glm_analysis.m** to perform the general linear model (GLM) analysis of the data. Example results are provided as **SickleUK_GLMResults_QSM.mat** and **SickleUK_GLMResults_R2s.mat**. The results are stored as MATLAB tables. This script also contains cells which perform the "leave-one-out" analysis for quantifying variance explained by each variable, and cells for carrying out zero-order t-tests for (uncorrected) significance testing of two groups.

The figures in the paper can be generated using MATLAB scripts contained here:
- **figure_ROI_boxplots.m** reproduces **Fig. 2** in the paper, which shows the distribution of $\chi$ and $R_2^*$ in each ROI, comparing SCA patients with healthy controls. This does not require GLM results, and can be run after **sickle_extract_averages.m**.
- **figure_variance_explained.m** reproduces **Fig. 3** in the paper, which is a bar chart showing the proportion of inter-subject variance explained by each model parameter. This is based on leave-one-out analysis of QSM data, hence it requires those results from **sickle_glm_analysis.m**.
- **figure_scatterROI.m** and **figure_scatter2.m** reproduce the scatter graphs that are shown in **Fig. 4** in the paper. This does not require GLM results, and can be run after **sickle_extract_averages.m**, or using the data in **SickleUK_QSMData.mat**.

The file **SickleUK_SubjectData.csv** contains demographic and cognitive testing information about the subjects included in the study. These data are also contained in the first columns of **SickleUK_QSMData.csv**.

