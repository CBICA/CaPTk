@echo off
rem This is a script for executing the default mode for the HGG data of the Brats18 challenge
rem Run this batch file for every subject of the dataset
rem You need to provide the dataset name as an argument (for example: Brats18_CBICA_AAB_1)
set dataset=%1

rem -------------------------------
rem Change these values accordingly
rem -------------------------------
set exe_path="C:/GeodesicTraining/build/Debug/GeodesicTraining.exe"
set dataset_folder="C:/GeodesicTraining/datasets/MICCAI_BraTS_2018_Data_Training_HGG_only"
set outputDir="C:/GeodesicTraining/output/Brats18"
set input_flair="%dataset_folder%/%dataset%/%dataset%_flair.nii.gz"
set input_t1="%dataset_folder%/%dataset%/%dataset%_t1.nii.gz"
set input_t1ce="%dataset_folder%/%dataset%/%dataset%_t1ce.nii.gz"
set input_t2="%dataset_folder%/%dataset%/%dataset%_t2.nii.gz"
set labels="%dataset_folder%/%dataset%/mask.nii.gz"
set ground_truth="%dataset_folder%/%dataset%/%dataset%_seg.nii.gz" 
set gtmode=reversegeotrain
set extra_run=--reportseconds

rem ---------------------------
rem Don't change anything below
rem ---------------------------
set str_run_def=%exe_path% -m %gtmode% -d %dataset% -o %outputDir% -l %labels%
set mri_images_run=-flair %input_flair% -t1 %input_t1% -t1ce %input_t1ce% -t2 %input_t2%
set gt_run=-g %ground_truth% -gs 0
set default_run=%str_run_def% %mri_images_run% %gt_run% %extra_run%

rem Run 
echo. & echo %dataset% & echo. & %default_run%
