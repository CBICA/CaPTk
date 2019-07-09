/*
@ @licstart  The following is the entire license notice for the
JavaScript code in this file.

Copyright (C) 1997-2017 by Dimitri van Heesch

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

@licend  The above is the entire license notice
for the JavaScript code in this file
*/
var NAVTREE =
[
  [ "Cancer Imaging Phenomics Toolkit (CaPTk)", "index.html", [
    [ "Overview", "index.html", null ],
    [ "Installation Instructions", "Installation.html", null ],
    [ "Getting Started", "Getting_Started.html", "Getting_Started" ],
    [ "How To Guides", "How_To_Guides.html", [
      [ "Pre-processing", "How_To_Guides.html#ht_Preprocessing", [
        [ "DICOM to NIfTI conversion", "How_To_Guides.html#preprocessing_dcm2nii", null ],
        [ "Image Co-registration", "How_To_Guides.html#preprocessing_reg", null ],
        [ "Denoise-SUSAN (ITK filter)", "How_To_Guides.html#preprocessing_susan", null ],
        [ "N4 Bias Correction (ITK filter)", "How_To_Guides.html#preprocessing_biasN4", null ],
        [ "Histogram Matching", "How_To_Guides.html#preprocessing_histoMatch", null ],
        [ "Z-Scoring Normalization", "How_To_Guides.html#preprocessing_zScoreNorm", null ]
      ] ],
      [ "Segmentation", "How_To_Guides.html#ht_Segmentation", [
        [ "Geodesic Training Segmentation", "How_To_Guides.html#seg_GeoTrain", null ],
        [ "Geodesic Distance Transform-based Segmentation", "How_To_Guides.html#seg_Geodesic", null ],
        [ "ITK-SNAP", "How_To_Guides.html#seg_SNAP", null ],
        [ "Deep Learning Segmentation", "How_To_Guides.html#seg_DL", null ]
      ] ],
      [ "Feature Extraction", "How_To_Guides.html#ht_FeatureExtraction", null ],
      [ "Specialized Applications (SAs) Usage", "How_To_Guides.html#ht_SpecialApps", [
        [ "Brain Cancer", "How_To_Guides.html#Glioblastoma", [
          [ "Glioblastoma Molecular Subtype Prediction", "How_To_Guides.html#Glioblastoma_Molecular", null ],
          [ "WhiteStripe Normalization", "How_To_Guides.html#Glioblastoma_WhiteStripe", null ],
          [ "Glioblastoma EGFRvIII Surrogate Index (PHI Estimator)", "How_To_Guides.html#Glioblastoma_PHI", null ],
          [ "Glioblastoma Infiltration Index (Recurrence)", "How_To_Guides.html#Glioblastoma_Recurrence", null ],
          [ "Glioblastoma Survival Prediction Index", "How_To_Guides.html#Glioblastoma_Survival", null ],
          [ "Glioblastoma EGFRvIII SVM Index", "How_To_Guides.html#Glioblastoma_EGFRvIII", null ],
          [ "Pseudoprogression Infiltration Index", "How_To_Guides.html#Glioblastoma_Pseudoprogression", null ],
          [ "Population Atlas", "How_To_Guides.html#Glioblastoma_Atlas", null ],
          [ "Confetti", "How_To_Guides.html#Glioblastoma_Confetti", null ],
          [ "Directionality Estimator", "How_To_Guides.html#Glioblastoma_Directionality", null ]
        ] ],
        [ "Breast Cancer", "How_To_Guides.html#BreastCancer", [
          [ "Breast Density Estimation (LIBRA)", "How_To_Guides.html#BreastCancer_LIBRA", null ],
          [ "Texture Feature Extraction", "How_To_Guides.html#BreastCancer_texture", null ],
          [ "Breast Segmentation", "How_To_Guides.html#BreastCancer_breastSegmentation", null ]
        ] ],
        [ "Lung Cancer", "How_To_Guides.html#LungCancer", [
          [ "Radiomics Analysis of Lung Cancer", "How_To_Guides.html#LungCancer_SBRT", null ]
        ] ],
        [ "Miscellaneous Applications", "How_To_Guides.html#Miscellaneous", [
          [ "Perfusion Alignment", "How_To_Guides.html#misc_Perfusion_Alignment", null ],
          [ "Perfusion Derivatives", "How_To_Guides.html#misc_Perfusion_Derivatives", null ],
          [ "Diffusion Derivatives", "How_To_Guides.html#misc_Diffusion_Derivatives", null ],
          [ "PCA Volume Extraction", "How_To_Guides.html#misc_PCA_Extraction", null ],
          [ "Training Module", "How_To_Guides.html#misc_Training_Module", null ]
        ] ]
      ] ]
    ] ],
    [ "Scientific Findings using CaPTk", "Science.html", [
      [ "Non-invasive Imaging Biomarker of EGFRvIII in Glioblastoma Patients", "Science.html#phiEstimator", null ],
      [ "Prediction of Overall Survival in Glioblastoma Patients", "Science.html#survivalPredictor", null ],
      [ "Probability Maps of Potential Recurrence of Glioblastoma Tumors", "Science.html#recurrencePredictor", null ],
      [ "Imaging Biomarkers Related to Cancer Risk and Development of Breast Cancer", "Science.html#libraPapers", null ]
    ] ],
    [ "Technical Reference", "Technical_Reference.html", [
      [ "Further Application Details and Assumptions", "Technical_Reference.html#tr_Apps", [
        [ "Image Visualization", "Technical_Reference.html#appsVisualization", null ],
        [ "Extracted Features", "Technical_Reference.html#appsFeatures", null ]
      ] ],
      [ "Build CaPTk from Source", "Technical_Reference.html#tr_buildFromSource", [
        [ "Prerequisites", "Technical_Reference.html#prerequisites", null ],
        [ "Build", "Technical_Reference.html#actualBuild", null ],
        [ "Documentation & Tests", "Technical_Reference.html#optionalBuilds", null ],
        [ "Linux Build Guide", "Technical_Reference.html#linuxBuild", null ]
      ] ],
      [ "For Developers", "Technical_Reference.html#tr_forDeves", [
        [ "General Information", "Technical_Reference.html#generalInfo", null ],
        [ "Dependencies", "Technical_Reference.html#dependencies", null ],
        [ "Integrating your C++ application into CaPTk", "Technical_Reference.html#cppIntegration", null ],
        [ "Integrating your Python application into CaPTk", "Technical_Reference.html#pyIntegration", null ]
      ] ]
    ] ],
    [ "Download Instructions", "Download.html", null ],
    [ "Changelog: Release Notes", "ReleaseNotes.html", null ],
    [ "People (Credits)", "People.html", null ]
  ] ]
];

var NAVTREEINDEX =
[
"Download.html"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';