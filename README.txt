/**

\mainpage Overview

CaPTk is a software platform, written in C++, for analysis of radiographic images of cancer, currently focusing on brain, breast, and lung cancer. CaPTk integrates advanced, validated tools performing various aspects of medical image analysis that have been developed in the context of active clinical research studies and collaborations toward addressing real clinical needs. With emphasis given on being a very lightweight and efficient image viewer and eliminating the prerequisite for a substantial computational background, CaPTk aims to facilitate the swift translation of advanced computational algorithms into routine clinical quantification, analysis, decision making, and reporting workflow.

Its long-term goal is to provide widely used technology that makes use of advanced imaging analytics in cancer prediction, diagnosis and prognosis, as well as in better understanding the biological mechanisms of cancer development.

The package leverages the value of quantitative imaging analytics along with machine learning to derive phenotypic imaging signatures, based on two-level functionality. At the first level, image analysis algorithms are used to extract a rich panel of diverse and complementary features, such as multi-parametric intensity histograms, textural, morphologic, and kinetic variables, connectomics, and spatial patterns. At the second level, these radiomic features are fed into multi-variable machine learning models to produce diagnostic,  prognostic and predictive  biomarkers (specific examples are given in the “<a href="Science.html">Scientific Findings</a>” section). Fig. 1 shows results from clinical studies in three areas:

-# Computational neuro-oncology of brain gliomas for precision diagnostics, prediction of outcome and subsequent treatment planning.
-# Prediction of treatment response for breast and lung cancer.
-# Risk assessment for breast cancer.

CaPTk is developed and maintained by the Center for Biomedical Image Computing and Analytics (CBICA - https://www.cbica.upenn.edu) at the University of Pennsylvania, and draws upon research from several groups within the Center and beyond. 

\image html 0_overview_resize.png "Fig.1. Overview of all functions and applications of CaPTk, in its two-level architecture"
\image latex 0_overview_resize.png "Fig.1. Overview of all functions and applications of CaPTk, in its two-level architecture"

## Bug Tracker and Feature Request
 
We coordinate our bugs and feature requests via out GitHub page: https://github.com/CBICA/CaPTk/issues

## Frequently Asked Questions (FAQ)

Please see our [FAQ Section](gs_FAQ.html).

## Supporting Grant
This work is supported by the NIH/NCI/ITCR* grant U24-CA189523.
<br>* National Institutes of Health / National Cancer Institute / Informatics Technology for Cancer Research

## Disclaimer
- The software has been designed for research purposes only and has neither been reviewed nor approved for clinical use by the Food and Drug Administration (FDA) or by any other federal/state agency.
- This code (excluding dependent libraries) is governed by the license provided in https://www.med.upenn.edu/sbia/software-agreement.html unless otherwise specified.
- The minimum recommended resolution is 1024x768. We have seen some issues with high DPI screens and bug reports related to it will be appreciated.

## Contact
For more information, please contact <a href="mailto:software@cbica.upenn.edu">software@cbica.upenn.edu</a>.

--------------------------------------------------------------------

--------------------------------------------------------------------

\htmlonly
<div align="right"><a href="Getting_Started.html"><b>Next (Getting Started)<b></a>
\endhtmlonly

--------------------------------------------------------------------
*/
