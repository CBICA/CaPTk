..
  Computational Breast Imaging Group
  Department of Radiology
  University of Pennsylvania
  Richards Building
  3700 Hamilton Walk, 7th floor
  Philadelphia, PA 19104

  Web:   http://www.uphs.upenn.edu/radiology/research/labs/cbig/

  Copyright (c) 2014-2016 University of Pennsylvania. All rights reserved.
  See http://www.cbica.upenn.edu/sbia/software/license.html or COPYING file.


  Introduction
  ============

  The amount of fibroglandular tissue content in the breast as estimated mammographically, commonly referred to as breast percent density (PD%), is one of the most significant risk factors for developing breast cancer. Approaches to quantify breast density commonly focus on either semiautomated methods or visual assessment, both of which are highly subjective. This software package was developed to be a fully-automated density estimation method that works on both raw (i.e., "FOR PROCESSING") and vendor postprocessed (i.e., "FOR PRESENTATION") digital mammography images,and has thus far been validated to work on GE Healthcare and Hologic digital mammography systems.

  Briefly, the software first applies an edge-detection algorithm to delineate the boundary of the breast and the boundary of the pectoral muscle. Following the segmentation of the breast, an adaptive multi-class fuzzy c-means algorithm is applied to identify and partition the mammographic breast tissue area, into multiple regions (i.e., clusters) of similar x-ray attenuation. These clusters are then aggregated by a support-vector machine classifier to a final dense tissue area, segmentation. The ratio of the segmented absolute dense area to the total breast area is then used to obtain a measure of breast percent density (PD%). 

  The software saves both the physical segmentations of the breast and dense tissues (as .mat files with binary image matrices), as well as the quantitative estimates of breast area, dense area and PD% to a comma seperated text file (.csv, openable by Excel) to a user-defined directory, as described below.

  LIBRA has two modes of operations: 

  1. An easy-to-use interactive mode with Graphical-User-Interface where the user is prompted to select either a single DICOM image or a folder of DICOM images, an output folder for the results, and whether they wish to save intermediate files. 
  2. A command-line interface amenable to batch processing and scripting, where the user can explicitly define the input and output paths.
  
  Requirements
  ============

  This package can be run and is fully tested on both Windows and Unix platforms. It has been widely tested and validated with MATLAB R2012a, R2012b, R2013a, R2014a, and R2014b. CAUTION: There is a known bug and incompatibility to MATLAB R2013b (version: 8.2.0.705) on Windows that the dicomread.m function would randomly crash and terminate the matlab routine.

  Several MATLAB toolboxes are required to run the software. Please make sure you have valid licenses for all of the following:
	  - Curve Fitting Toolbox
	  - Image Processing Toolbox
	  - Signal Processing Toolbox
	  - Statistics Toolbox
	  - Symbolic Math Toolbox
	  - MATLAB Compiler (optional for compiling source into an executable)
..

Source Package Overview
=======================

The package contains::

  Source/
  |--->	Code/			# Main source codes
  |	|-> libra.m
  |	|-> libra_exper.m
  |	|-> libra_run.m
  |	|-> libra_version.m
  |	`-> private/
  |
  |--->	Model/			# SVM models used by libra.m
  |	|-> GE_Processed_SVM.mat
  |	|-> GE_Raw_SVM.mat
  |	|-> Hologic_Processed_SVM.mat
  |	`-> Hologic_Raw_SVM.mat
  |
  |--->	libra_demo.m		# Demo/testing program
  |--->	libra_startup.m		# Libra start up command
  `--->	libra_compile.m		# Function to comiple source code into an executable

  Sample_Data/			# Data used in libra_demo.m
    |-> Case1.dcm
    |-> Case2.dcm
    |-> Case3.dcm
    |-> Case4.dcm
    |-> Case5.dcm
    `-> ground_truth.mat
    
    Executable/
      `-> libra(.exe)       # Pre-compiled standalone executable 
                            # (Windows platform extension in parentheses.)

  AUTHORS.txt         		# A list of the people who contributed to this software.
  COPYING.txt         		# The copyright and license notices.
  ChangeLog.txt       		# History and evolution of the software
  README.txt          		# This readme file.
  LIBRA_Software_Manual.pdf	# Software manual

  doc/                    # Documentation files for the software pacakge (For Developer) 

The Source/Code/ and Source/Code/private/ directories include all the core functions written in MATLAB. The models in Source/Model/ that are vendor-specific are GE_Processed_SVM.mat, GE_Raw_SVM.mat, Hologic_Processed_SVM.mat and Hologic_Raw_SVM.mat. Note: If a mammogram being analyzed is derived from a vendor for which no training currently exists (or the vendor is unknown), the system defaults to using Hologic-trained models. If the presentation intent type (i.e., 'For Processing' or 'For Presentation') is not defined in the DICOM header, the system defaults to using 'For Presentation' training. It is also worth noting that although training has only currently been completed for GE Healthcare and Hologic systems, Siemens 'For Processing' mammograms has also been evaluated with the system and have been found to have an acceptable level of performance.

Source/Code/libra.m is a wrapper script that performs density estimation 1) directly on one given mammogram or 2) loops through (NOT recursively) a given input directory to find out all the image to be processed. libra will output the result in a given output directory with a Density.csv and a sub-directory Result_Images/ containing final segmentation and breast mask (optional.) There are two modes: interactive mode with GUI and batch processing mode at command-line.

Source/Code/libra_exper.m is the main code that executes the breast density estimation pipeline given one input DICOM image and one output directory and the output text file.

Sample_Data/ contains sample images and testing ground truth for the sample script libra_demo.m.

doc/ contains Sphinx compatible documentation rst files that can be used to automatically generate website and pdf documentation. It is reserved for developers.

Usage
=====

**The following sections are served as a guideline for users intending to use the LIBRA package from the Source code. Please be advised to follow the instruction and execute the commands in MATLAB environment. For the instruciton on Graphical-User-Interface mode, please check out the section "Running the Software - Interactive Mode
" in LIBRA_Software_Manual.pdf.**

In Matlab environment, run::

  >> libra_startup

to initialize the environment for the software.

To run breast density estimation, run::

  >> RESULTS = libra(PATH_INPUT, PATH_OUTPUT_DIR, SAVE_INTERMEDIATES);

where PATH_INPUT should be a mammogram with dcm as file extension or a directory containing dicom images (*.dcm) to be processed. PATH_OUTPUT_DIR is a path to the output directory. For example,

in Windows::

  >> RESULTS = libra('D:\LIBRA\Sample_Data\Case1.dcm', 'D:\LIBRA\Result\', 1);

or in Linux::

  >> RESULTS = libra('/path/to/LIBRA/Sample_Data/Case1.dcm', '/path/to/LIBRA/Result/', 1);

Program returns a cell array RESULTS that has Filename, BreastArea(sqcm), DenseArea(sqcm), and BreastDensity(%) as four columns, as below::

  RESULTS =
  
      'File Analyzed'  'BreastArea(sqcm)'    'DenseArea(sqcm)'    'BreastDensity(%)'
      'Case3'          [        142.1792]    [        31.4047]    [         22.0881]

You can expect the outputs to be::

	PATH_OUTPUT_DIR/
		|-> Result_Images/
		|  	|-> Masks_<PREFIX>.mat
		|  	|-> <PREFIX>_density_segmentation.jpg
		|  	|-> <PREFIX>_Windowed_Original.jpg (optional)
		|  	`-> <PREFIX>_density_imagesc.jpg (optional)
		|-> Density.csv
		`-> SkippedImg.csv (if any failed)

Optional outputs can be saved by setting a boolean variable SAVE_INTERMEDIATES to TRUE (or 1.) By default, these optional files are not saved.

Output segmentation images will be generated in PATH_OUTPUT_DIR/Result_Images/ while the estimated desity number will be written in comma seperated file, Density.csv, and failure information, if any, due to lack of required imaging information in the dicom header to be used in the classification will be written into SkippedImg.csv in PATH_OUTPUT_DIR. You can find the breast mask and the estimated values in Masks_*.mat file.

By doing::

  >> help libra

you will get a full description and usage of the program.

Testing
=======

To test the package with the sample data to see if the package runs correctly on your machine, please change the directory to the Source/ directory, simply run the following commands in MATLAB environment::

  >> libra_startup
  >> libra_demo

Output directories "Demo_Test" and "Demo_Test_Batch" will be created inside the source package. For more information on the output results please refer to the USAGE section. 

The demo not only demonstrates breast and dense tissue segmentations but also tests whether your computing environment produces comparable results against precomputed "ground truth". You may find the ground truth in Sample_Data/ground_truth.mat and manually compare the numbers with the results in Demo_Test/Density.csv::

  >> load ground_truth.mat
  >> whos
    Name     Size   Bytes  Class     Attributes
    truth    4x3       96  double

  >> truth

  truth =

     70.6779   52.8680   74.7996      % Ground truth for Case1.dcm
    165.6580   12.8664    7.7645      % Case2.dcm
    143.0126   32.2381   22.5421      % Case3.dcm
     92.6134   21.0520   22.7311      % Case4.dcm
     80.3807   26.1074   32.4397      % Case5.dcm

Each row contains three numbers. They are BreastArea(sqcm), DenseArea(sqcm), and BreastDensity(%), respectively.

Compilation
===========

To compile the libra.m into a executable, change the directory to the Source/ directory within the main directory of the LIBRA package and run the following command in MATLAB environment::

  >> libra_compile('D:\path\to\output\executable\folder\') % In Windows
  
or::

  >> libra_compile('/path/to/output/executable/directory/') % In Unix

A stand-alone executable libra (libra.exe in Windows) will be created in the specified directory. When libra_compile is called without an input argument, an "Executable" directory will be created in Source directory and "libra" executable will be generated there. This compilation requires MATLAB Compiler (mcc).

After compilation, you can use the stand-alone executable without matlab environment. [*]_
::

  $ libra PATH_INPUT_DIR PATH_OUTPUT_DIR

.. [*] Note: Matlab Compiler Runtime (MCR) is required. Download a version of MCR that corresponds to the MCC used to compile the executable from http://www.mathworks.com/products/compiler/mcr/.

..
  References
  ==========
  When publishing results or data derived from this software, please include a bibliographical reference to the 2012 Keller et al. Medical Physics paper (PubMed ID: 22894417) describing the methodology as well as to the 2015 Keller et al. Breast Cancer Research paper (PubMed ID: 26303303) evaluating the association between LIBRA breast density and breast cancer risk in a case-control study.

  Keller, B.M., Nathan, D.L., Wang, Y., Zheng, Y., Gee, J.C., Conant, E.F., and Kontos, D., "Estimation of breast percent density in raw and processed full field digital mammography images via adaptive fuzzy c-means clustering and support vector machine segmentation," Medical physics 39 (8), 4903-4917.

  B.M. Keller, J. Chen, D. Daye, E.F. Conant, and D. Kontos. "Preliminary evaluation of the publicly available Laboratory for Breast Radiodensity Assessment (LIBRA) software tool: comparison of fully automated area and volumetric density measures in a case–control study with digital mammography." Breast Cancer Research 17(1), 1-17
..
